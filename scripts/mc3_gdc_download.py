"""
GDC MC3 publication supplement downloader (NCI MC3 / Pan-Cancer Atlas).

Uses official GDC manifest files (UUID, filename, md5, size) from:
  https://gdc.cancer.gov/about-data/publications/mc3-2017

Downloads via: https://api.gdc.cancer.gov/data/<file_id>

Open-access files work without credentials. Controlled-access files require a GDC
download token (same as GDC Data Transfer Tool). Set environment variable
GDC_TOKEN or pass --gdc-token.

Note: `vep_supporting_files.tar.gz` is listed on the publication page but does not
appear in PanCan-MC3_Open_GDC-Manifest_1.txt or PanCan-MC3_Controlled_GDC-Manifest_0.txt;
it is skipped here with a warning—fetch manually if needed.

If you meant ACC.tar.gz (TCGA adrenocortical), not "CC.tar.gz", the target list uses ACC.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
from pathlib import Path

import requests

# Official manifest URLs (tab-separated: id, filename, md5, size)
GDC_MC3_OPEN_MANIFEST_URL = (
    "https://gdc.cancer.gov/system/files/public/file/PanCan-MC3_Open_GDC-Manifest_1.txt"
)
GDC_MC3_CONTROLLED_MANIFEST_URL = (
    "https://gdc.cancer.gov/system/files/public/file/PanCan-MC3_Controlled_GDC-Manifest_0.txt"
)
GDC_DATA_ENDPOINT = "https://api.gdc.cancer.gov/data"

# Filenames requested for the class pipeline (ACC = TCGA cohort; not "CC.tar.gz").
MC3_TARGET_FILENAMES: tuple[str, ...] = (
    "mc3.v0.2.8.PUBLIC.maf.gz",
    "ACC.tar.gz",
    "BLCA.tar.gz",
    "BLCA2.tar.gz",
    "BRCA.tar.gz",
    "CESC.tar.gz",
    "CHOL.tar.gz",
    "COAD.tar.gz",
    "DLBC.tar.gz",
    "ESCA.tar.gz",
    "GBM.tar.gz",
    "HNSC.tar.gz",
    "KICH.tar.gz",
    "KIRC.tar.gz",
    "KIRP.tar.gz",
    "LAML.tar.gz",
    "LGG.tar.gz",
    "LIHC.tar.gz",
    "LUAD.tar.gz",
    "LUSC.tar.gz",
    "MESO.tar.gz",
    "OV.tar.gz",
    "PAAD.tar.gz",
    "PCPG.tar.gz",
    "PRAD.tar.gz",
    "READ.tar.gz",
    "SARC.tar.gz",
    "SKCM.tar.gz",
    "STAD.tar.gz",
    "TGCT.tar.gz",
    "THCA.tar.gz",
    "THYM.tar.gz",
    "UCEC.tar.gz",
    "UCS.tar.gz",
    "UVM.tar.gz",
    "broad_mc3_vcfs.tar.gz",
    "ref_data_for_oxog_all_permissions.tar.gz",
    "ref_data_for_oxog.tar",
    "gencode.v19.basic.exome.bed",
    "gaf_20111020Plusbroad_wex_1.1_hg19.bed",
    "contestkeys.txt",
    "mark_bitgt.txt",
    "mark_nonexonic.txt",
    "ndp.mark.txt",
    "nonpreferredpair_maf_keys.txt",
    "oxog.annotation",
    "strandBias.filter_v2.txt.gz",
    "wga_maf_keys.txt",
    "pcadontusekeys.txt",
    "pancan.merged.v0.2.exac_pon_tagged.txt.gz",
    "pancan.merged.v0.2.pfiltered.broad_pon_tagged_v2.txt.gz",
    "readme.txt",
    "mc3_qc_scores_2016-04-26a.txt",
    "final_summed_tokens.hist.bin",
    "oxoG_docker_with_patch_and_minibam.tar",
    "vep_supporting_files.tar.gz",
    "pairing_of_best_bams_2016-04-05.txt",
)


def _parse_manifest_text(text: str) -> dict[str, tuple[str, int, str]]:
    """Map filename -> (file_id, size_bytes, md5_hex)."""
    out: dict[str, tuple[str, int, str]] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            parts = line.split()
        if len(parts) < 4:
            continue
        file_id, filename, md5_hex, size_s = parts[0], parts[1], parts[2], parts[3]
        # Header row: id filename md5 size
        if file_id.lower() == "id":
            continue
        out[filename] = (file_id, int(size_s), md5_hex.lower())
    return out


def _fetch_manifest(session: requests.Session, url: str) -> dict[str, tuple[str, int, str]]:
    r = session.get(url, timeout=120)
    r.raise_for_status()
    return _parse_manifest_text(r.text)


def _open_file_id_set(open_map: dict[str, tuple[str, int, str]]) -> set[str]:
    return {t[0] for t in open_map.values()}


def _download_one(
    session: requests.Session,
    file_id: str,
    dest: Path,
    expected_size: int,
    expected_md5: str,
    token: str | None,
    verify_md5: bool,
) -> None:
    url = f"{GDC_DATA_ENDPOINT}/{file_id}"
    headers = {}
    if token:
        headers["X-Auth-Token"] = token

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")

    with session.get(url, headers=headers, stream=True, timeout=600) as r:
        if r.status_code == 403:
            raise PermissionError(
                f"GDC returned 403 for {file_id}. Controlled-access files need a valid "
                "GDC token (export GDC_TOKEN=... from https://portal.gdc.cancer.gov )."
            )
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=8 * 1024 * 1024):
                if chunk:
                    f.write(chunk)

    got = tmp.stat().st_size
    if got != expected_size:
        tmp.unlink(missing_ok=True)
        raise OSError(f"Size mismatch for {dest.name}: got {got}, expected {expected_size}")

    if verify_md5:
        h = hashlib.md5()
        with open(tmp, "rb") as f:
            for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
                h.update(chunk)
        if h.hexdigest() != expected_md5:
            tmp.unlink(missing_ok=True)
            raise OSError(f"MD5 mismatch for {dest.name}")

    tmp.replace(dest)


def download_mc3_publication_supplements(
    out_dir: str | Path,
    *,
    gdc_token: str | None = None,
    skip_controlled: bool = False,
    verify_md5: bool = False,
    overwrite: bool = False,
    session: requests.Session | None = None,
) -> None:
    """
    Download MC3 publication files listed in MC3_TARGET_FILENAMES into out_dir.

    Parameters
    ----------
    out_dir
        Directory under which files are stored (flat; e.g. data/raw/mc3).
    gdc_token
        GDC download token for controlled-access files.
    skip_controlled
        If True, only download open-access manifest files that appear in the target list.
    verify_md5
        If True, verify MD5 after download (slow for multi-GB files).
    overwrite
        If False, skip download when dest exists and size matches expected.
    """
    out_path = Path(out_dir).expanduser().resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    own_session = session is None
    sess = session or requests.Session()

    try:
        try:
            open_map = _fetch_manifest(sess, GDC_MC3_OPEN_MANIFEST_URL)
            controlled_map = _fetch_manifest(sess, GDC_MC3_CONTROLLED_MANIFEST_URL)
        except requests.RequestException as e:
            print(f"Failed to fetch GDC manifests: {e}", file=sys.stderr)
            raise

        open_ids = _open_file_id_set(open_map)
        combined: dict[str, tuple[str, int, str]] = {}
        combined.update(controlled_map)
        combined.update(open_map)

        for name in MC3_TARGET_FILENAMES:
            if name not in combined:
                print(
                    f"[skip] {name!r} not found in open or controlled MC3 manifests "
                    "(see https://gdc.cancer.gov/about-data/publications/mc3-2017 ).",
                    file=sys.stderr,
                )
                continue

            file_id, size_b, md5_hex = combined[name]
            controlled = file_id not in open_ids
            if controlled and skip_controlled:
                print(f"[skip-controlled] {name} (requires GDC token)", file=sys.stderr)
                continue
            if controlled and not (gdc_token or "").strip():
                print(
                    f"[skip] {name} is controlled-access; set GDC_TOKEN or pass --gdc-token "
                    "(or use --mc3-open-only).",
                    file=sys.stderr,
                )
                continue

            dest = out_path / name
            if (
                not overwrite
                and dest.is_file()
                and dest.stat().st_size == size_b
            ):
                print(f"[exists] {name}")
                continue

            print(f"[download] {name} ({size_b} bytes) …")
            try:
                _download_one(
                    sess,
                    file_id,
                    dest,
                    size_b,
                    md5_hex,
                    gdc_token if controlled else None,
                    verify_md5,
                )
            except PermissionError as e:
                print(f"[error] {name}: {e}", file=sys.stderr)
                continue
            print(f"[done] {name}")
    finally:
        if own_session:
            sess.close()


def _cli() -> None:
    parser = argparse.ArgumentParser(description="Download MC3 GDC publication supplements.")
    parser.add_argument(
        "--out",
        default="data/raw/mc3",
        help="Output directory (default: data/raw/mc3).",
    )
    parser.add_argument(
        "--gdc-token",
        default=os.environ.get("GDC_TOKEN"),
        help="GDC download token (or set env GDC_TOKEN).",
    )
    parser.add_argument(
        "--mc3-open-only",
        action="store_true",
        help="Only download open-access targets (skips controlled files).",
    )
    parser.add_argument(
        "--verify-md5",
        action="store_true",
        help="Verify MD5 after each file (slow for large archives).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-download even if file exists with correct size.",
    )
    args = parser.parse_args()
    repo = Path(__file__).resolve().parent.parent
    out = Path(args.out)
    if not out.is_absolute():
        out = repo / out
    download_mc3_publication_supplements(
        out,
        gdc_token=args.gdc_token,
        skip_controlled=args.mc3_open_only,
        verify_md5=args.verify_md5,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    _cli()
