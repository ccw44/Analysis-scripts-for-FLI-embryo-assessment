#!/usr/bin/env python3
"""Stage the vendor protein result package and write the sample sheet.

    python prepare_inputs.py results/proteomics [--archive proteinomics.zip]

With --archive the tabular files are extracted from Result_N0011/ and hashed; without it
an already extracted tree is used and the same manifest is rebuilt. Fields the vendor file
does not state are written as ``unknown``; nothing is inferred.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import zipfile
from pathlib import Path

SAMPLES = ["H1", "H2", "H3", "L1", "L2", "L3"]
SUFFIXES = (".txt", ".csv", ".json")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def from_archive(archive: Path, base: Path):
    manifest = []
    with zipfile.ZipFile(archive) as zf:
        for name in zf.namelist():
            if name.startswith("Result_N0011/") and name.endswith(SUFFIXES):
                data = zf.read(name)
                target = base / "inputs" / name
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(data)
                manifest.append({"file": name, "bytes": len(data),
                                 "sha256": hashlib.sha256(data).hexdigest()})
        for name in zf.namelist():
            if name.endswith("\u6837\u672c\u4fe1\u606f.txt"):
                (base / "inputs" / "sample_info_source.txt").write_bytes(zf.read(name))
    return manifest, sha256_file(archive)


def from_tree(root: Path, base: Path):
    if not root.is_dir():
        raise SystemExit(f"vendor result directory not found: {root}")
    manifest = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.suffix.lower() in SUFFIXES:
            manifest.append({
                "file": str(path.relative_to(base / "inputs")),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            })
    return manifest, "not-applicable (extracted tree used in place)"


def write_sample_sheet(base: Path) -> Path:
    """HS = high score, LS = low score; independent single-embryo samples, not paired.

    A shared numeric suffix is a label, not evidence of pairing, so pair_id is
    ``not_paired``.
    """
    out = base / "inputs" / "sample_metadata.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "group", "biological_unit", "pair_id",
                    "developmental_stage", "source", "batch"])
        for s in SAMPLES:
            w.writerow([s, "HS" if s.startswith("H") else "LS",
                        "independent_single_embryo_user_confirmed", "not_paired",
                        "blastocyst_filename_only_unverified", "unknown", "unknown"])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Stage the vendor protein result package.")
    ap.add_argument("base", type=Path, help="output base, e.g. results/proteomics")
    ap.add_argument("--archive", type=Path, default=None,
                    help="vendor archive; omit to use an extracted tree under <base>/inputs")
    a = ap.parse_args()

    base = a.base.resolve()
    (base / "inputs").mkdir(parents=True, exist_ok=True)

    if a.archive is not None:
        manifest, archive_sha = from_archive(a.archive.resolve(), base)
        mode = "archive"
    else:
        manifest, archive_sha = from_tree(base / "inputs" / "Result_N0011", base)
        mode = "in place"

    (base / "inputs" / "manifest.json").write_text(
        json.dumps({"archive_sha256": archive_sha, "files": manifest}, indent=2),
        encoding="utf-8")
    sheet = write_sample_sheet(base)

    print(f"mode         : {mode}")
    print(f"source files : {len(manifest)}")
    print(f"manifest     : {base / 'inputs' / 'manifest.json'}")
    print(f"sample sheet : {sheet}  ({len(SAMPLES)} samples)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
