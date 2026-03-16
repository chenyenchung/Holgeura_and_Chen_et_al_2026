#!/usr/bin/env python3
"""Generate a demo-sized preprocessing input dataset.

This script subsamples raw synapses per neuropil from the original FlyWire
input table, then keeps matching annotation rows and relevant cell-type rows.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import time
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple


SYN_REQUIRED_COLUMNS = (
    "pre_root_id_720575940",
    "post_root_id_720575940",
    "neuropil",
)
ANN_REQUIRED_COLUMNS = ("pre_root_id", "post_root_id", "neuropil")
TYPE_REQUIRED_COLUMNS = ("root_id",)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create demo preprocessing input files by neuropil subsampling."
    )
    parser.add_argument(
        "--syn",
        default="data/fafb_v783_princeton_synapse_table.csv.gz",
        help="Path to raw synapse coordinate CSV(.gz)",
    )
    parser.add_argument(
        "--ann",
        default="data/connections_princeton_no_threshold.csv.gz",
        help="Path to synapse-count annotation CSV(.gz)",
    )
    parser.add_argument(
        "--types",
        default="data/consolidated_cell_types.csv.gz",
        help="Path to cell-type annotation CSV(.gz)",
    )
    parser.add_argument(
        "--outdir",
        default="data/demo_data",
        help="Output directory for demo input files",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Random seed for reproducible subsampling",
    )
    parser.add_argument(
        "--me-count",
        type=int,
        default=50000,
        help="Rows to sample per medulla side (ME_L/ME_R)",
    )
    parser.add_argument(
        "--lo-count",
        type=int,
        default=25000,
        help="Rows to sample per lobula side (LO_L/LO_R)",
    )
    parser.add_argument(
        "--lop-count",
        type=int,
        default=25000,
        help="Rows to sample per lobula plate side (LOP_L/LOP_R)",
    )
    parser.add_argument(
        "--keep-full-types",
        action="store_true",
        help="Write the full type table instead of keeping only referenced root IDs",
    )
    return parser.parse_args()


def strip_root_prefix(root_id: str) -> str:
    prefix = "720575940"
    if root_id.startswith(prefix):
        return root_id[len(prefix):]
    return root_id


def open_reader(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("r", newline="")


def open_writer(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        return gzip.open(path, "wt", newline="")
    return path.open("w", newline="")


def require_columns(header: Sequence[str], required: Sequence[str], file_label: str) -> None:
    missing = [col for col in required if col not in header]
    if missing:
        raise ValueError(
            f"{file_label} is missing required columns: {', '.join(missing)}"
        )


def reservoir_sample_synapses(
    syn_path: Path,
    target_counts: Dict[str, int],
    seed: int,
) -> Tuple[List[str], Dict[str, List[dict]], Counter]:
    rng = random.Random(seed)
    seen = Counter()
    reservoirs = {np: [] for np in target_counts}

    with open_reader(syn_path) as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Empty synapse input file: {syn_path}")
        syn_header = list(reader.fieldnames)
        require_columns(syn_header, SYN_REQUIRED_COLUMNS, str(syn_path))

        for row in reader:
            npil = row["neuropil"]
            if npil not in target_counts:
                continue

            seen[npil] += 1
            n_seen = seen[npil]
            k = target_counts[npil]
            bucket = reservoirs[npil]

            if len(bucket) < k:
                bucket.append(row)
                continue

            j = rng.randrange(n_seen)
            if j < k:
                bucket[j] = row

    for npil, k in target_counts.items():
        if seen[npil] < k:
            raise ValueError(
                f"Not enough rows for {npil}: requested {k}, found {seen[npil]}"
            )

    return syn_header, reservoirs, seen


def build_selection_sets(
    reservoirs: Dict[str, List[dict]],
) -> Tuple[Set[Tuple[str, str, str]], Set[str], Dict[str, int]]:
    pair_keys: Set[Tuple[str, str, str]] = set()
    root_ids: Set[str] = set()
    sampled_counts: Dict[str, int] = {}

    for npil, rows in reservoirs.items():
        sampled_counts[npil] = len(rows)
        for row in rows:
            pre = row["pre_root_id_720575940"]
            post = row["post_root_id_720575940"]
            pair_keys.add((pre, post, npil))
            root_ids.add(pre)
            root_ids.add(post)

    return pair_keys, root_ids, sampled_counts


def write_synapse_subset(
    output_path: Path,
    syn_header: Sequence[str],
    reservoirs: Dict[str, List[dict]],
    neuropil_order: Iterable[str],
) -> Counter:
    per_neuropil = Counter()
    with open_writer(output_path) as handle:
        writer = csv.DictWriter(handle, fieldnames=syn_header)
        writer.writeheader()

        for npil in neuropil_order:
            for row in reservoirs[npil]:
                writer.writerow(row)
                per_neuropil[npil] += 1
    return per_neuropil


def write_annotation_subset(
    ann_path: Path,
    output_path: Path,
    pair_keys: Set[Tuple[str, str, str]],
    target_counts: Dict[str, int],
) -> Tuple[int, Counter]:
    n_written = 0
    per_neuropil = Counter()

    with open_reader(ann_path) as in_handle, open_writer(output_path) as out_handle:
        reader = csv.DictReader(in_handle)
        if reader.fieldnames is None:
            raise ValueError(f"Empty annotation input file: {ann_path}")
        ann_header = list(reader.fieldnames)
        require_columns(ann_header, ANN_REQUIRED_COLUMNS, str(ann_path))

        writer = csv.DictWriter(out_handle, fieldnames=ann_header)
        writer.writeheader()

        for row in reader:
            npil = row["neuropil"]
            if npil not in target_counts:
                continue

            key = (
                strip_root_prefix(row["pre_root_id"]),
                strip_root_prefix(row["post_root_id"]),
                npil,
            )
            if key in pair_keys:
                writer.writerow(row)
                n_written += 1
                per_neuropil[npil] += 1

    return n_written, per_neuropil


def write_type_subset(
    types_path: Path,
    output_path: Path,
    root_ids: Set[str],
    keep_full_types: bool,
) -> int:
    n_written = 0

    with open_reader(types_path) as in_handle, open_writer(output_path) as out_handle:
        reader = csv.DictReader(in_handle)
        if reader.fieldnames is None:
            raise ValueError(f"Empty type input file: {types_path}")
        type_header = list(reader.fieldnames)
        require_columns(type_header, TYPE_REQUIRED_COLUMNS, str(types_path))

        writer = csv.DictWriter(out_handle, fieldnames=type_header)
        writer.writeheader()

        if keep_full_types:
            for row in reader:
                writer.writerow(row)
                n_written += 1
            return n_written

        for row in reader:
            rid = strip_root_prefix(row["root_id"])
            if rid in root_ids:
                writer.writerow(row)
                n_written += 1

    return n_written


def main() -> None:
    args = parse_args()
    start = time.time()

    syn_path = Path(args.syn)
    ann_path = Path(args.ann)
    types_path = Path(args.types)
    outdir = Path(args.outdir)

    for in_path in (syn_path, ann_path, types_path):
        if not in_path.exists():
            raise FileNotFoundError(f"Input file not found: {in_path}")

    target_counts = {
        "ME_L": args.me_count,
        "ME_R": args.me_count,
        "LO_L": args.lo_count,
        "LO_R": args.lo_count,
        "LOP_L": args.lop_count,
        "LOP_R": args.lop_count,
    }
    neuropil_order = ["ME_L", "ME_R", "LO_L", "LO_R", "LOP_L", "LOP_R"]

    print("Sampling synapse rows...")
    syn_header, reservoirs, n_seen = reservoir_sample_synapses(
        syn_path=syn_path,
        target_counts=target_counts,
        seed=args.seed,
    )

    pair_keys, root_ids, sampled_counts = build_selection_sets(reservoirs)

    syn_out = outdir / "fafb_v783_princeton_synapse_table.csv.gz"
    ann_out = outdir / "connections_princeton_no_threshold.csv.gz"
    types_out = outdir / "consolidated_cell_types.csv.gz"

    print(f"Writing synapse subset to {syn_out}...")
    syn_written = write_synapse_subset(
        output_path=syn_out,
        syn_header=syn_header,
        reservoirs=reservoirs,
        neuropil_order=neuropil_order,
    )

    print(f"Writing annotation subset to {ann_out}...")
    ann_written_total, ann_written_by_neuropil = write_annotation_subset(
        ann_path=ann_path,
        output_path=ann_out,
        pair_keys=pair_keys,
        target_counts=target_counts,
    )

    print(f"Writing type subset to {types_out}...")
    types_written = write_type_subset(
        types_path=types_path,
        output_path=types_out,
        root_ids=root_ids,
        keep_full_types=args.keep_full_types,
    )

    manifest = {
        "seed": args.seed,
        "targets": target_counts,
        "source_counts": dict(n_seen),
        "sampled_synapse_counts": dict(sampled_counts),
        "written_synapse_counts": dict(syn_written),
        "annotation_rows_written": ann_written_total,
        "annotation_rows_by_neuropil": dict(ann_written_by_neuropil),
        "types_rows_written": types_written,
        "unique_sampled_pairs": len(pair_keys),
        "unique_sampled_root_ids": len(root_ids),
        "files": {
            "syn": str(syn_out),
            "ann": str(ann_out),
            "types": str(types_out),
        },
    }

    manifest_path = outdir / "manifest.json"
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    elapsed = time.time() - start
    print("Done.")
    print(json.dumps(manifest, indent=2))
    print(f"Elapsed: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
