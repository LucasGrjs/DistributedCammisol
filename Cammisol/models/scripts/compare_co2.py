#!/usr/bin/env python3
"""Compare cumulative CO2 emissions between a centralized and a distributed CAMMISOL run.

Usage (from anywhere, after running both simulations):

    python3 compare_co2.py

By default it reads:
  - ../output.log/results_central/CO2.csv                  (centralized run)
  - ../output.log/results/CO2/CO2_<rank>.csv (one per rank) (distributed run)

and writes a merged CSV + a plot under ../output.log/comparison/.
"""

import argparse
import csv
import re
from pathlib import Path


def read_co2_series(path: Path):
    """Read a CO2 csv produced by cammisol.gaml / Distribution_CAMMISOL.gaml.

    The first one or two lines are header text (not numeric); every other line
    holds one cumulative CO2 value (grams), one per saved cycle.
    """
    values = []
    with path.open(newline="") as f:
        for row in csv.reader(f):
            if not row:
                continue
            cell = row[0].strip().strip("'\"")
            try:
                values.append(float(cell))
            except ValueError:
                continue  # header line
    return values


def collect_distributed_series(directory: Path):
    """Read every CO2_<rank>.csv in directory, return {rank: [values...]}."""
    series = {}
    for file in sorted(directory.glob("CO2_*.csv")):
        match = re.search(r"CO2_(\d+)\.csv$", file.name)
        if not match:
            continue
        series[int(match.group(1))] = read_co2_series(file)
    return series


def sum_distributed(series: dict):
    """Sum the per-rank series index by index.

    Ranks whose nematodes all migrated away (or that finished a bit earlier)
    can end up with fewer saved rows than the others; their last known
    cumulative value is carried forward so shorter series don't silently
    undercount the total once they run out of rows.
    """
    if not series:
        return []
    length = max(len(values) for values in series.values())
    total = [0.0] * length
    for values in series.values():
        if not values:
            continue
        last = values[-1]
        for i in range(length):
            total[i] += values[i] if i < len(values) else last
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    default_output_log = Path(__file__).resolve().parent.parent / "output.log"
    parser.add_argument("--central", type=Path, default=default_output_log / "results_central" / "CO2.csv",
                         help="Path to the centralized run's CO2.csv")
    parser.add_argument("--distributed-dir", type=Path, default=default_output_log / "results" / "CO2",
                         help="Directory containing the distributed run's CO2_<rank>.csv files")
    parser.add_argument("--step", type=int, default=10,
                         help="Number of cycles between two saved rows (must match the 'cycle mod N = 0' condition in the .gaml files)")
    parser.add_argument("--out-dir", type=Path, default=default_output_log / "comparison",
                         help="Where to write the comparison CSV/plot")
    parser.add_argument("--no-plot", action="store_true", help="Skip the matplotlib plot, only print/save the numeric comparison")
    args = parser.parse_args()

    if not args.central.exists():
        parser.error(f"central CO2 file not found: {args.central}")
    if not args.distributed_dir.exists():
        parser.error(f"distributed CO2 directory not found: {args.distributed_dir}")

    central = read_co2_series(args.central)
    per_rank = collect_distributed_series(args.distributed_dir)
    if not per_rank:
        parser.error(f"no CO2_<rank>.csv files found in {args.distributed_dir}")
    distributed = sum_distributed(per_rank)

    n = min(len(central), len(distributed))
    if n == 0:
        parser.error("no comparable data points (one of the two series is empty)")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "CO2_comparison.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cycle", "central_CO2", "distributed_CO2", "abs_diff", "relative_diff_pct"])
        for i in range(n):
            c, d = central[i], distributed[i]
            diff = d - c
            rel = (diff / c * 100) if c != 0 else float("nan")
            writer.writerow([i * args.step, c, d, diff, rel])

    print(f"ranks found          : {sorted(per_rank)} ({len(per_rank)} processes)")
    print(f"central samples      : {len(central)}")
    print(f"distributed samples  : {len(distributed)} (summed across ranks)")
    print(f"compared samples     : {n} (cycle 0 .. {(n - 1) * args.step})")
    print()
    final_c, final_d = central[n - 1], distributed[n - 1]
    final_diff = final_d - final_c
    final_rel = (final_diff / final_c * 100) if final_c != 0 else float("nan")
    print(f"final central CO2    : {final_c:.6e} g")
    print(f"final distributed CO2: {final_d:.6e} g")
    print(f"final abs diff       : {final_diff:.6e} g")
    print(f"final relative diff  : {final_rel:.4f} %")
    print()
    print(f"comparison csv written to: {csv_path}")

    if args.no_plot:
        return

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed, skipping plot (pip install matplotlib)")
        return

    cycles = [i * args.step for i in range(n)]
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(cycles, central[:n], marker="o", markersize=3, label="Centralized")
    ax1.plot(cycles, distributed[:n], marker="s", markersize=3, label="Distributed (sum of ranks)")
    ax1.set_xlabel("Cycle")
    ax1.set_ylabel("Cumulative nematode CO2 emissions (g)")
    ax1.set_title("Centralized vs Distributed CAMMISOL - CO2 emissions")
    ax1.grid(True, linestyle=":", linewidth=0.7, alpha=0.7)
    ax1.legend(loc="upper left")

    ax2 = ax1.twinx()
    relative = [((distributed[i] - central[i]) / central[i] * 100) if central[i] != 0 else 0.0 for i in range(n)]
    ax2.plot(cycles, relative, color="gray", linestyle="--", linewidth=1, alpha=0.6, label="Relative diff (%)")
    ax2.set_ylabel("Relative difference (%)", color="gray")
    ax2.legend(loc="upper right")

    png_path = args.out_dir / "CO2_comparison.png"
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    print(f"comparison plot written to: {png_path}")


if __name__ == "__main__":
    main()
