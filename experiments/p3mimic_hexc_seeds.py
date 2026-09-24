"""het_excess on the committed MIMIC-IV matrix over five permutation seeds.

The heterogeneity paper quotes the divergence on this matrix as a mean and sd over permutation
seeds. That figure was measured on 2026-09-02 and written into prereg/P3MIMIC.md, but no results
file held it, so it could not be recomputed. This writes results/p3mimic_hexc_seeds.csv from the
same code path as the registered run: run_real_suite._mimic4 (40,000-record subsample, seed 0)
and run_real_suite.diagnostics, whose seed argument sets only the mask permutation. Seed 0 must
reproduce the registered het_excess in results/p3mimic_suite.csv, or the file is not written.
"""
import csv
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_real_suite import _mimic4, diagnostics  # noqa: E402

R = HERE.parent / "results"
SEEDS = range(5)


def main():
    reg = next(csv.DictReader(open(R / "p3mimic_suite.csv")))
    X, y = _mimic4(max_n=40_000, seed=0)
    rows = []
    for s in SEEDS:
        d = diagnostics(X, y, seed=s)
        rows.append(dict(seed=s, het=d["het"], het_null=d["het_null"], het_excess=d["het_excess"]))
        print(f"  seed {s}: het {d['het']:.4f}  null {d['het_null']:.4f}  excess {d['het_excess']:+.4f}", flush=True)
    # the fits are not bit-reproducible across runs (differences near 1e-7 were seen); the paper
    # quotes three decimals, so agreement to 1e-6 is the same computation
    if abs(rows[0]["het_excess"] - float(reg["het_excess"])) > 1e-6:
        raise SystemExit(f"seed 0 gives {rows[0]['het_excess']:+.6f}, the registered run "
                         f"{float(reg['het_excess']):+.6f}: not the same computation, nothing written")
    with open(R / "p3mimic_hexc_seeds.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    e = np.array([r["het_excess"] for r in rows])
    print(f"wrote results/p3mimic_hexc_seeds.csv: mean {e.mean():+.3f}, sd {e.std():.3f} "
          f"(ddof=0, the convention of the paper's other +/- figures; ddof=1 gives {e.std(ddof=1):.3f})")


if __name__ == "__main__":
    main()
