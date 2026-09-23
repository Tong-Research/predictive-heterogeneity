"""Score GENDESIGN against prereg/GENDESIGN.md.

The predictions are fixed in the prereg; this only reads them off the results. G3 and G4 are
instrument checks: if either fails nothing is claimed from the run, so they are reported first
and the script says so rather than printing a verdict on the rest.
"""
import glob
import pathlib
import sys

import pandas as pd

R = pathlib.Path(__file__).resolve().parents[1] / "results/gendesign"
HET_HI = 1.5


def load():
    fs = sorted(glob.glob(str(R / "*.csv")))
    if not fs:
        print(f"no results in {R} -- nothing to score, which is NOT a pass")
        sys.exit(2)
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    return d


def main() -> int:
    d = load()
    mix = d[d.outcome_from == "mix"].copy()
    flat = d[d.outcome_from == "shared_flat"].copy()
    rhos = sorted(mix.rho.dropna().unique())
    print(f"{len(d)} configurations; mix arms rho={rhos}; shared_flat rows={len(flat)}\n")

    gain = mix.pivot_table(index="het", columns="rho", values="delta_vs_tuned_ind", aggfunc="mean")
    print("safe-adaptive gain over the tuned indicator, by slope heterogeneity and dose")
    print(gain.round(4).to_string(), "\n")

    spread = mix.pivot_table(index=["mechanism", "het", "info", "concentration", "n"],
                             columns="rho", values="coef_spread", aggfunc="mean")
    rel = ((spread.max(axis=1) - spread.min(axis=1)) / spread.max(axis=1).abs()).max()
    print(f"largest relative spread of coef_spread across rho within a cell: {rel:.2e}\n")

    checks, instrument_ok = {}, True
    # --- instrument -------------------------------------------------------------------
    g4 = bool(rel < 1e-6)
    checks["G4 coefficient spread constant across rho (<1e-6)"] = g4
    instrument_ok &= g4
    if HET_HI in gain.index and 0.0 in gain.columns:
        g3 = bool(gain.loc[:, 0.0].abs().max() <= 0.01)
    else:
        g3 = False
    checks["G3 control (rho=0) within +/-0.01 of zero at every het"] = g3
    instrument_ok &= g3
    # --- predictions ------------------------------------------------------------------
    if HET_HI in gain.index:
        row = gain.loc[HET_HI].sort_index()
        checks["G1 gain monotone non-decreasing in rho at het=1.5"] = bool(row.is_monotonic_increasing)
        checks["G2 gain(rho=1) - gain(rho=0) > +0.10 at het=1.5"] = bool(row.iloc[-1] - row.iloc[0] > 0.10)
        if 0.5 in row.index:
            checks["G6 gain(rho=0.5) strictly between the endpoints at het=1.5"] = bool(
                row.iloc[0] < row.loc[0.5] < row.iloc[-1])
    if len(flat):
        f = flat.groupby("het").delta_vs_tuned_ind.mean()
        checks["G5 shared_flat within +/-0.01 of zero at every het"] = bool((f.abs() <= 0.01).all())
        checks["G5b shared_flat does not rise with het"] = bool(
            not f.is_monotonic_increasing or f.max() - f.min() <= 0.01)
        print("shared_flat gain by het"); print(f.round(4).to_string(), "\n")

    for k, v in checks.items():
        print(("  HELD    " if v else "  FAILED  ") + k)
    print()
    if not instrument_ok:
        print("INSTRUMENT FAILED — per prereg/GENDESIGN.md nothing is claimed from this run.")
        return 1
    failed = [k for k, v in checks.items() if not v]
    print("ALL HELD" if not failed else f"NOT ALL HELD: {failed}\n-> apply the withdrawal condition in prereg/GENDESIGN.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
