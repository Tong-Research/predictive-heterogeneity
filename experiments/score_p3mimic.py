"""Score prereg/P3MIMIC.md (Paper 3's MIMIC-IV test replicated on mimic4_sites.npz). Committed before any
result row. Reads results/p3mimic_suite.csv (run_real_suite --only mimic4 --seeds 2) and
results/p3mimic_ceiling.csv (oracle_ceiling_b --datasets mimic4). Refuses missing files.
Predictions: 1 tie (|delta| < 0.002 or p >= 0.05); 2 safe_oracle == tuned_indicator within 0.002;
3 H_exc > 0; 4 ceiling excess <= 0."""
from __future__ import annotations
import pathlib, sys
import pandas as pd
R = pathlib.Path(__file__).resolve().parent.parent / "results"

def main() -> int:
    s, c = R / "p3mimic_suite.csv", R / "p3mimic_ceiling.csv"
    for f in (s, c):
        if not f.exists(): raise SystemExit(f"  {f.name} MISSING -- refusing to score")
    a = pd.read_csv(s).iloc[0]; b = pd.read_csv(c).iloc[0]
    print(f"  mimic4 n={int(a.n)} d={int(a.d)}  iota {a.iota:+.3f}  H {a.het:.3f}  H_exc {a.het_excess:+.3f}  c30 {a.cov30:.3f}")
    print(f"  tuned_impute {a.tuned_impute:.4f}  tuned_indicator {a.tuned_indicator:.4f}  safe_adaptive {a.safe_adaptive:.4f}  safe_oracle {a.safe_oracle:.4f}")
    print(f"  delta (safe_adaptive - tuned_indicator) {a.delta:+.4f}  p={a.p:.3g}  wins {int(a.wins)}/{int(a.n_folds)}")
    v = {}
    v["P1 tie"] = (abs(a.delta) < 0.002) or (a.p >= 0.05)
    v["P2 oracle=indicator"] = abs(a.safe_oracle - a.tuned_indicator) < 0.002
    v["P3 H_exc>0"] = a.het_excess > 0
    v["P4 ceiling excess<=0"] = b.excess <= 0
    print(f"  ceiling real {b.ceiling_real:+.4f}  permuted mean {b.perm_mean:+.4f} [{b.perm_lo:+.4f}, {b.perm_hi:+.4f}]  excess {b.excess:+.4f}  p={b.p:.3f}")
    for k, ok in v.items(): print(f"    {k:24s} {'HOLDS' if ok else 'FAILS'}")
    if a.delta >= 0.002 and a.p < 0.05: print("  WITHDRAWAL: the method beats the indicator here -- Paper 3's decisive case is not decisive on this cohort; rewrite both extracts honestly.")
    if a.het_excess <= 0: print("  WITHDRAWAL: diagnostic did not flag this cohort; withdraw 'flagged in advance' for this matrix.")
    print("  Always: state the extract (42 labs, n, prevalence) beside every MIMIC-IV number.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
