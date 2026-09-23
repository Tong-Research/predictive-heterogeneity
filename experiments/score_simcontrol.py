"""Score SIMCONTROL against prereg/SIMCONTROL.md."""
import glob, pathlib, sys, numpy as np, pandas as pd
R = pathlib.Path(__file__).resolve().parents[1] / "results/simcontrol"
d = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(R / "*.csv")))], ignore_index=True)
print(f"{len(d)} configurations, generators {sorted(d.outcome_from.unique())}, "
      f"mechanisms {sorted(d.mechanism.unique())}\n")
g = d.pivot_table(index="het", columns="outcome_from", values="delta_vs_tuned_ind",
                  aggfunc=["mean", "sem"])
print("safe-adaptive gain over the tuned indicator baseline, by slope heterogeneity")
print(g.round(4).to_string(), "\n")
sp = d.pivot_table(index="het", columns="outcome_from", values="coef_spread", aggfunc="mean")
sp["rel_diff"] = (sp["pattern"] - sp["shared"]).abs() / sp["pattern"]
print("realised coefficient spread (the control must match the standard generator)")
print(sp.round(4).to_string(), "\n")
m = d.pivot_table(index=["mechanism", "het"], columns="outcome_from", values="delta_vs_tuned_ind")
print("by mechanism"); print(m.round(4).to_string(), "\n")
P = d[d.outcome_from == "pattern"].groupby("het").delta_vs_tuned_ind.mean()
S = d[d.outcome_from == "shared"].groupby("het").delta_vs_tuned_ind.mean()
hi = max(P.index)
ok = {
 f"C1 standard gain rises with het and exceeds +0.05 at het={hi}":
   bool(P.is_monotonic_increasing and P.loc[hi] > 0.05),
 "C2 control gain within +/-0.01 of zero at every het": bool((S.abs() <= 0.01).all()),
 f"C3 difference exceeds +0.10 at het={hi}": bool(P.loc[hi] - S.loc[hi] > 0.10),
 "C4 coefficient spread agrees to within 1%": bool((sp["rel_diff"] <= 0.01).all()),
 "C5 control gain does not rise with het": bool(not S.is_monotonic_increasing or S.max() - S.min() <= 0.01),
}
for k, v in ok.items(): print(("  HELD    " if v else "  FAILED  ") + k)
print("\n" + ("ALL FIVE HELD" if all(ok.values()) else "NOT ALL HELD -- see prereg/SIMCONTROL.md withdrawal"))
sys.exit(0)
