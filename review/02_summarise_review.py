"""Summarise the trial-level coding in review_trials.csv (one row per trial)."""
import pandas as pd

t = pd.read_csv("review_trials.csv")
rec = pd.read_csv("review_records.csv")
N = len(t)
binary = ["historical_benchmark", "external_patient_level_control", "ae_discontinuation_reported",
          "losses_counted", "handling_of_lost_described", "informative_censoring_sensitivity"]
rows = [dict(item=c, n=int(t[c].sum()), N=N, percent=round(100 * t[c].mean(), 1)) for c in binary]
rows.append(dict(item="benchmark_or_external_control", n=int(((t.historical_benchmark + t.external_patient_level_control) > 0).sum()), N=N,
                 percent=round(100 * ((t.historical_benchmark + t.external_patient_level_control) > 0).mean(), 1)))
for k, v in t.survival_followup_after_discontinuation.value_counts().items():
    rows.append(dict(item="survival_followup_after_discontinuation: " + k, n=int(v), N=N, percent=round(100 * v / N, 1)))
for k, v in t.primary_endpoint.value_counts().items():
    rows.append(dict(item="primary_endpoint: " + k, n=int(v), N=N, percent=round(100 * v / N, 1)))
for col in ["setting", "histology", "phase"]:
    for k, v in t[col].value_counts().items():
        rows.append(dict(item=f"{col}: {k}", n=int(v), N=N, percent=round(100 * v / N, 1)))
a = t[t.ae_discontinuation_pct.notna()]
nz = a[a.ae_discontinuation_pct > 0].ae_discontinuation_pct
rows += [dict(item="ae_disc_pct_determinable", n=len(a), N=N, percent=None),
         dict(item="ae_disc_pct_zero", n=int((a.ae_discontinuation_pct == 0).sum()), N=N, percent=None),
         dict(item="ae_disc_pct_nonzero_min", n=None, N=len(nz), percent=nz.min()),
         dict(item="ae_disc_pct_nonzero_median", n=None, N=len(nz), percent=round(nz.median(), 1)),
         dict(item="ae_disc_pct_nonzero_max", n=None, N=len(nz), percent=nz.max()),
         dict(item="sample_size_median", n=None, N=N, percent=float(t.n.median())),
         dict(item="sample_size_q1", n=None, N=N, percent=float(t.n.quantile(0.25))),
         dict(item="sample_size_q3", n=None, N=N, percent=float(t.n.quantile(0.75))),
         dict(item="sample_size_min", n=None, N=N, percent=float(t.n.min())),
         dict(item="sample_size_max", n=None, N=N, percent=float(t.n.max())),
         dict(item="year_min", n=None, N=N, percent=float(t.year.min())),
         dict(item="year_max", n=None, N=N, percent=float(t.year.max())),
         dict(item="html_full_text_only", n=int(t.text_source.str.startswith("PMC HTML").sum()), N=N, percent=None)]
for (st, dc), v in rec.groupby(["screening_stage", "decision"]).size().items():
    rows.append(dict(item=f"records: {st} | {dc}", n=int(v), N=len(rec), percent=None))
for (st, rs), v in rec[rec.decision == "Excluded"].groupby(["screening_stage", "reason"]).size().items():
    rows.append(dict(item=f"excluded: {st} | {rs}", n=int(v), N=len(rec), percent=None))
out = pd.DataFrame(rows)
out.to_csv("review_summary.csv", index=False)
print(out.to_string(index=False))
