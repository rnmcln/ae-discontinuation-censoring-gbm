"""Summarise the manually verified coding in review_coded.csv."""
import pandas as pd

d = pd.read_csv("review_coded.csv")
e = d[d.eligible == 1]
N = len(e)
items = ["ae_disc_reported", "postdisc_survival_fu_stated", "os_ltfu_reported",
         "os_censoring_of_lost_described", "informative_censoring_sensitivity",
         "historical_benchmark", "external_patient_level_control"]
rows = [dict(item=c, n=int(e[c].sum()), N=N, percent=round(100 * e[c].mean(), 1)) for c in items]
for k, v in e.primary_endpoint.value_counts().items():
    rows.append(dict(item="primary_endpoint: " + k, n=int(v), N=N, percent=round(100 * v / N, 1)))
nz = e[e.ae_disc_pct > 0].ae_disc_pct
rows += [dict(item="ae_disc_pct_nonzero_min", n=None, N=len(nz), percent=nz.min()),
         dict(item="ae_disc_pct_nonzero_median", n=None, N=len(nz), percent=round(nz.median(), 1)),
         dict(item="ae_disc_pct_nonzero_max", n=None, N=len(nz), percent=nz.max()),
         dict(item="ae_disc_explicit_zero", n=int((e.ae_disc_pct == 0).sum()), N=N, percent=None),
         dict(item="ae_disc_reported_in_table_only", n=int((e.ae_disc_source == "table").sum()), N=N, percent=None)]
out = pd.DataFrame(rows)
out.to_csv("review_summary.csv", index=False)
print(out.to_string(index=False))
print("excluded:", int((d.eligible == 0).sum()), dict(d[d.eligible == 0].exclusion_reason.value_counts()))
