import sys, csv, math
from collections import defaultdict
import matplotlib.pyplot as plt

path = sys.argv[1] if len(sys.argv) > 1 else "operator_multistart.csv"
rows = []
with open(path, newline="") as fp:
    for r in csv.DictReader(fp):
        def val(k):
            try: return float(r[k])
            except Exception: return float("nan")
        r["_sx0"] = val("Theta_sum_before")
        r["_sy0"] = val("Theta_shape_before")
        r["_sx1"] = val("Theta_sum_after")
        r["_sy1"] = val("Theta_shape_after")
        rows.append(r)

ops = sorted(set(r["operator_name"] for r in rows))
colors = {op: plt.cm.tab20(i % 20) for i, op in enumerate(ops)}
plt.figure(figsize=(8, 7))
for r in rows:
    if not all(math.isfinite(r[k]) for k in ["_sx0", "_sy0", "_sx1", "_sy1"]):
        continue
    c = colors[r["operator_name"]]
    plt.plot([r["_sx0"], r["_sx1"]], [r["_sy0"], r["_sy1"]], color=c, alpha=0.25)
    plt.scatter([r["_sx1"]], [r["_sy1"]], color=c, s=10)
for op in ops:
    sub = [r for r in rows if r["operator_name"] == op and all(math.isfinite(r[k]) for k in ["_sx0","_sy0","_sx1","_sy1"])]
    if not sub: continue
    x0 = sum(r["_sx0"] for r in sub)/len(sub); y0 = sum(r["_sy0"] for r in sub)/len(sub)
    x1 = sum(r["_sx1"] for r in sub)/len(sub); y1 = sum(r["_sy1"] for r in sub)/len(sub)
    plt.arrow(x0, y0, x1-x0, y1-y0, color=colors[op], width=0.0005, head_width=0.004, length_includes_head=True, label=op)
plt.axvline(0.05, color="k", ls="--", lw=1)
plt.axhline(0.05, color="k", ls="--", lw=1)
plt.xlabel("Theta_sum")
plt.ylabel("Theta_shape")
plt.title("Operator vectors in defect plane")
plt.legend(fontsize=7, loc="best")
plt.tight_layout()
plt.savefig("operator_vectors.png", dpi=180)

metrics = defaultdict(lambda: [0,0,0,0])
for r in rows:
    op = r["operator_name"]
    metrics[op][0] += 1
    metrics[op][1] += int(r.get("success_left","0"))
    metrics[op][2] += int(r.get("success_down","0"))
    metrics[op][3] += int(r.get("success_diagonal","0"))
labels = list(metrics)
x = range(len(labels))
plt.figure(figsize=(max(8, len(labels)*0.5), 4))
for j, name in enumerate(["left","down","diagonal"], start=1):
    vals = [metrics[op][j]/max(metrics[op][0],1) for op in labels]
    plt.bar([i + 0.25*(j-2) for i in x], vals, width=0.25, label=name)
plt.xticks(list(x), labels, rotation=45, ha="right")
plt.ylim(0,1)
plt.ylabel("success fraction")
plt.legend()
plt.tight_layout()
plt.savefig("operator_success_fractions.png", dpi=180)
