import sys
import csv
import math
from collections import defaultdict
import matplotlib.pyplot as plt

path = sys.argv[1] if len(sys.argv) > 1 else "defect_plane_atlas.csv"
rows = []
with open(path, newline="") as fp:
    for r in csv.DictReader(fp):
        def val(k):
            try:
                return float(r[k])
            except Exception:
                return float("nan")
        r["_sum"] = val("Theta_sum")
        r["_shape"] = val("Theta_shape")
        r["_theta"] = val("Theta_NM")
        r["_du"] = val("du_rms")
        r["_phi"] = val("phi")
        rows.append(r)

classes = sorted(set(r["sample_class"] for r in rows))
markers = {"voro_radical": "o", "post_surgery_nonradical": "s", "abstract": "^", "not_run": "x"}

plt.figure(figsize=(7, 6))
for c in classes:
    xs = [r["_sum"] for r in rows if r["sample_class"] == c and math.isfinite(r["_sum"])]
    ys = [r["_shape"] for r in rows if r["sample_class"] == c and math.isfinite(r["_shape"])]
    plt.scatter(xs, ys, label=c, s=36)
plt.axvline(0.05, color="k", lw=1, ls="--")
plt.axhline(0.05, color="k", lw=1, ls="--")
plt.xlabel("Theta_sum")
plt.ylabel("Theta_shape")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig("defect_plane.png", dpi=200)

def scatter(xkey, ykey, out, xlabel, ylabel):
    plt.figure(figsize=(6, 5))
    for c in classes:
        xs = [r[xkey] for r in rows if r["sample_class"] == c and math.isfinite(r[xkey]) and math.isfinite(r[ykey])]
        ys = [r[ykey] for r in rows if r["sample_class"] == c and math.isfinite(r[xkey]) and math.isfinite(r[ykey])]
        plt.scatter(xs, ys, label=c, s=30)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out, dpi=200)

scatter("_theta", "_du", "theta_vs_du_rms.png", "Theta_NM", "du_rms")
scatter("_sum", "_du", "theta_sum_vs_du_rms.png", "Theta_sum", "du_rms")
scatter("_shape", "_du", "theta_shape_vs_du_rms.png", "Theta_shape", "du_rms")

for cls, out in [("voro_soft_phi_scan", "phi_trajectory.png"), ("hex_coord_noise", "coord_noise_trajectory.png"),
                 ("hex_radius_noise", "radius_noise_trajectory.png")]:
    pts = [r for r in rows if r["sample_class"] == cls and math.isfinite(r["_sum"]) and math.isfinite(r["_shape"])]
    pts.sort(key=lambda r: r["_phi"] if cls == "voro_soft_phi_scan" else float(r["perturb_amp"]))
    if pts:
        plt.figure(figsize=(6, 5))
        plt.plot([r["_sum"] for r in pts], [r["_shape"] for r in pts], marker="o")
        for r in pts:
            label = f"{r['_phi']:.2f}" if cls == "voro_soft_phi_scan" else f"{float(r['perturb_amp']):.2f}"
            plt.text(r["_sum"], r["_shape"], label, fontsize=8)
        plt.axvline(0.05, color="k", lw=1, ls="--")
        plt.axhline(0.05, color="k", lw=1, ls="--")
        plt.xlabel("Theta_sum")
        plt.ylabel("Theta_shape")
        plt.tight_layout()
        plt.savefig(out, dpi=200)
