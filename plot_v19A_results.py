import sys
import pandas as pd
import matplotlib.pyplot as plt

traj_path = sys.argv[1] if len(sys.argv) > 1 else "overlap_continuation_trajectory.csv"
cycle_path = "cycle_F_decomposition.csv"
pareto_path = "pareto_archive.csv"

def savefig(name):
    plt.tight_layout()
    plt.savefig(name, dpi=180)
    plt.close()

traj = pd.read_csv(traj_path)
plt.figure()
for col in ["Theta_sum", "Theta_shape", "Theta_NM", "U_overlap"]:
    if col in traj:
        plt.plot(traj["cycle_global"], traj[col], label=col)
plt.xlabel("cycle")
plt.legend()
savefig("v19A_cycle_trajectory.png")

plt.figure()
sc = plt.scatter(traj["Theta_sum"], traj["Theta_shape"], c=traj["w_U"], cmap="viridis", s=24)
plt.xlabel("Theta_sum")
plt.ylabel("Theta_shape")
plt.colorbar(sc, label="w_U")
savefig("v19A_defect_plane.png")

try:
    pareto = pd.read_csv(pareto_path)
    plt.figure()
    plt.scatter(pareto["U_overlap"], pareto["Theta_NM"], s=24)
    plt.xlabel("U_overlap")
    plt.ylabel("Theta_NM")
    savefig("v19A_pareto.png")
except FileNotFoundError:
    pass

try:
    cyc = pd.read_csv(cycle_path)
    plt.figure()
    x = cyc["cycle"]
    bottoms = None
    for col in ["F_sum_after", "F_shape_after", "F_U_after", "F_radius_after"]:
        vals = cyc[col]
        if bottoms is None:
            plt.bar(x, vals, label=col)
            bottoms = vals.copy()
        else:
            plt.bar(x, vals, bottom=bottoms, label=col)
            bottoms = bottoms + vals
    plt.xlabel("cycle")
    plt.legend()
    savefig("v19A_F_decomposition.png")
except FileNotFoundError:
    pass

plt.figure()
plt.scatter(traj["cycle_global"], traj["accepted"], s=12)
plt.xlabel("cycle")
plt.ylabel("accepted")
savefig("v19A_acceptance.png")

try:
    bdu = pd.read_csv("adaptive_bounded_du_trajectory.csv")
    plt.figure()
    plt.plot(bdu["cycle"], bdu["bounded_du"])
    plt.xlabel("cycle")
    plt.ylabel("bounded_du")
    savefig("v19A_bounded_du.png")
except FileNotFoundError:
    pass
