# Usage: python plot.py  (reads results/*.json, writes docs/src/assets/b_nl_cd.png and b_nl_cv.png)
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
cfg = json.load(open(os.path.join(here, "common_sizes.json")))
res = {}
for t in (1, 8):
    res[t] = {}
    for lang in ("julia", "python"):
        f = os.path.join(here, "results", f"{lang}_{t}t.json")
        for name, r in json.load(open(f)).items():
            r["lang"] = lang
            res[t][name] = r

# fixed color per package (identity), line style per language
colors = {"CellListMap": "#2a78d6", "NearestNeighbors": "#eb6834", "PointNeighbors": "#1baf7a",
          "scipy": "#eda100", "vesin": "#e87ba4", "freud": "#008300", "matscipy": "#6250d6"}
order = list(colors)

def label(name, r):
    return f"{name}{'.jl' if r['lang'] == 'julia' else ''} {r['version']}"

titles = {"cd": f"Constant density (ρ = {cfg['density']} atoms/Å³, water), cutoff = {cfg['cutoff']:g} Å",
          "cv": f"Constant volume (V = {cfg['cv_side']:g}³ Å³), cutoff = {cfg['cutoff']:g} Å"}

plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#888", "axes.labelcolor": "#222",
                     "xtick.color": "#444", "ytick.color": "#444"})
for case in ("cd", "cv"):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    for ax, t in zip(axes, (1, 8)):
        for name in order:
            # scipy has no parallel pair search (only full lists through query_ball_point)
            if name not in res[t] or (name == "scipy" and t > 1):
                continue
            r = res[t][name]
            data = r[case]
            N = [d[0] for d in data]; tm = [d[1] for d in data]
            lw = 2.6 if name == "CellListMap" else 1.6
            ax.plot(N, tm, color=colors[name], lw=lw, ls="-" if r["lang"] == "julia" else "--",
                    marker="o", ms=4.5, label=label(name, r), zorder=3 if name == "CellListMap" else 2)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.grid(True, which="major", color="#e4e4e2", lw=0.8); ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.set_title("Serial" if t == 1 else "8 threads", fontsize=11)
        ax.set_xlabel("Number of particles")
        if t > 1:
            ax.text(0.98, 0.03, "scipy: no parallel pair search", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=8.5, color="#666")
    axes[0].set_ylabel("Time to build the neighbor list / s")
    # legend from the serial panel, which has every package
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(titles[case], fontsize=12)
    fig.tight_layout(rect=(0, 0.11, 1, 1))
    fig.savefig(os.path.join(here, "..", "..", "docs", "src", "assets", f"b_nl_{case}.png"), dpi=130)
    print("wrote", f"b_nl_{case}.png")
