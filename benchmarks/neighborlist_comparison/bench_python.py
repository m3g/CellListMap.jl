# Neighbor list benchmarks (periodic cubic boxes) of Python packages.
# Usage: python bench_python.py <nthreads> <output.json> [package,package,...]
# Results are merged into <output.json> after each package.
import json, sys, time, os
import numpy as np

nthreads = int(sys.argv[1])
outfile = sys.argv[2]
cfg = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "common_sizes.json")))
cutoff = cfg["cutoff"]
TMAX = 60.0  # do not run larger systems once a single run exceeds this

import scipy, vesin, freud, matscipy
from scipy.spatial import cKDTree
from matscipy.neighbours import neighbour_list
freud.parallel.set_num_threads(nthreads)

def run_scipy(x, L):
    t = cKDTree(x, boxsize=L)
    if nthreads == 1:
        p = t.query_pairs(cutoff, output_type="ndarray")
        return len(p)
    # parallel queries are only available through query_ball_point (full lists)
    nb = t.query_ball_point(x, cutoff, workers=nthreads, return_sorted=False)
    return (sum(len(v) for v in nb) - len(x)) // 2

def run_vesin(x, L):
    nl = vesin.NeighborList(cutoff=cutoff, full_list=False, n_threads=nthreads)
    i, j, d = nl.compute(points=x, box=np.diag([L, L, L]), periodic=True, quantities="ijd")
    return len(i)

def run_freud(x, L):
    box = freud.box.Box.cube(L)
    xc = x - L / 2
    lc = freud.locality.AABBQuery(box, xc)
    nl = lc.query(xc, dict(r_max=cutoff, exclude_ii=True)).toNeighborList()
    return len(nl) // 2  # freud returns full lists

def run_matscipy(x, L):
    i, j, d = neighbour_list("ijd", positions=x, cutoff=cutoff, cell=np.diag([L, L, L]), pbc=[True] * 3)
    return len(i) // 2  # full lists

competitors = {"scipy": (run_scipy, scipy.__version__),
               "vesin": (run_vesin, vesin.__version__),
               "freud": (run_freud, freud.__version__)}
if nthreads == 1:
    competitors["matscipy"] = (run_matscipy, matscipy.__version__)

def timeit(f, x, L):
    best = np.inf
    for _ in range(3):
        t0 = time.perf_counter(); n = f(x, L); t = time.perf_counter() - t0
        best = min(best, t)
        if t > 5.0:
            break
    return best, n

if len(sys.argv) > 3:
    competitors = {k: v for k, v in competitors.items() if k in sys.argv[3].split(",")}

rng = np.random.default_rng(321)
results = json.load(open(outfile)) if os.path.exists(outfile) else {}
for name, (f, ver) in competitors.items():
    xw = rng.random((3000, 3)) * 30.0
    f(xw, 30.0)  # warmup (box side must exceed 2*cutoff)
    res = {"version": ver, "cd": [], "cv": []}
    for case in ("cd", "cv"):
        for N in cfg[case + "_N"]:
            L = (N / cfg["density"]) ** (1 / 3) if case == "cd" else cfg["cv_side"]
            x = rng.random((N, 3)) * L
            t, npairs = timeit(f, x, L)
            print(f"{name} {nthreads}t {case} N={N} t={t:.3f} npairs={npairs}", flush=True)
            res[case].append([N, t, int(npairs)])
            # predict the time of the next size (linear for cd, quadratic for cv)
            p = 1 if case == "cd" else 2
            k = cfg[case + "_N"].index(N)
            if k + 1 < len(cfg[case + "_N"]) and t * (cfg[case + "_N"][k + 1] / N) ** p > TMAX:
                break
    results[name] = res
    json.dump(results, open(outfile, "w"), indent=1)
