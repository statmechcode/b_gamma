#!/usr/bin/env python3
"""
Author: Giuseppe Petrillo (giuseppe51289@gmail.com)
"""

import numpy as np
import math
from scipy.optimize import minimize

# Global Parameters
BINP = 1.0
AINP = 0.5
DMTH = 0.0
ITER_BOOT = 1000
EPS = 1e-15
R_EARTH = 6370.0  # km


def geodist(lat1, lon1, lat2, lon2):
    """Geodetic distance between two points in km"""
    # converti in radianti
    phi1, lam1 = math.radians(lat1), math.radians(lon1)
    phi2, lam2 = math.radians(lat2), math.radians(lon2)
    dlam = abs(lam1 - lam2)
    cos_d = (math.sin(phi1)*math.sin(phi2)
             + math.cos(phi1)*math.cos(phi2)*math.cos(dlam))
    # correggi errore numerico
    cos_d = min(1.0, max(-1.0, cos_d))
    return math.acos(cos_d) * R_EARTH


def neg_log_likelihood(params, dm):
    """Negative log likelihood"""
    b, a = params
    b = max(b, EPS)
    a = max(a, EPS)
    gamma = b / (b + a)
    dmm_mean = np.mean(np.abs(dm))
    # somma del termine log
    sum_term = np.sum(np.log(np.maximum(1.0 - gamma * np.exp(-a * np.abs(dm)), EPS)))
    logL = (
        math.log(b) - math.log(2.0) - math.log(a)
        + 2.0 * math.log(a + b) - math.log(2.0 * b + a)
        - b * dmm_mean + sum_term / len(dm)
    )
    return -logL


def fit_likelihood(dm):
    res = minimize(
        neg_log_likelihood,
        x0=[BINP, AINP],
        args=(dm,),
        bounds=[(0.8, 4.0), (0.3, 10.0)],
        method='L-BFGS-B'
    )
    return res.x


def boot(npt, dm, iter_boot=ITER_BOOT):
    bs, as_ = [], []
    for _ in range(iter_boot):
        sample = np.random.choice(dm, size=npt, replace=True)
        bb, aa = fit_likelihood(sample)
        bs.append(bb)
        as_.append(aa)
    return (
        np.mean(bs), np.std(bs, ddof=1),
        np.mean(as_), np.std(as_, ddof=1)
    )


def main():
    # Reading catalog
    t, q, y, x, z = [], [], [], [], []
    with open('input_catalog.dat') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            ti, qi, yi, xi, zi = map(float, parts[:5])
            t.append(ti)
            x.append(xi)
            y.append(yi)
            q.append(qi + 0.1 * (np.random.rand() - 0.5))

    x = np.array(x)
    y = np.array(y)
    q = np.array(q)
    npt = len(q)
    print(f"Number of points: {npt}")

    # output file
    with open('output.dat', 'w') as out:
        dth = 20.0 
        print(f"\n--- dth = {dth:.2f} km ---")

        # dm calculus
        dm_list = []
        for i in range(npt - 100):
            for j in range(1, 101):
                k = i + j
                if k >= npt:
                    break
                d = geodist(y[i], x[i], y[k], x[k])
                if d < dth:
                    dmm = q[k] - q[i]
                    if dmm > DMTH:
                        dm_list.append(dmm - DMTH)
                    break  # esci dal loop j
            # fine loop j
        k_pairs = len(dm_list)
        if k_pairs == 0:
            print(f"No valid couples for dth={dth:.2f}")
            return

        print(f"Find {k_pairs} couples, executing fit likelihood...")
        b_opt, a_opt = fit_likelihood(dm_list)
        print(f"b = {b_opt:.4f}, a = {a_opt:.4f}")

        print("Executing bootstrapping...")
        bm, sb, am, sa = boot(k_pairs, dm_list)
        print(f"Final result: b = {bm:.4f} ± {sb:.4f}, a = {am:.4f} ± {sa:.4f}")

        out.write(f"{dth} {bm:.6f} {sb:.6f} {am:.6f} {sa:.6f}\n")
        out.flush()


if __name__ == '__main__':
    main()
