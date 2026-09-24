#!/usr/bin/env python3
"""Plot -P12 from the HARP2 look-up table and check its numerical convergence.

Figure 1 (p12_lut_overview.png): -P12(Theta) for several effective radii (CER) and
effective variances (CEV). The oscillations after the main cloudbow are the
supernumerary bows; they fade as CEV increases.

Figure 2 (p12_lut_convergence.png, only with --check): for three size distributions,
the LUT compared with an independent integration with a 4x finer size-parameter step
(dx = 0.005) on a 0.05 degree grid, and with a deliberately under-averaged integration
(dx = 0.08). Residual interference ripple from an insufficient size average shows up as
fast oscillations (period about 1 to 2 degrees) in the difference panels. The check
takes about 10 minutes with MIEPYTHON_USE_JIT=1.

Example:
    MIEPYTHON_USE_JIT=1 python3 plot_harp2_lut.py --check
"""
import argparse
import os
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_lut(fname):
    """Return wavelength (um), theta (deg), re (um), ve and -P12(theta, re, ve)."""
    vals = []
    with open(fname) as f:
        for line in f:
            if not line.startswith('#'):
                vals += line.split()
    nt, nr, nv = map(int, vals[:3])
    v = np.array(vals[3:], float)
    wl, o = v[0], 1
    theta = v[o:o + nt]; o += nt
    re = v[o:o + nr]; o += nr
    ve = v[o:o + nv]; o += nv
    o += 2 * nr * nv                                        # ssa, qext
    if o + nt * nr * nv != v.size:
        raise ValueError(f'{fname}: unexpected number of values')
    mp12 = v[o:].reshape(nv, nr, nt).transpose(2, 1, 0)     # theta varies fastest
    return wl, theta, re, ve, mp12


def lut_curve(re, ve, mp12, r, e):
    i, j = np.argmin(abs(re - r)), np.argmin(abs(ve - e))
    if abs(re[i] - r) > 1e-6 or abs(ve[j] - e) > 1e-6:
        raise ValueError(f'CER {r}, CEV {e} is not a LUT node')
    return mp12[:, i, j]


def overview(lut, fname):
    wl, theta, re, ve, mp12 = lut
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.5), sharex=True)

    def panel(a, pairs, labels, title):
        for k, ((r, e), lab) in enumerate(zip(pairs, labels)):
            a.plot(theta, lut_curve(re, ve, mp12, r, e),
                   color=plt.cm.viridis(0.9 * k / (len(pairs) - 1)), lw=1.4, label=lab)
        a.axvspan(135, 165, color='0.92', zorder=0)
        a.set_title(title, fontsize=11)
        a.grid(alpha=0.3)
        a.legend(fontsize=8.5, ncol=2)
        a.set_ylabel(r'$-P_{12}$')

    radii = [5, 8, 12, 16, 20, 25]
    variances = [0.01, 0.02, 0.05, 0.10, 0.20, 0.40]
    panel(ax[0, 0], [(r, 0.01) for r in radii], [f'CER {r} µm' for r in radii],
          'CEV = 0.01 (narrowest in the LUT)')
    panel(ax[0, 1], [(r, 0.05) for r in radii], [f'CER {r} µm' for r in radii], 'CEV = 0.05')
    panel(ax[1, 0], [(10, e) for e in variances], [f'CEV {e:.2f}' for e in variances],
          'CER = 10 µm')
    panel(ax[1, 1], [(20, e) for e in variances], [f'CEV {e:.2f}' for e in variances],
          'CER = 20 µm')
    for a in ax[1]:
        a.set_xlabel('Scattering angle (deg)')
    fig.suptitle(f'HARP2 LUT: $-P_{{12}}$ at {wl * 1000:.0f} nm, gamma size distribution '
                 '(grey band = 135-165° fit window)', fontsize=12)
    fig.tight_layout()
    fig.savefig(fname, dpi=130)
    plt.close(fig)


def convergence(lut, fname, m):
    import miepython
    from scipy.special import gammaincinv
    import harp2_lut_generator as gen

    wl, theta, re, ve, mp12 = lut
    th_f = np.arange(theta[0], theta[-1] + 1e-4, 0.05)
    mu_f = np.cos(np.radians(th_f))

    def integrate(r_eff, v_eff, dx):
        # The area-weighted distribution is gamma(1/ve, re*ve): integrate over all but
        # 2e-10 of the geometric cross section.
        k, sc = 1.0 / v_eff, r_eff * v_eff
        rlo, rhi = gammaincinv(k, 1e-10) * sc, gammaincinv(k, 1 - 1e-10) * sc
        x = np.arange(2 * np.pi * rlo / wl, 2 * np.pi * rhi / wl + dx, dx)
        r = x * wl / (2 * np.pi)
        w = gen.gamma_number_distribution(r, r_eff, v_eff) * np.pi * r ** 2
        s12, csca = np.zeros(th_f.size), 0.0
        for wi, xi in zip(w, x):
            s1, s2 = miepython.S1_S2(m, xi, mu_f, norm='qsca')
            s12 += wi * 0.5 * (abs(s2) ** 2 - abs(s1) ** 2)
            csca += wi * miepython.efficiencies_mx(m, xi)[1]
        return -4 * np.pi * s12 / csca

    def single_droplet(r):
        s1, s2 = miepython.S1_S2(m, 2 * np.pi * r / wl, mu_f, norm='qsca')
        return -4 * np.pi * 0.5 * (abs(s2) ** 2 - abs(s1) ** 2)

    cases = [(10, 0.01), (25, 0.01), (10, 0.10)]
    fig, ax = plt.subplots(3, 2, figsize=(13, 11), gridspec_kw={'width_ratios': [1.6, 1]})
    for row, (r, e) in enumerate(cases):
        ref, coarse = integrate(r, e, 0.005), integrate(r, e, 0.08)
        tab = lut_curve(re, ve, mp12, r, e)
        pk = ref.max()
        a = ax[row, 0]
        if row == 0:
            a.plot(th_f, single_droplet(r), color='0.75', lw=0.6,
                   label=f'single droplet r = {r} µm (no averaging)')
        a.plot(th_f, coarse, color='tab:red', lw=0.8, label='size step dx = 0.08 (under-averaged)')
        a.plot(th_f, ref, color='k', lw=1.6, label='reference: dx = 0.005, 0.05° grid')
        a.plot(theta, tab, 'o', color='tab:blue', ms=2.2, label='LUT: dx = 0.02, 0.25° grid')
        a.axvspan(135, 165, color='0.93', zorder=0)
        if row == 0:
            a.set_ylim(min(ref.min(), 0) - 0.3 * pk, 1.6 * pk)
        a.set_title(f'CER = {r} µm, CEV = {e}', fontsize=11)
        a.grid(alpha=0.3)
        a.set_ylabel(r'$-P_{12}$')
        a.legend(fontsize=8)

        b = ax[row, 1]
        ref_i, coarse_i = np.interp(theta, th_f, ref), np.interp(theta, th_f, coarse)
        b.plot(theta, 100 * (coarse_i - ref_i) / pk, color='tab:red', lw=0.8,
               label='dx = 0.08 minus reference')
        b.plot(theta, 100 * (tab - ref_i) / pk, color='tab:blue', lw=1.4, label='LUT minus reference')
        b.axhline(0, color='k', lw=0.5)
        b.grid(alpha=0.3)
        b.set_ylabel('Difference (% of peak)')
        b.legend(fontsize=8)
        err_lut = 100 * np.abs(tab - ref_i).max() / pk
        err_coarse = 100 * np.abs(coarse_i - ref_i).max() / pk
        err_interp = 100 * np.abs(np.interp(th_f, theta, tab) - ref).max() / pk
        b.set_title(f'max |LUT - ref| = {err_lut:.2f}%,  dx = 0.08: {err_coarse:.1f}%', fontsize=10)
        print(f'CER {r} CEV {e}: max |LUT - ref| = {err_lut:.3f}% of peak, dx = 0.08: '
              f'{err_coarse:.2f}%, 0.25 deg linear interpolation: {err_interp:.3f}%', flush=True)
    for a in ax[-1]:
        a.set_xlabel('Scattering angle (deg)')
    fig.suptitle('Is the size averaging sufficient?  LUT versus an independent 4x finer size '
                 'integration', fontsize=12)
    fig.tight_layout()
    fig.savefig(fname, dpi=130)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--lut', default=os.path.join(HERE, 'harp2_lut_670nm.txt'))
    p.add_argument('--outdir', default=os.path.join(HERE, 'figures'))
    p.add_argument('--check', action='store_true', help='also run the convergence check (slow)')
    p.add_argument('--m-real', type=float, default=1.331, help='must match the LUT')
    p.add_argument('--m-imag', type=float, default=1.9e-8, help='must match the LUT')
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    lut = read_lut(args.lut)
    overview(lut, os.path.join(args.outdir, 'p12_lut_overview.png'))
    if args.check:
        convergence(lut, os.path.join(args.outdir, 'p12_lut_convergence.png'),
                    complex(args.m_real, -args.m_imag))


if __name__ == '__main__':
    main()
