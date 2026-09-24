#!/usr/bin/env python3
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Copyright (c) 2026, University of Maryland Baltimore County
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are
# permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of
#    conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list
#    of conditions and the following disclaimer in the documentation and/or other
#    materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be
#    used to endorse or promote products derived from this software without specific prior
#    written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
# THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
# OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# History
# Sep 2026 - Initial version, HARP2 simulator for COSP2
# Sep 2026 - ve range extended to 0.40 (broad CAM6/MG2 distributions); tail check
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
Generate the look-up table (LUT) of the polarized phase function -P12(Theta; re, ve) of
liquid water clouds used by the COSP HARP2 simulator (mod_harp2_sim).

The droplet size distribution is the two-parameter gamma distribution of Hansen and
Travis (1974, eq. 2.56),

    n(r) ~ r**((1 - 3 ve)/ve) * exp(-r / (re ve)),

where re is the effective radius and ve the effective variance. Mie calculations use
miepython (https://github.com/scottprahl/miepython). Set MIEPYTHON_USE_JIT=1 (requires
numba) for a large speed-up. scipy is used to check, before the Mie calculations, that the
radius range covers every size distribution of the table (tail fraction < 1e-6).

Normalization: the phase matrix is normalized so that (1/4pi) * integral(P11 dOmega) = 1,
so that the single-scattering polarized reflectance of a layer is

    Rp = ssa * (-P12(Theta)) / (4 (mu + mu0)) * [1 - exp(-tau (1/mu + 1/mu0))].

Sign convention: P12 = (|S2|^2 - |S1|^2)/2 (Bohren and Huffman 1983), so -P12 > 0 means
polarization perpendicular to the scattering plane (as for Rayleigh scattering and for the
primary cloudbow).

Output is a plain-text file read by cosp_harp2_init (no netCDF dependency in COSP):
    lines beginning with '#' are comments
    nTheta nRe nVe
    wavelength (microns)
    theta(1:nTheta)                     scattering angle (degrees)
    re(1:nRe)                           effective radius (microns)
    ve(1:nVe)                           effective variance (-)
    ssa(1:nRe,1:nVe)                    single-scattering albedo, re varies fastest
    qext(1:nRe,1:nVe)                   mean extinction efficiency, re varies fastest
    mP12(1:nTheta,1:nRe,1:nVe)          -P12, theta varies fastest, then re, then ve

Example:
    MIEPYTHON_USE_JIT=1 python3 harp2_lut_generator.py -o harp2_lut_670nm.txt
"""
import argparse
import datetime

import numpy as np
import miepython

DEFAULT_VE = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18,
              0.20, 0.22, 0.24, 0.27, 0.30, 0.35, 0.40]
MAX_TAIL_FRACTION = 1.0e-6   # largest acceptable neglected fraction of scattering cross section


def truncated_fraction(re, ve, r_min, r_max):
    """Fraction of the geometric cross section outside [r_min, r_max].

    The area-weighted gamma distribution r**2 n(r) is itself a gamma distribution with
    shape 1/ve and scale re*ve, so the neglected fractions are regularized incomplete
    gamma functions. This bounds the relative error in the scattering integrals.
    """
    from scipy.special import gammainc, gammaincc
    k, scale = 1.0 / ve, re * ve
    return gammainc(k, r_min / scale) + gammaincc(k, r_max / scale)


def gamma_number_distribution(r, re, ve):
    """Unnormalized gamma size distribution n(r) of Hansen and Travis (1974)."""
    lnn = ((1.0 - 3.0 * ve) / ve) * np.log(r) - r / (re * ve)
    return np.exp(lnn - lnn.max())


def monodisperse_tables(m, wavelength, theta, dx, r_max):
    """Mie S11, S12 (qsca normalization), Qext and Qsca on a uniform size-parameter grid."""
    x = np.arange(dx, 2.0 * np.pi * r_max / wavelength + dx, dx)
    mu = np.cos(np.radians(theta))
    s11 = np.empty((x.size, theta.size))
    s12 = np.empty((x.size, theta.size))
    qext = np.empty(x.size)
    qsca = np.empty(x.size)
    for i, xi in enumerate(x):
        S1, S2 = miepython.S1_S2(m, xi, mu, norm='qsca')
        a1, a2 = np.abs(S1) ** 2, np.abs(S2) ** 2
        s11[i] = 0.5 * (a1 + a2)
        s12[i] = 0.5 * (a2 - a1)
        qext[i], qsca[i], _, _ = miepython.efficiencies_mx(m, xi)
    r = x * wavelength / (2.0 * np.pi)
    return r, s11, s12, qext, qsca


def polydisperse(r, s11, s12, qext, qsca, re, ve):
    """Size-distribution averaged ssa, qext and 4pi-normalized P11, P12."""
    w = gamma_number_distribution(r, re, ve) * np.pi * r ** 2   # area weighting
    csca = np.dot(w, qsca)
    cext = np.dot(w, qext)
    p11 = 4.0 * np.pi * np.dot(w, s11) / csca
    p12 = 4.0 * np.pi * np.dot(w, s12) / csca
    return csca / cext, cext / w.sum(), p11, p12


def build_lut(wavelength, m, theta, re_grid, ve_grid, dx, r_max):
    r, s11, s12, qe, qs = monodisperse_tables(m, wavelength, theta, dx, r_max)
    ssa = np.empty((re_grid.size, ve_grid.size))
    qext = np.empty_like(ssa)
    mp12 = np.empty((theta.size, re_grid.size, ve_grid.size))
    for j, ve in enumerate(ve_grid):
        for i, re in enumerate(re_grid):
            ssa[i, j], qext[i, j], _, p12 = polydisperse(r, s11, s12, qe, qs, re, ve)
            mp12[:, i, j] = -p12
    return ssa, qext, mp12


def write_lut(fname, wavelength, m, theta, re_grid, ve_grid, ssa, qext, mp12, dx, r_max, tail):
    def block(values, per_line=8, fmt='{:.4e}'):
        values = np.ravel(values, order='F')
        return '\n'.join(' '.join(fmt.format(v) for v in values[k:k + per_line])
                         for k in range(0, values.size, per_line)) + '\n'

    with open(fname, 'w') as f:
        f.write('# COSP HARP2 simulator: liquid-cloud polarized phase function LUT, -P12(Theta; re, ve)\n')
        f.write('# Generated {} by harp2_lut_generator.py with miepython {}\n'.format(
            datetime.date.today().isoformat(), miepython.__version__))
        f.write('# Gamma size distribution (Hansen and Travis 1974); refractive index m = {:.4f} - {:.3e}i\n'
                .format(m.real, -m.imag))
        f.write('# Size-parameter step dx = {}, maximum radius = {} microns, largest neglected\n'
                '# fraction of cross section = {:.1e}\n'.format(dx, r_max, tail))
        f.write('# P11 normalized to (1/4pi) int(P11 dOmega) = 1; -P12 > 0 is perpendicular polarization\n')
        f.write('# Layout: nTheta nRe nVe / wavelength (um) / theta (deg) / re (um) / ve / ssa(re,ve) /\n')
        f.write('#         qext(re,ve) / -P12(theta,re,ve); first index varies fastest\n')
        f.write('{} {} {}\n'.format(theta.size, re_grid.size, ve_grid.size))
        f.write('{:.4f}\n'.format(wavelength))
        f.write(block(theta, fmt='{:.3f}'))
        f.write(block(re_grid, fmt='{:.3f}'))
        f.write(block(ve_grid, fmt='{:.4f}'))
        f.write(block(ssa, fmt='{:.8f}'))
        f.write(block(qext, fmt='{:.5f}'))
        f.write(block(mp12))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-o', '--output', default='harp2_lut_670nm.txt')
    p.add_argument('--wavelength', type=float, default=0.670, help='microns')
    p.add_argument('--m-real', type=float, default=1.331, help='real refractive index of water')
    p.add_argument('--m-imag', type=float, default=1.9e-8, help='imaginary refractive index (positive)')
    p.add_argument('--theta-min', type=float, default=125.0)
    p.add_argument('--theta-max', type=float, default=170.0)
    p.add_argument('--dtheta', type=float, default=0.25)
    p.add_argument('--re-min', type=float, default=4.0)
    p.add_argument('--re-max', type=float, default=30.0)
    p.add_argument('--dre', type=float, default=0.5)
    p.add_argument('--ve', type=float, nargs='+', default=DEFAULT_VE)
    p.add_argument('--dx', type=float, default=0.02, help='size-parameter step of the Mie integration')
    p.add_argument('--r-max', type=float, default=300.0, help='largest droplet radius (microns)')
    a = p.parse_args()

    m = complex(a.m_real, -a.m_imag)   # miepython convention: m = n - ik
    theta = np.arange(a.theta_min, a.theta_max + 0.5 * a.dtheta, a.dtheta)
    re_grid = np.arange(a.re_min, a.re_max + 0.5 * a.dre, a.dre)
    ve_grid = np.array(a.ve)
    if np.any(np.diff(ve_grid) <= 0) or ve_grid.min() <= 0 or ve_grid.max() >= 0.5:
        raise SystemExit('ve must be increasing and within (0, 0.5)')

    # Check, before the expensive Mie calculations, that the radius range covers every
    # size distribution of the table
    r_min = a.dx * a.wavelength / (2.0 * np.pi)
    tail = max(truncated_fraction(re, ve, r_min, a.r_max) for re in re_grid for ve in ve_grid)
    print('Largest neglected fraction of cross section: {:.2e}'.format(tail))
    if tail > MAX_TAIL_FRACTION:
        raise SystemExit('Increase --r-max: neglected fraction {:.2e} > {:.1e}'.format(
            tail, MAX_TAIL_FRACTION))

    ssa, qext, mp12 = build_lut(a.wavelength, m, theta, re_grid, ve_grid, a.dx, a.r_max)
    write_lut(a.output, a.wavelength, m, theta, re_grid, ve_grid, ssa, qext, mp12, a.dx, a.r_max,
              tail)
    print('Wrote {} ({} x {} x {})'.format(a.output, theta.size, re_grid.size, ve_grid.size))


if __name__ == '__main__':
    main()
