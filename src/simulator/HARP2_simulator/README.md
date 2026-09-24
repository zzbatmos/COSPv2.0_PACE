# HARP2 simulator

The HARP2 simulator mimics the polarimetric ("cloudbow") retrieval of the liquid-cloud
droplet effective radius (re) and effective variance (ve) from the multi-angle polarized
reflectance measured by the Hyper-Angular Rainbow Polarimeter 2 (HARP2) on NASA's PACE
mission. It uses the hyper-angular 670 nm band.

**This copy is the back-port to COSP `v2.1.4cesm`**, the COSP used by the CESM2.1 and
CESM2.2 releases of CAM. The physics (`harp2_simulator.F90`), the look-up table and its
generator are identical to the COSP v2.2 version (branch `claude/great-lovelace-wsm2hk`);
only the plumbing follows the v2.1.4 layout. See "Differences from the COSP v2.2 version"
below, and `CAM_INTEGRATION.md` (section 5) for use in a CESM2.1/2.2 CAM.

Like the MODIS simulator, it works as a pseudo-retrieval. For every subcolumn it first
computes what the instrument would measure, then retrieves re and ve from that signal
with an algorithm of the same form as the real one. The vertical weighting of the
retrieved quantities therefore comes from the radiative transfer; it is not prescribed.

## Files

| File | Content |
|---|---|
| `harp2_simulator.F90` | Module `mod_harp2_sim`: forward model, cloudbow fit, grid-box statistics |
| `../cosp_harp2_interface.F90` | Module `mod_cosp_harp2_interface`: parameters, LUT reader, viewing geometry, input type |
| `../../cosp.F90` | Input assignment and daylight mask (inline, as for the other v2.1.4 simulators), calls, output fields |
| `harp2_lut_generator.py` | Generates the Mie look-up table with [miepython](https://github.com/scottprahl/miepython) |
| `harp2_lut_670nm.txt` | Look-up table used by default (670 nm) |
| `../../../unit_testing/harp2/` | Standalone unit tests (`make test`) |
| `CAM_INTEGRATION.md` | Guide for adding HARP2 to CESM/CAM (both COSP versions) |
| `cam_release_Makefile.in.patch` | Tested patch for CAM's `src/physics/cosp2/Makefile.in` (CESM2.1/2.2) |

## Method

**Forward model.** The polarized reflectance of each subcolumn is computed in the
single-scattering approximation for the stack of model layers, over a black surface:

```
Rp(j) = 1/(4(mu_j+mu0)) * sum_k [S_k(Theta_j)/dtau_k] * exp(-tau_k^top m_j) * [1 - exp(-dtau_k m_j)]
S_k   = dtauLiq_k * ssa_k * (-P12_liq(Theta; re_k, ve_k)) + dtauRay_k * (-P12_ray(Theta))
m_j   = 1/mu_j + 1/mu0
```

`-P12_liq` comes from the Mie look-up table for a gamma size distribution
(Hansen and Travis 1974). Molecular (Rayleigh) scattering is included, with the optical
depth of each layer taken from its pressure thickness. Ice crystals attenuate the signal
but are assumed to add no polarized signal. Polarization by multiple scattering is
neglected, which is a good approximation in the cloudbow region
(Breon and Goloub 1998).

**Retrieval.** The geometry-normalized signal `4(mu_j+mu0) Rp(j)` is fitted over the
scattering-angle window 135 to 165 degrees with the parametric model of Breon and Goloub
(1998):

```
A * (-P12_liq(Theta; re, ve)) + B cos^2(Theta) + C
```

A, B and C are found by linear least squares for every (re, ve) node of the look-up
table. The node with the smallest residual is kept and refined by parabolic
interpolation along re and ve. The smooth terms B and C absorb the Rayleigh signal. A
subcolumn gets one of these outcomes:

| Flag | Meaning |
|---|---|
| 0 | Clear (column optical thickness < 0.3) |
| 1 | Successful liquid cloudbow retrieval |
| 2 | Cloudy, but no cloudbow detected (A < 0.2, e.g. ice topped or thick ice above liquid) |
| 3 | Cloudy, but the cloudbow window is not sampled by the viewing geometry |
| 4 | Cloudbow detected, but the fit failed (poor fit, or re at the edge of the table) |

Two further per-subcolumn diagnostics make table saturation visible instead of silent:
**input clamped** (liquid layers with re or ve outside the table carry at least 1% of the
polarized signal; the forward model then uses the nearest table edge) and **ve at a table
limit** (a successful retrieval whose ve is 0.01 or 0.40). Retrievals at a ve limit are
kept in the means and histograms but counted separately.

Because the polarized signal is single-scattered, it comes from roughly the top
`1/m` (a few tenths) of optical depth of the cloud. The retrieved re is therefore more
cloud-top weighted than the MODIS 3.7 micron retrieval, and vertical variations of re
within that layer broaden the retrieved ve.

**Viewing geometry.** 60 along-track views between -57 and +57 degrees in the solar
principal plane, the geometry near the sub-satellite track of an early-afternoon
sun-synchronous orbit. The solar zenith angle is taken from `cospgridIN%sza` when it is
allocated (30 degrees otherwise). There are no retrievals for solar zenith angles above
75 degrees. `harp2_subcolumn` accepts any set of view angles, so other geometries (e.g.
cross-track position) can be supplied by the caller.

All retrieval and geometry parameters are set in `COSP_HARP2_INIT`. COSP v2.1.4 has no
orbit swathing, so HARP2 observes every sunlit point (`cospgridIN%sunlit > 0`) with a
solar zenith angle of at most 75 degrees.

## Inputs

| Input | Content |
|---|---|
| `cospIN%tau_067` | Subcolumn optical thickness at 0.67 microns (shared with ISCCP, MISR, MODIS) |
| `cospIN%fracLiq` | Liquid fraction of the optical thickness (shared with MODIS) |
| `cospIN%reffLiq` | Subcolumn liquid effective radius (microns) |
| `cospIN%veffLiq` | Subcolumn liquid effective variance |
| `cospgridIN%phalf`, `cospgridIN%sunlit`, `cospgridIN%sza` | Pressure at layer edges, daylight flag, solar zenith angle |

`cospgridIN%sza` (degrees) is a new, optional component of `cosp_column_inputs` added by
this back-port (v2.1.4 had none). Only HARP2 uses it. Hosts should allocate it and fill it
with the actual solar zenith angle, e.g. in CAM `acos(coszrs)*180/pi`; if it is not
allocated, HARP2 assumes 30 degrees. The offline driver sets it to 0, as the COSP v2.2
driver does.

The example optics in the offline driver use the model effective radius and a constant
effective variance (`harp2_veffLiq` in the input namelist, default 0.10). Radii and
variances outside the table are clamped to its edges in the forward model.

## Outputs

| cospOUT field | netCDF name | Units | Content |
|---|---|---|---|
| `harp2_Cloud_Fraction_Liquid_Mean` | `clwharp2` | % | Fraction of subcolumns with a successful cloudbow retrieval |
| `harp2_Cloud_Particle_Size_Liquid_Mean` | `reffclwharp2` | m | Mean retrieved effective radius |
| `harp2_Effective_Variance_Liquid_Mean` | `veffclwharp2` | 1 | Mean retrieved effective variance |
| `harp2_Reff_vs_Veff_Liquid` | `clharp2reffveff` | % | Joint histogram of re and ve |
| `harp2_Retrieval_Flag_Fraction` | `harp2_flag_fraction` | % | Fraction of subcolumns per outcome flag 0 to 4 (sums to 100) |
| `harp2_Veff_Limit_Fraction` | `harp2_veff_limit_fraction` | % | Subcolumns whose successful retrieval has ve at a table limit |
| `harp2_Input_Clamped_Fraction` | `harp2_input_clamped_fraction` | % | Subcolumns whose liquid re or ve lie outside the table |

`clwharp2` is the fraction of subcolumns with a successful retrieval, not the liquid cloud
cover. The histogram bins are explicit (`mod_cosp_config`, independent of the v2.1.4
MODIS bins): re edges 0, 4, 8, 10, 12.5, 15,
20, 30 microns and above; ve edges 0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.125, 0.15, 0.175, 0.20,
0.25, 0.30, 0.35, 1.0. All outputs are `R_UNDEF` at night and for solar zenith angles
above 75 degrees; on observed
points without any retrieval, the histogram and fractions are 0 and the mean re and ve
are `R_UNDEF`. Hosts that time-average cloud-fraction-weighted means must use the same
missing-value mask for numerator and denominator (see `CAM_INTEGRATION.md`, 3.5).

## Using it

Call `COSP_INIT` with the two optional arguments

```fortran
call COSP_INIT(..., Lharp2=.true., harp2_lut_file='path/to/harp2_lut_670nm.txt')
if (.not. harp2_lut_loaded()) stop 'HARP2 look-up table not loaded'   ! or the host's abort
```

and associate any of the output fields above. The two arguments are optional and
trailing, so existing positional `COSP_INIT` calls (such as CAM's) compile unchanged.
Without the check, a missing table only prints an error, and COSP continues with the
HARP2 outputs set to `R_UNDEF`. Because `cosp.F90` now uses the HARP2 modules, a host
build must compile `harp2_simulator.F90` and `cosp_harp2_interface.F90` even when HARP2
is switched off (for CAM: `cam_release_Makefile.in.patch`).

In the offline driver, set the flags `Lclwharp2`, `Lreffclwharp2`, `Lveffclwharp2`,
`Lclharp2reffveff`, `Lharp2flagfrac`, `Lharp2vefflimfrac` and `Lharp2clampfrac` in
`driver/run/cosp2_output_nl.txt`, and `harp2_lut_file` (and optionally `harp2_veffLiq`) in
`driver/run/cosp2_input_nl.txt`. The flags are off by default because the reference
outputs do not contain HARP2 fields.

## Differences from the COSP v2.2 version

| Item | COSP v2.2 branch | This back-port (v2.1.4cesm) |
|---|---|---|
| `harp2_simulator.F90`, LUT, generator, unit test | reference | identical files |
| Input and output types | `mod_cosp_stats` (`cosp_stats.F90`) | `mod_cosp` (`cosp.F90`) |
| `cospIN%reffLiq`, `cospIN%veffLiq` | new | new (same names, units) |
| `cospgridIN%sza` | existing field (shared with RTTOV) | new field, HARP2 only |
| Input assignment | `COSP_ASSIGN_harp2IN` in the interface module | inline in `cosp.F90` |
| Orbit swaths | `cospIN%harp2_swathIN` | none (v2.1.4 has no swathing) |
| `COSP_INIT` | trailing optional `Lharp2`, `harp2_lut_file` | the same |
| Output pointers, netCDF names, bins | reference | identical |

Verification of the back-port (COSP offline driver, UKMO test case, gfortran):

* With HARP2 off, all 135 output variables are bit-identical to unmodified
  `v2.1.4cesm`. With HARP2 on, all pre-existing variables are identical except
  `calipso_tau`, which v2.1.4 allocates but never fills (uninitialized memory; one value
  of 1.9e-40 appeared).
* With HARP2 on, every HARP2 output is bit-identical to the COSP v2.2 branch, unchunked
  and with `NPOINTS_IT=40` (4 chunks).
* A missing look-up table disables HARP2 (outputs `R_UNDEF`) and leaves the other outputs
  unchanged.
* The unit tests pass in single and double precision.
* The COSP library builds with CAM release's `Makefile.in` (patched), CAM's
  `cosp_kinds.F90` (double precision), `cosp_errorHandling.F90` and optics, also with
  `-fdefault-real-8`; without the patch the build fails, as expected.
* Pre-existing v2.1.4 quirk, unrelated to HARP2: in chunked offline runs, `cloudsatpia`
  holds uninitialized values for points outside the first chunk (the whole array, not the
  chunk section, is passed to `quickbeam_column`; fixed in COSP v2.2). CAM calls COSP with
  one chunk (`start_idx=1, stop_idx=ncol`), so CAM is not affected.

## Regenerating the look-up table

```bash
pip install miepython numba scipy
MIEPYTHON_USE_JIT=1 python3 harp2_lut_generator.py -o harp2_lut_670nm.txt
```

The defaults are those of the distributed table: 670 nm, m = 1.331 - 1.9e-8i,
scattering angles 125 to 170 degrees every 0.25 degrees, re from 4 to 30 microns every
0.5 microns, 20 values of ve from 0.01 to 0.40 (denser above 0.1, where broad model
distributions such as CAM6/MG2's lie), a size-parameter step of 0.02 and radii up to
300 microns. Accuracy of -P12, as a fraction of its peak: about 0.2% from the size
integration for the narrowest distributions and 0.1% at ve = 0.40; at most 0.6% from
linear interpolation in scattering angle and 0.08% from interpolation between ve nodes.
The generator computes the neglected fraction of geometric cross section outside the
radius range analytically (incomplete gamma function) and refuses to run if it exceeds
1e-6 (it is 2e-8 for the distributed table). Other bands (e.g. 870 nm) can be produced
with `--wavelength`, `--m-real` and `--m-imag`. Generation takes about 11 minutes.

## Limitations and possible extensions

* Single scattering only: no surface, aerosol, or multiply scattered polarized signal,
  and no 3D radiative effects.
* No measurement noise and no spatial aggregation: real cloudbow retrievals often
  aggregate several pixels, which also broadens ve through horizontal variability of re.
* Ice is non-polarizing in the forward model.
* Fixed principal-plane view geometry (no cross-track position).
* Not yet validated against vector radiative transfer or the real HARP2 product; treat
  the output as an idealized pseudo-retrieval.
* Cloud optical thickness and cloud-top height are not retrieved; the MODIS and ISCCP
  simulators provide those.

## References

* Breon, F.-M., and P. Goloub, 1998: Cloud droplet effective radius from spaceborne
  polarization measurements. *Geophys. Res. Lett.*, 25, 1879-1882.
* Hansen, J. E., and L. D. Travis, 1974: Light scattering in planetary atmospheres.
  *Space Sci. Rev.*, 16, 527-610.
