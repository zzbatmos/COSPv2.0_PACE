# Integrating the COSP HARP2 simulator into CESM/CAM: version guide

Audience: an agent (or developer) working in a CESM/CAM checkout that wants to add the
HARP2 simulator from this repository. Read section 2 first: the right integration path
depends on which COSP version your CAM uses, and the two paths are not interchangeable.

Revision 3. New: the back-port for CESM2.1/2.2 exists and was verified (branch
`claude/harp2-cosp-v2.1.4cesm`); Path A (section 5) now uses it, with a tested patch for
CAM's `Makefile.in` and release line numbers for every CAM edit. This file is identical on
both HARP2 branches.

Revision 2 (after an independent review from a CESM2.1 installation). Changes from
revision 1: corrected history-averaging recipe (3.5), corrected MG2 `MU` facts (3.2),
the ve range of the LUT extended to 0.40 with saturation diagnostics, HARP2 swath input
made non-breaking (no CAM change needed for it any more), explicit histogram bins,
`harp2_lut_loaded()`, new outcome diagnostics, validation wording (6).

All facts below were checked against these commits (line numbers refer to them; `grep`
in your own checkout before editing, since line numbers drift):

| Code | Branch / tag | Commit |
|---|---|---|
| HARP2 COSP branch (v2.2) | `zzbatmos/COSPv2.0_PACE`, branch `claude/great-lovelace-wsm2hk` | base `5eb05e5` (= CFMIP COSP `v2.2.1`) + HARP2 commits |
| HARP2 back-port (v2.1.4) | `zzbatmos/COSPv2.0_PACE`, branch `claude/harp2-cosp-v2.1.4cesm` | base `34d8eef` (= CFMIP COSP `v2.1.4cesm`) + HARP2 commit |
| CAM for CESM2.1 | `ESCOMP/CAM` `cam_cesm2_1_rel` (branch head) | `405a3f2` |
| CAM for CESM2.2 | `ESCOMP/CAM` `cam_cesm2_2_rel` | `c47c4d7` |
| CAM development | `ESCOMP/CAM` `cam_development` | `e9e28e4` (2026-09-18) |
| COSP used by CESM2.1/2.2 | `CFMIP/COSPv2.0` tag `v2.1.4cesm` | `34d8eef` (2019-10-15) |

The reviewed installation runs `cam_cesm2_1_rel_60` (commit `a03b84b`) with COSP
`v2.1.4cesm`; its line numbers can differ by a few lines from `405a3f2` (e.g. the `MU`
reset is at line 2551 there and 2553 here).

## 1. What the HARP2 branch provides

* Physics: `src/simulator/HARP2_simulator/harp2_simulator.F90` (module `mod_harp2_sim`).
  Single-scattering polarized reflectance at 670 nm for 60 HARP2 view angles, then a
  Breon and Goloub (1998) cloudbow fit `A*(-P12(re,ve)) + B cos^2 + C` over 135 to 165
  degrees. See `README.md` in the same directory for the method.
* Interface: `src/simulator/cosp_harp2_interface.F90` (module `mod_cosp_harp2_interface`):
  parameters, LUT reader, viewing geometry, `COSP_ASSIGN_harp2IN` (daylight and swath
  masks), and `harp2_lut_loaded()`.
* Look-up table: `src/simulator/HARP2_simulator/harp2_lut_670nm.txt` (text, generated with
  miepython by `harp2_lut_generator.py`): 181 scattering angles (125 to 170 degrees),
  53 radii (4 to 30 microns), **20 variances (0.01 to 0.40)**. Numerical accuracy (section
  3.2) was checked at the new broad end.
* New COSP inputs (in `cosp_optical_inputs`, v2.2.1 layout):
  `cospIN%reffLiq(npoints,ncolumns,nlevels)` liquid effective radius in **microns**,
  `cospIN%veffLiq(npoints,ncolumns,nlevels)` liquid effective variance.
  Also used: `cospIN%tau_067`, `cospIN%fracLiq`, `cospgridIN%phalf`, `cospgridIN%sunlit`,
  `cospgridIN%sza` (degrees; defaults to 30 if not allocated).
* HARP2 orbit swaths come from a **separate** component, `cospIN%harp2_swathIN`
  (`type(swath_inputs)`, off by default). `cospIN%cospswathsIN` keeps its original
  `dimension(6)`, so v2.2 hosts that copy a six-element array compile unchanged.
* New `cosp_outputs` pointers (all `R_UNDEF` where HARP2 does not observe: night,
  SZA > 75 degrees, outside the HARP2 swath):

  | Pointer | Units | Content |
  |---|---|---|
  | `harp2_Cloud_Fraction_Liquid_Mean` | % | fraction of subcolumns with a successful cloudbow retrieval (**not** liquid cloud cover) |
  | `harp2_Cloud_Particle_Size_Liquid_Mean` | m | mean retrieved re over successful retrievals (`R_UNDEF` if none) |
  | `harp2_Effective_Variance_Liquid_Mean` | 1 | mean retrieved ve over successful retrievals (`R_UNDEF` if none) |
  | `harp2_Reff_vs_Veff_Liquid(npoints,numHARP2ReffBins,numHARP2VeffBins)` | % | joint re-ve histogram; sums to the retrieval fraction |
  | `harp2_Retrieval_Flag_Fraction(npoints,numHARP2Flags)` | % | fraction of subcolumns per outcome, index = flag + 1: 0 clear (tau < 0.3), 1 retrieval, 2 no cloudbow, 3 cloudbow not sampled by the geometry, 4 fit failed; sums to 100 |
  | `harp2_Veff_Limit_Fraction` | % | subcolumns whose successful retrieval has ve on a table limit (0.01 or 0.40) |
  | `harp2_Input_Clamped_Fraction` | % | subcolumns where liquid layers with (re, ve) outside the table carry at least 1% of the polarized signal |

* Histogram bins (in `mod_cosp_config`, explicit, identical in any COSP version):
  re edges 0, 4, 8, 10, 12.5, 15, 20, 30, 10^4 microns (8 bins, the COSP v2.2 MODIS liquid
  edges, stored in meters); ve edges 0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.125, 0.15, 0.175,
  0.20, 0.25, 0.30, 0.35, 1.0 (13 bins; the last bin holds ve at the 0.40 limit).
* `COSP_INIT` gained two trailing optional arguments: `Lharp2`, `harp2_lut_file`.
  Existing positional calls are unaffected.
* Tested only in the COSP offline driver (unit tests pass; with HARP2 off or on, every
  pre-existing output is bit-identical to COSP `v2.2.1`, also with chunking, model levels
  and swathing). The reviewer also ran the unit tests with CAM's double-precision
  `cosp_kinds`. Not yet built or run inside CAM.
* The **back-port** branch `claude/harp2-cosp-v2.1.4cesm` provides the same simulator for
  COSP `v2.1.4cesm`: identical physics files, LUT, output pointers and bins; types in
  `mod_cosp`; `cospgridIN%sza` added as a new optional field; no swaths. Its offline
  HARP2 outputs are bit-identical to the v2.2 branch, and its COSP library builds with
  CAM release's own `Makefile.in` (patched), `cosp_kinds.F90` and optics (section 5).

## 2. Step 1: find out which COSP your CAM uses

Run from the CESM root (CESM2.x places CAM in `components/cam`; adjust if you use a
standalone CAM checkout):

```bash
CAM=components/cam
git -C $CAM describe --tags --always; git -C $CAM rev-parse --abbrev-ref HEAD
grep -A7 '\[cosp2\]' $CAM/Externals_CAM.cfg 2>/dev/null        # CESM2.1/2.2 (manage_externals)
grep -A7 'submodule "cosp2"' $CAM/.gitmodules 2>/dev/null      # cam_development (git-fleximod)
C=$CAM/src/physics/cosp2/src                                   # the COSP checkout itself
git -C $C describe --tags --always 2>/dev/null
grep -l "type cosp_optical_inputs" $C/src/cosp.F90 $C/src/cosp_stats.F90
grep -c "type swath_inputs" $C/src/cosp_stats.F90
grep -n "SUBROUTINE COSP_INIT" -A5 $C/src/cosp.F90
```

The source fingerprints are authoritative (a checkout may be locally modified):

| Fingerprint | COSP `v2.1.4cesm` | COSP `v2.2.1` |
|---|---|---|
| `type cosp_optical_inputs` defined in | `src/cosp.F90` (module `mod_cosp`) | `src/cosp_stats.F90` (module `mod_cosp_stats`) |
| `type swath_inputs` / `compute_orbitmasks` | absent | present |
| `COSP_INIT` last arguments | `..., Nlevels, cloudsat_micro_scheme)` | `..., rttov_Ninstruments, rttov_instrument_namelists, rttov_configs, unitn, debug)` |
| `cosp_column_inputs%sza` | absent (only scalar `zenang` for RTTOV) | present |
| `COSP_ASSIGN_*IN` routines in `src/simulator/cosp_*_interface.F90` | absent (inputs assigned inline in `cosp.F90`) | present |
| MODIS liquid re histogram bins (`nReffLiq`) | 6 | 8 |

| CAM-side fingerprint | `cam_cesm2_1_rel` / `cam_cesm2_2_rel` | `cam_development` |
|---|---|---|
| COSP pin | `Externals_CAM.cfg`: `tag = v2.1.4cesm` | `.gitmodules`: `fxtag = v2.2.1` |
| `src/physics/cam/cospsimulator_intr.F90` | 3748 lines, identical in 2.1 and 2.2 | 4379 lines |
| `cospsimulator_nl` has `COSP_N_SWATHS_*` | no | yes (line 347 onward) |
| `cospstateIN%sza` filled | no (field does not exist) | yes, line 2360: `acos(coszrs)*180/pi` |
| Levels passed to COSP | `1:pver` | `ktop:pver`, `ktop = trop_cloud_top_lev` (line 16), `nlay = pver-ktop+1` (line 621) |
| Microphysics module with pbuf `'MU'` | `micro_mg_cam.F90` line 459 (MG2) | `micro_pumas_cam.F90` line 608 (PUMAS) |

Decision:

* COSP is **v2.2.x** (typically `cam_development`, which still supports `-phys cam6`):
  use **Path B** (section 4). The HARP2 branch drops in.
* COSP is **v2.1.4cesm** (CESM2.1 or CESM2.2 releases, the CMIP6-era CAM6): use
  **Path A** (section 5) with the back-port branch `claude/harp2-cosp-v2.1.4cesm`. Do
  **not** use the v2.2 branch in a release CAM: the release `cospsimulator_intr.F90` is
  written for the v2.1.4 API (types in `mod_cosp`, shorter `COSP_INIT`, no swath type) and
  will not compile against v2.2.1.
* Any other tag (v2.1.5 to v2.1.9, v2.2.0): use the fingerprints above to pick the
  nearer path.

## 3. Facts common to both paths (CAM6 physics)

### 3.1 Effective radius

Both CAM interfaces build the MODIS inputs inside their own `subsample_and_optics`
(release: routine at line 2787, MODIS block at line 3213; development: MODIS block
around lines 3610 to 3665). The liquid size on subcolumns, `MODIS_waterSize` (meters),
comes from `REL` (stratiform, pbuf) and `CV_REFFLIQ` (convective) through
`cosp_simulator_optics`. Fill HARP2 with

```fortran
cospIN%reffLiq = MODIS_waterSize*1.0e6_wp   ! microns
```

`cospIN%fracLiq` is computed by `modis_optics` only inside `if (lmodis_sim)`. For the
first integration, **require `cosp_lmodis_sim=.true.` whenever HARP2 is on**.

In the CESM2.1/2.2 release, CAM's own `optics/cosp_optics.F90` includes stratiform snow:
`cospIN%tau_067` contains the snow optical depth (release lines 3195 to 3206) and
`modis_optics` returns `fracLiq = tau_liq/(tau_liq + tau_ice + tau_snow)` (line 3257). HARP2
therefore treats ice and snow alike, as attenuating but non-polarizing. No change needed.

### 3.2 Effective variance from the microphysics

MG2/PUMAS use a gamma droplet size distribution `n(D) ~ D^mu exp(-lambda D)` and store
`mu` and `lambda` in the physics buffer as `'MU'` and `'LAMBDAC'` (grid-box, stratiform
liquid, `(pcols,pver)`). Where there is prognosed cloud water, CAM computes
`rel = (mu+3)/(2*lambdac)` (release `micro_mg_cam.F90` line 2547; development
`micro_pumas_cam.F90` lines 2946 and 2971). For this distribution

```
ve = 1/(mu + 3)
```

What `MU` actually contains in CESM2.1 `micro_mg_cam.F90` (verified in `405a3f2`; for
`cam_development` verify the same in `micro_pumas_cam.F90` and the PUMAS submodule under
`src/physics/pumas-frozen`, which was not inspected):

1. `size_dist_param_liq` (`micro_mg_utils.F90`) gives the Martin et al. (1994) value with a
   floor: `pgam = 1/(1-0.7 exp(-0.008 Nc))^2 - 1`, `pgam = max(pgam, 2)`, `Nc` the
   in-cloud droplet number in cm^-3. So prognosed ve <= 0.2, and ve = 0.2 whenever
   Nc > about 63 cm^-3. Where in-cloud water is below `qsmall` it returns the sentinel
   -100.
2. The sentinel is then **reset to 0** (`micro_mg_cam.F90` line 2553, `elsewhere` branch).
   So in the pbuf, `MU = 0` means no stratiform cloud water.
3. Where the stratiform cloud fraction `ast < 1e-4`, `MU` is overwritten with the
   convective value `mucon = 5.3` (line 2641; `LAMBDAC` too), i.e. **ve = 1/8.3 = 0.1205**.
   `REL` is **not** changed by this fallback, so in those cells `REL` and `MU` are not a
   consistent pair. `MU >= 2` therefore does not identify prognosed microphysics.

Recommended mapping (stratiform liquid):

```fortran
! mu, ast: pbuf 'MU' and 'AST', sliced like the other COSP inputs
where (mu > 0._r8)              ! prognosed (>= 2) or fallback (5.3)
   ve_ls = 1._r8/(mu + 3._r8)
elsewhere                       ! MU = 0: no stratiform cloud water
   ve_ls = cosp_harp2_veff_default
end where
! Optional: treat the low-cloud-fraction fallback like the no-cloud case
! where (ast < 1.e-4_r8) ve_ls = cosp_harp2_veff_default
```

Other points:

* Convective liquid (`frac_out == 2`) has no `mu`: keep its ve separately configurable
  (suggested namelist `cosp_harp2_veff_conv`).
* If microphysics runs on subcolumns (`use_subcol_microp`), `'MU'` is registered on
  subcolumns; the default CAM6 configuration does not do this.
* Slice `mu` exactly like the other COSP inputs: `(1:ncol,1:pver)` in the release,
  `(1:ncol,ktop:pver)` in development. Both pass levels in CAM's native top-to-bottom order.
* Map to COSP subcolumns with CAM's own `cosp_simulator_optics` (in
  `src/physics/cosp2/optics/cosp_optics.F90`). Argument order: the **first** 2D input goes
  to `frac_out == 2` (convective), the **second** to `frac_out == 1` (stratiform):

```fortran
call cosp_simulator_optics(nPoints, nColumns, nLevels, cospIN%frac_out, &
                           ve_conv, ve_ls, cospIN%veffLiq)   ! ve_conv, ve_ls: (nPoints,nLevels)
```

Look-up table range and saturation. The table now covers ve from 0.01 to 0.40
(20 nodes, denser above 0.1). Checks at the broad end: the neglected fraction of
geometric cross section outside [r_min, 300 microns] is at most 2e-8 over the whole table
(the generator refuses to run if it exceeds 1e-6); halving the Mie size-parameter step
changes -P12 by at most 0.10% of its peak at ve = 0.40; linear interpolation between ve
nodes is accurate to 0.08% of the peak. Inputs outside the table are still clamped to its
edges in the forward model, but this is no longer silent: `harp2_Input_Clamped_Fraction`
and `harp2_Veff_Limit_Fraction` report it. If either is non-negligible in CAM output,
widen the table (`--ve ...` in `harp2_lut_generator.py`; the gamma distribution requires
ve < 0.5) rather than trusting the saturated values. Evidence that this matters: in the
COSP offline UKMO test (input ve = 0.10 everywhere) the old 0.25-limited table returned a
largest ve of exactly 0.250 (saturated); the extended table returns 0.349. The clamping
diagnostic also revealed that 3.5% of subcolumns there have model radii outside 4 to 30
microns. Vertical variation of re near cloud top broadens the retrieved ve, so retrieved ve near or somewhat above the input ve is
expected, but it is not a validation requirement.

### 3.3 Solar zenith angle and daylight

* HARP2 uses `cospgridIN%sza` (degrees) for the viewing geometry and skips points with
  SZA > 75 or `sunlit == 0`. Both CAM interfaces set `sunlit` from `coszrs > 0`
  (and `cosp_runall`). The SZA must be CAM's actual solar geometry.
* `cam_development` already fills `sza`. The release has no such field (see Path A).

### 3.4 LUT file, namelist, and run validity

* Put `harp2_lut_670nm.txt` in the inputdata tree and pass its path through a new
  `cospsimulator_nl` variable (suggested `cosp_harp2_lut_file`, plus the switch
  `cosp_lharp2_sim` and the variance defaults `cosp_harp2_veff_default`,
  `cosp_harp2_veff_conv`). New namelist variables must also be added to
  `bld/namelist_files/namelist_definition.xml` (group `cospsimulator_nl`) and broadcast in
  `cospsimulator_intr_readnl` like the existing `cosp_l*_sim` switches.
* Every MPI task reads the file once in `COSP_INIT`. On failure COSP prints
  `ERROR (HARP2 simulator): cannot open look-up table ...`, disables HARP2 and fills its
  outputs with `R_UNDEF`, and the run would otherwise continue. **After `COSP_INIT`, call
  `harp2_lut_loaded()` (module `mod_cosp_harp2_interface`) and `endrun` if it is false
  while `cosp_lharp2_sim` is true.**

### 3.5 History output: averaging must use a shared mask

How CAM averages fields registered with `flag_xyfill=.true.`: `hbuf_accum_add`
(`control/cam_history_buffers.F90`, starting at line 60) adds a value only when it is not
the fill value and increments a per-point, per-field counter `nacs`; the time mean divides
each field by its own counter (`control/cam_history.F90` around lines 4530 to 4546).
Consequently a weighted numerator must be missing on exactly the same samples as its
denominator. If the numerator is `R_UNDEF` on samples where the retrieval fraction is 0
(observed but no retrieval), it is averaged over fewer samples and the ratio of the two
means is biased high (one fully retrieved 10 micron scene plus one observed clear scene
gives 1000 / 50 = 20 microns instead of 10).

HARP2 recipe (retrieval-fraction-weighted means):

```fortran
! cf  = cospOUT%harp2_Cloud_Fraction_Liquid_Mean   (%; R_UNDEF where not observed)
! re  = cospOUT%harp2_Cloud_Particle_Size_Liquid_Mean (m; R_UNDEF where cf == 0)
where (cf(:ncol) == R_UNDEF)           ! not observed: missing in numerator AND denominator
   re_w(:ncol) = R_UNDEF
   ve_w(:ncol) = R_UNDEF
elsewhere (cf(:ncol) > 0._r8)          ! observed, with retrievals
   re_w(:ncol) = re(:ncol)*cf(:ncol)
   ve_w(:ncol) = ve(:ncol)*cf(:ncol)
elsewhere                              ! observed, no retrieval: zero, not missing
   re_w(:ncol) = 0._r8
   ve_w(:ncol) = 0._r8
end where
call outfld('REFFCLWHARP2', re_w, pcols, lchnk)   ! time mean re = mean(re_w)/mean(CLWHARP2)
```

* The joint histogram and the three diagnostics from COSP are already consistent: they are
  0 (not missing) on observed points without retrievals and `R_UNDEF` only where HARP2 does
  not observe.
* Keep the conditional re and ve themselves missing where there are no retrievals if you
  also output them unweighted.
* **The same bias affects CAM's existing MODIS fields.** `REFFCLWMODIS` is set to
  `R_UNDEF` wherever `CLWMODIS` is 0 (release lines 2720 to 2725), so the time mean of
  `REFFCLWMODIS/CLWMODIS` is biased high. For a HARP2 versus MODIS comparison, apply the
  same shared-mask treatment to the MODIS fields used, or compare instantaneous output.
* Registration: follow the MODIS pattern `addfld(..., horiz_only, 'A', units, ...,
  flag_xyfill=.true., fill_value=R_UNDEF)` (release line 986, development line 990).
* New history coordinates, analogous to `cosp_reffliq` (development line 687):
  `cosp_harp2_re` (centers `harp2_histReffCenters`, edges `harp2_histReffEdges`, m),
  `cosp_harp2_ve` (`harp2_histVeffCenters`, `harp2_histVeffEdges`), and
  `cosp_harp2_flag` (`numHARP2Flags` values 0 to 4). All sizes are in `mod_cosp_config`.
* Suggested field names: `CLWHARP2` (%), `REFFCLWHARP2` (m, weighted), `VEFFCLWHARP2`
  (1, weighted), `CLHARP2REFFVEFF` (%), `HARP2FLAGFRAC` (%), `HARP2VEFFLIM` (%),
  `HARP2CLAMP` (%).

### 3.6 Interpretation caveats

This is an idealized cloudbow pseudo-retrieval: single scattering, black surface,
non-polarizing ice, fixed principal-plane view geometry, no measurement noise, no pixel
aggregation. `CLWHARP2` is the fraction of subcolumns with a successful retrieval, not
the physical liquid cloud cover; use `HARP2FLAGFRAC` to separate changes in cloud
occurrence from changes in retrieval selection. Independent validation against vector
radiative transfer and the real HARP2 geometry and product is still needed before the
output is interpreted as a reproduction of the HARP2 product.

### 3.7 Cost

About 13 to 14 microseconds per sunlit subcolumn (gfortran -O3, one core; re-measured
with the 20-node table, unchanged), only on COSP steps (`cosp_nradsteps`).

## 4. Path B: `cam_development` with COSP v2.2.1

### 4.1 Get the code into CAM's COSP checkout

```bash
cd $CAM/src/physics/cosp2/src
git fetch https://github.com/zzbatmos/COSPv2.0_PACE claude/great-lovelace-wsm2hk
git checkout FETCH_HEAD
```

The sparse checkout (`../.cosp_sparse_checkout`, content `/src/`) keeps only `src/`; all
HARP2 files, including the LUT, are under `src/`, so this is fine. The HARP2 changes
touch `src/cosp.F90`, `src/cosp_config.F90`, `src/cosp_stats.F90` and add the files in
section 1. To see exactly what changed: `git diff 5eb05e5 FETCH_HEAD -- src`.

### 4.2 Swath inputs

No change is needed: `cospIN%cospswathsIN` is still `dimension(6)`, so CAM's
`type(swath_inputs),dimension(6) :: cospswathsIN` (line 278) and the whole-array copy
(line 2390) remain valid. To restrict HARP2 to an overpass, set
`cospIN%harp2_swathIN` (e.g. from new `COSP_N_SWATHS_HARP2`, `COSP_SWATH_LOCALTIMES_HARP2`,
`COSP_SWATH_WIDTHS_HARP2` namelist variables) after that copy.

### 4.3 Build: `src/physics/cosp2/Makefile.in`

No change to `bld/configure` is needed; add to `Makefile.in` (recipes are tab-indented,
like the existing rules):

```make
# add to OBJS
        harp2_simulator.o cosp_harp2_interface.o

# add to the dependency list of cosp.o
        cosp_harp2_interface.o harp2_simulator.o

harp2_simulator.o      : cosp_kinds.o cosp_config.o cosp_constants.o cosp_stats.o
cosp_harp2_interface.o : cosp_kinds.o cosp_constants.o cosp_errorHandling.o \
                         harp2_simulator.o cosp_stats.o

harp2_simulator.o : $(COSP_PATH)/src/src/simulator/HARP2_simulator/harp2_simulator.F90
	$(F90) $(F90FLAGS) -c $<

cosp_harp2_interface.o : $(COSP_PATH)/src/src/simulator/cosp_harp2_interface.F90
	$(F90) $(F90FLAGS) -c $<
```

`cosp_harp2_interface.F90` uses `mod_cosp_error` (`errorMessage`), which CAM provides in
`src/physics/cosp2/cosp_errorHandling.F90`.

### 4.4 Edits in `src/physics/cam/cospsimulator_intr.F90`

1. Namelist: `cosp_lharp2_sim`, `cosp_harp2_lut_file`, `cosp_harp2_veff_default`,
   `cosp_harp2_veff_conv` (+ broadcast + `namelist_definition.xml`). Refuse
   `cosp_lharp2_sim` without `cosp_lmodis_sim` (3.1).
2. `COSP_INIT` call (line 1334): append
   `Lharp2=cosp_lharp2_sim, harp2_lut_file=trim(cosp_harp2_lut_file)`; then check
   `harp2_lut_loaded()` (3.4).
3. `construct_cospIN` (allocations at line 3703): allocate `y%reffLiq`, `y%veffLiq`
   `(npoints,ncolumns,nlevels)` when HARP2 is on (`tau_067` and `fracLiq` are already
   allocated there); deallocate in `destroy_cospIN`.
4. `subsample_and_optics`: after the MODIS optics, fill `reffLiq` (3.1) and `veffLiq`
   (3.2). Get `'MU'` (and `'AST'` if used) with `pbuf_get_index`/`pbuf_get_field` in the run
   routine and pass the `(1:ncol,ktop:pver)` slice down.
5. `construct_cosp_outputs`: allocate the `harp2_*` pointers you output (sizes
   `numHARP2ReffBins`, `numHARP2VeffBins`, `numHARP2Flags` from `mod_cosp_config`);
   `destroy_cosp_outputs`: deallocate them.
6. History: `add_hist_coord` and `addfld` in `cospsimulator_intr_init`, `outfld` in the
   run routine, with the shared-mask weighting of 3.5.
7. `sza` is already filled (line 2360); nothing to do.
8. **Which columns COSP runs on**: COSP is run only on columns where at least one output
   of an active simulator is in an active history file (`hist_fld_col_active`, lines 1992
   to 2069; the sunlit flag at line 2295 requires `run_cosp`), unless `cosp_runall` is
   true. Add a `fname_harp2` list (the HARP2 field names) and a `run_harp2` array, set it
   from `hist_fld_col_active` when HARP2 is on, and include `any(run_harp2(:,i))` in the
   `run_cosp` condition. Otherwise a run whose history holds HARP2 fields but no MODIS
   fields gets `R_UNDEF` everywhere.

## 5. Path A: CESM2.1/2.2 release (COSP v2.1.4cesm): use the back-port branch

The back-port is done: branch `claude/harp2-cosp-v2.1.4cesm` of `zzbatmos/COSPv2.0_PACE`,
based on `v2.1.4cesm` (`34d8eef`). What it changes, and how it was verified, is listed in
its `src/simulator/HARP2_simulator/README.md` ("Differences from the COSP v2.2 version").
Line numbers below are for `cam_cesm2_1_rel` `405a3f2`; `cospsimulator_intr.F90` and
`src/physics/cosp2/Makefile.in` are byte-identical in `cam_cesm2_2_rel` `c47c4d7`. Work in
a separate sandbox (or at least a separate case and build), not in one used by queued
experiments.

### 5.1 Get the code into CAM's COSP checkout

Either check the branch out by hand (quick; a later `checkout_externals` would restore
`v2.1.4cesm`):

```bash
cd $CAM/src/physics/cosp2/src
git fetch https://github.com/zzbatmos/COSPv2.0_PACE claude/harp2-cosp-v2.1.4cesm
git checkout FETCH_HEAD
git diff --stat 34d8eef HEAD -- src      # the HARP2 changes under src/
```

or make it durable in `$CAM/Externals_CAM.cfg`, section `[cosp2]`: set
`repo_url = https://github.com/zzbatmos/COSPv2.0_PACE` and replace `tag = v2.1.4cesm` by
`hash = <commit>` (from `git ls-remote https://github.com/zzbatmos/COSPv2.0_PACE
claude/harp2-cosp-v2.1.4cesm`; `branch = claude/harp2-cosp-v2.1.4cesm` also works but is
not reproducible). Keep `sparse = ../.cosp_sparse_checkout`. Then, from `$CAM`:
`./manage_externals/checkout_externals -e Externals_CAM.cfg` and check with `-S`.
If the repository is private, the machine needs GitHub credentials for either way.

The sparse checkout keeps only `src/`, which holds everything CAM needs, including the
LUT, this guide and the `Makefile.in` patch.

### 5.2 Build: patch `src/physics/cosp2/Makefile.in` (required even with HARP2 off)

The back-port's `cosp.F90` uses the HARP2 modules, so CAM's COSP library must compile the
two HARP2 files whether or not HARP2 is switched on (otherwise the build stops with
`Cannot open module file 'mod_cosp_harp2_interface.mod'`). Apply the tested patch from
the CAM root:

```bash
cd $CAM
patch -p1 < src/physics/cosp2/src/src/simulator/HARP2_simulator/cam_release_Makefile.in.patch
```

It adds `cosp_harp2_interface.o harp2_simulator.o` to `OBJS` and to the dependencies of
`cosp.o`, the two dependency lines, and the two compile rules (paths
`$(COSP_PATH)/src/src/simulator/...`, same style as the MODIS interface rule). CAM's
`configure` regenerates the COSP `Makefile` from this template, so rebuild the case from
clean (`./case.build --clean-all`, then `./case.build`). Checked offline: with the patch,
`libcosp.a` builds from CAM's template, CAM's `cosp_kinds.F90` (`wp = dp`),
`cosp_errorHandling.F90` (it provides the `errorMessage` that HARP2 uses) and optics,
with gfortran debug flags and also with `-fdefault-real-8`; without it the build fails
as described.

### 5.3 Edits in `src/physics/cam/cospsimulator_intr.F90` (release)

1. **Namelist** (`cospsimulator_nl`, line 443; broadcasts from line 470; switch logic
   around line 502): add `cosp_lharp2_sim`, `cosp_harp2_lut_file`,
   `cosp_harp2_veff_default`, `cosp_harp2_veff_conv`, broadcast them, add them to
   `bld/namelist_files/namelist_definition.xml`, and `endrun` if `cosp_lharp2_sim` is
   true while `cosp_lmodis_sim` is false (3.1). Broadcast **before** the
   `call setcosp2values(...)` at line 577, because that routine calls `COSP_INIT`.
2. **`COSP_INIT`** (line 324, in `setcosp2values`): append
   `Lharp2=cosp_lharp2_sim, harp2_lut_file=trim(cosp_harp2_lut_file)` to the positional
   call, then
   `if (cosp_lharp2_sim .and. .not. harp2_lut_loaded()) call endrun('HARP2 LUT not loaded')`
   with `use mod_cosp_harp2_interface, only: harp2_lut_loaded`.
3. **Solar zenith angle** (new field): in `construct_cospstateIN` (line 3311) add
   `y%sza(npoints)` to the allocation; in `destroy_cospstateIN` deallocate it; in the run
   routine after `cospstateIN%sunlit = cam_sunlit(1:ncol)` (line 2077) set
   `cospstateIN%sza = acos(max(-1._r8, min(1._r8, coszrs(1:ncol))))*180._r8/pi`
   (`coszrs` is an argument of `cospsimulator_intr_run`, line 1214; `pi` is already
   imported from `physconst` at line 1197). Night points get SZA > 90 and are masked.
4. **Inputs**: `construct_cospIN` (line 3270) already allocates `tau_067` and `fracLiq`
   unconditionally; add `y%reffLiq` and `y%veffLiq` `(npoints,ncolumns,nlevels)` (when
   HARP2 is on) and deallocate them in `destroy_cospIN` (line 3441).
5. **Optics** (`subsample_and_optics`, line 2787; MODIS block lines 3213 to 3262): after
   `call modis_optics(...)` fill `cospIN%reffLiq = MODIS_waterSize*1.0e6_wp` and
   `cospIN%veffLiq` (3.2). `subsample_and_optics` has no pbuf access: in
   `cospsimulator_intr_init` get `pbuf_get_index('MU')` (and `'AST'` if used) next to the
   other indices (lines 1157 onward), in the run routine get the fields with
   `pbuf_get_field`, compute `ve_ls` and `ve_conv` on `(1:ncol,1:pver)`, and pass them as
   two new arguments (like `dtau_s`). The 0.67 micron block (line 3195) already runs when
   MODIS is on.
6. **Outputs**: `construct_cosp_outputs` (line 3334): allocate the `harp2_*` pointers you
   write (sizes `numHARP2ReffBins`, `numHARP2VeffBins`, `numHARP2Flags` from
   `mod_cosp_config`); deallocate them in `destroy_cosp_outputs`.
7. **History**: register next to the MODIS fields (`addfld`, from line 986;
   `add_default`, line 1049) and the new coordinates (3.5); write after the MODIS
   `outfld` calls (around line 2726) with the shared-mask weighting of 3.5. The HARP2
   outputs are already `R_UNDEF` at night, so they need no entry in the `sunlit_passive`
   block (lines 2172 to 2255).
8. **Which columns COSP runs on** (lines 1695 to 1737; the sunlit flag at line 2043
   requires `run_cosp`): as in 4.4 item 8, add `fname_harp2` (next to `fname_modis`, line
   1384), `run_harp2`, the `hist_fld_col_active` loop and `any(run_harp2(:,i))` in the
   `run_cosp` condition.
9. Levels and ordering need no change: the release passes all levels `1:pver`, top to
   bottom, and `cospstateIN%phalf` (`pver+1` interfaces, top value set to 0, line 2081),
   which is what HARP2 expects. CAM calls COSP once per chunk (`start_idx=1,
   stop_idx=ncol`, line 2115).

## 6. Validation checklist inside CAM

Use a separate case and executable (do not share one with other queued experiments).

1. Build with HARP2 compiled in but switched off; run a short case (e.g. 5 days) and
   compare all existing COSP history fields with an unmodified control build: they must
   be bit-for-bit identical (same compiler and flags).
2. Switch HARP2 on (with MODIS on): `harp2_lut_loaded()` is true, the log has no `HARP2`
   error lines, the fields exist.
3. Closure and masks on instantaneous output (e.g. `nhtfrq=1` for a day):
   `HARP2FLAGFRAC` sums to 100 on observed points; `CLHARP2REFFVEFF` summed over bins
   equals `CLWHARP2`; all HARP2 fields are missing at night and for SZA > 75; retrieved
   re (`REFFCLWHARP2/CLWHARP2` per point) lies within 4 to 30 microns.
4. History averaging: construct or find a point with a cloudy, a clear and a night sample
   in one averaging period and check that the averaged `REFFCLWHARP2/CLWHARP2` equals the
   retrieval-weighted mean of the instantaneous values (3.5).
5. Saturation: `HARP2VEFFLIM` and `HARP2CLAMP` should be near zero; if not, widen the
   table (3.2).
6. Physics sanity (not pass/fail): retrieved ve should sit near the model's ve (about 0.2
   for prognosed MG2 clouds with Nc > 63 cm^-3, 0.12 for fallback cells, the namelist
   value for convective liquid), typically broadened a little by vertical re variation.
   HARP2 re versus MODIS re should follow the model's re profile near cloud top (in the
   COSP offline UKMO test HARP2 re was on average 1 micron smaller because the model re
   peaks mid-cloud); remember the MODIS averaging bias in 3.5.
7. Restart: the simulator has no internal state beyond the LUT, which `COSP_INIT` reads
   again on restart; check that a restarted run matches a continuous one.
8. Timing: check the COSP timers against section 3.7.

## 7. Status of earlier open items

1. Swath compatibility: done (separate `harp2_swathIN`, section 4.2).
2. ve range of the table: extended to 0.40 with diagnostics (section 3.2).
3. Back-port for CESM2.1/2.2 (COSP `v2.1.4cesm`): done, branch
   `claude/harp2-cosp-v2.1.4cesm` (section 5). The CAM-side edits of 5.3 are still to be
   made in the CAM checkout.
