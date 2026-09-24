# Integrating the COSP HARP2 simulator into CESM/CAM: version guide

Audience: an agent (or developer) working in a CESM/CAM checkout that wants to add the
HARP2 simulator from this repository. Read section 2 first: the right integration path
depends on which COSP version your CAM uses, and the two paths are not interchangeable.

All facts below were checked against these commits (line numbers refer to them; `grep`
in your own checkout before editing, since line numbers drift):

| Code | Branch / tag | Commit |
|---|---|---|
| HARP2 COSP branch | `zzbatmos/COSPv2.0_PACE`, branch `claude/great-lovelace-wsm2hk` | base `5eb05e5` (= CFMIP COSP `v2.2.1`) + HARP2 commits |
| CAM for CESM2.1 | `ESCOMP/CAM` `cam_cesm2_1_rel` (latest tag `cam_cesm2_1_rel_60`) | `405a3f2` |
| CAM for CESM2.2 | `ESCOMP/CAM` `cam_cesm2_2_rel` (latest tag `cam_cesm2_2_rel_09`) | `c47c4d7` |
| CAM development | `ESCOMP/CAM` `cam_development` | `e9e28e4` (2026-09-18) |
| COSP used by CESM2.1/2.2 | `CFMIP/COSPv2.0` tag `v2.1.4cesm` | `34d8eef` (2019-10-15) |

## 1. What the HARP2 branch provides

* Physics: `src/simulator/HARP2_simulator/harp2_simulator.F90` (module `mod_harp2_sim`).
  Single-scattering polarized reflectance at 670 nm for 60 HARP2 view angles, then a
  Breon and Goloub (1998) cloudbow fit `A*(-P12(re,ve)) + B cos^2 + C` over 135 to 165
  degrees. See `README.md` in the same directory for the method.
* Interface: `src/simulator/cosp_harp2_interface.F90` (module `mod_cosp_harp2_interface`):
  parameters, LUT reader, viewing geometry, `COSP_ASSIGN_harp2IN` (daylight and swath
  masks).
* Look-up table: `src/simulator/HARP2_simulator/harp2_lut_670nm.txt` (1.4 MB text,
  generated with miepython by `harp2_lut_generator.py`).
* New COSP inputs (in `cosp_optical_inputs`, v2.2.1 layout):
  `cospIN%reffLiq(npoints,ncolumns,nlevels)` liquid effective radius in **microns**,
  `cospIN%veffLiq(npoints,ncolumns,nlevels)` liquid effective variance.
  Also used: `cospIN%tau_067`, `cospIN%fracLiq`, `cospgridIN%phalf`, `cospgridIN%sunlit`,
  `cospgridIN%sza` (degrees; defaults to 30 if not allocated).
* New `cosp_outputs` pointers: `harp2_Cloud_Fraction_Liquid_Mean` (%),
  `harp2_Cloud_Particle_Size_Liquid_Mean` (m), `harp2_Effective_Variance_Liquid_Mean` (1),
  `harp2_Reff_vs_Veff_Liquid(npoints,numHARP2ReffBins,numHARP2VeffBins)` (%).
* `COSP_INIT` gained two trailing optional arguments: `Lharp2`, `harp2_lut_file`.
  Existing positional calls are unaffected.
* Swath array `cospIN%cospswathsIN` was enlarged from `dimension(6)` to `dimension(7)`
  (index 7 = HARP2). **This breaks hosts that copy a size-6 array into it (see 4.2).**
* Tested only in the COSP offline driver: unit tests pass; with HARP2 off or on, every
  pre-existing output is bit-identical to COSP `v2.2.1` (also with chunking, model
  levels and swathing). Not yet built or run inside CAM.

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
  **Path A** (section 5), a back-port. Do **not** replace COSP with the HARP2 branch in a
  release CAM: the release `cospsimulator_intr.F90` is written for the v2.1.4 API (types
  in `mod_cosp`, shorter `COSP_INIT`, no swath type) and will not compile against v2.2.1.
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

`cospIN%fracLiq` is computed by `modis_optics` only inside `if (lmodis_sim)`. Either
require `cosp_lmodis_sim=.true.` whenever HARP2 is on, or widen that condition.

### 3.2 Effective variance from the microphysics (recommended over a constant)

MG2/PUMAS use a gamma droplet size distribution `n(D) ~ D^mu exp(-lambda D)` and store
`mu` and `lambda` in the physics buffer as `'MU'` and `'LAMBDAC'` (grid-box, stratiform
liquid, `(pcols,pver)`). CAM computes `rel = (mu+3)/(2*lambdac)` from them
(release `micro_mg_cam.F90` line 2528; development `micro_pumas_cam.F90` lines 2946 and
2971), which confirms that for this distribution

```
ve = 1/(mu + 3)
```

Important details (verified in CESM2.1 `micro_mg_utils.F90`, `size_dist_param_liq`;
for `cam_development` verify the same in the PUMAS submodule under
`src/physics/pumas-frozen`, which was not inspected):

* `mu` follows Martin et al. (1994) with a floor: `pgam = 1/(1-0.7 exp(-0.008 Nc))^2 - 1`,
  `pgam = max(pgam, 2)`, with `Nc` the in-cloud droplet number in cm^-3. So **ve <= 0.2**,
  and ve = 0.2 whenever Nc > about 63 cm^-3 (most clouds).
* Where there is no liquid, `mu` is set to `-100` (sentinel), and above `top_lev` it stays
  0. Only `mu >= 2` is valid. Guard it:
  `where (mu >= 2._r8) ve_ls = 1._r8/(mu+3._r8) elsewhere ve_ls = ve_default`.
* Convective liquid (`frac_out == 2`) has no `mu`: use a namelist constant
  (suggested name `cosp_harp2_veff_conv`).
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

Scientific consequence: with ve near 0.2 in the model, retrieved ve will be near or
above 0.2 (vertical variation of re near cloud top broadens it further). The shipped LUT
stops at **ve = 0.25**, so retrievals may pile up at the edge. Regenerate the LUT with a
wider ve range before production runs (the gamma distribution is valid for ve < 0.5):

```bash
pip install miepython numba
MIEPYTHON_USE_JIT=1 python3 harp2_lut_generator.py -o harp2_lut_670nm.txt \
  --ve 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.10 0.12 0.15 0.20 0.25 0.30 0.35
```

(About 5 minutes. Check whether the branch already contains an extended LUT: see section 7.)

### 3.3 Solar zenith angle and daylight

* HARP2 uses `cospgridIN%sza` (degrees) for the viewing geometry and skips points with
  SZA > 75 or `sunlit == 0`. Both CAM interfaces set `sunlit` from `coszrs > 0`
  (and `cosp_runall`).
* `cam_development` already fills `sza`. The release has no such field (see Path A).

### 3.4 LUT file and namelist

* Put `harp2_lut_670nm.txt` in the inputdata tree and pass its path through a new
  `cospsimulator_nl` variable (suggested `cosp_harp2_lut_file`, plus the switch
  `cosp_lharp2_sim`). New namelist variables must also be added to
  `bld/namelist_files/namelist_definition.xml` (group `cospsimulator_nl`) and broadcast in
  `cospsimulator_intr_readnl` like the existing `cosp_l*_sim` switches.
* Every MPI task reads the file once in `COSP_INIT`. On failure COSP prints
  `ERROR (HARP2 simulator): cannot open look-up table ...`, disables HARP2 and fills its
  outputs with `R_UNDEF`; other simulators are unaffected.

### 3.5 History output conventions in CAM

* Follow the MODIS pattern: `addfld(..., horiz_only, 'A', units, ..., flag_xyfill=.true.,
  fill_value=R_UNDEF)` (release line 986, development line 990).
* CAM writes cloud-fraction-weighted means for MODIS so that time averages are correct,
  e.g. `REFFCLWMODIS` is "MODIS Liquid Cloud Particle Size*CLWMODIS" (release line 1014;
  weighting code at line 2720). Do the same for HARP2: output
  `REFFCLWHARP2 = re*CLWHARP2` and `VEFFCLWHARP2 = ve*CLWHARP2`, set to `R_UNDEF` where
  either is `R_UNDEF`.
* The joint histogram needs two new history coordinates, analogous to `cosp_reffliq`
  (development line 687): e.g. `cosp_harp2_re` (centers `harp2_histReffCenters`, edges
  `harp2_histReffEdges`, in m) and `cosp_harp2_ve` (`harp2_histVeffCenters`,
  `harp2_histVeffEdges`). All are in `mod_cosp_config`.
* Suggested field names: `CLWHARP2` (%), `REFFCLWHARP2` (m), `VEFFCLWHARP2` (1),
  `CLHARP2REFFVEFF` (%).

### 3.6 Cost

About 14 microseconds per sunlit subcolumn (gfortran -O3, one core), only on COSP
steps (`cosp_nradsteps`).

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

### 4.2 Known compile break and its fix

`cospsimulator_intr.F90` declares `type(swath_inputs),dimension(6) :: cospswathsIN`
(line 278) and copies it whole: `cospIN%cospswathsIN = cospswathsIN` (line 2390). With the
HARP2 branch the COSP component has `dimension(7)`, so this is a shape mismatch at
compile time, even with HARP2 off. First check the current branch:

```bash
grep -n "type(swath_inputs)" $C/src/cosp_stats.F90
```

If it says `dimension(7)`, change CAM's declaration to `dimension(7)` (leave element 7
with `N_inst_swaths = 0`, or add `COSP_N_SWATHS_HARP2` etc. to the namelist for a HARP2
overpass mask, index 7). If the branch has been changed back to `dimension(6)` with a
separate HARP2 swath component (planned, see section 7), no CAM change is needed for
compatibility.

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

1. Namelist: `cosp_lharp2_sim`, `cosp_harp2_lut_file`, `cosp_harp2_veff_conv`
   (+ broadcast + `namelist_definition.xml`).
2. `COSP_INIT` call (line 1334): append
   `Lharp2=cosp_lharp2_sim, harp2_lut_file=trim(cosp_harp2_lut_file)`.
3. `construct_cospIN` (allocations at line 3703): allocate `y%reffLiq`, `y%veffLiq`
   `(npoints,ncolumns,nlevels)` when HARP2 is on (`tau_067` and `fracLiq` are already
   allocated there); deallocate in `destroy_cospIN`.
4. `subsample_and_optics`: after the MODIS optics, fill `reffLiq` (3.1) and `veffLiq`
   (3.2). Get `'MU'` with `pbuf_get_index`/`pbuf_get_field` in the run routine and pass
   the `(1:ncol,ktop:pver)` slice down.
5. `construct_cosp_outputs`: allocate the four `harp2_*` pointers
   (sizes from `numHARP2ReffBins`, `numHARP2VeffBins` in `mod_cosp_config`);
   `destroy_cosp_outputs`: deallocate them.
6. History: `add_hist_coord` and `addfld` in `cospsimulator_intr_init`, `outfld` in the
   run routine, with the weighting in 3.5.
7. `sza` is already filled (line 2360); nothing to do.

## 5. Path A: CESM2.1/2.2 release (COSP v2.1.4cesm), back-port

The HARP2 physics is portable; the plumbing must follow the v2.1.4 layout. Work on a
local branch inside `$CAM/src/physics/cosp2/src` (a `v2.1.4cesm` checkout). Use the
HARP2 branch as the template: `git diff 5eb05e5 <HARP2 branch> -- src/cosp.F90`.

1. **Copy unchanged:** `harp2_simulator.F90`, the LUT, `harp2_lut_generator.py`. The
   core only needs `hist2D` (`mod_cosp_stats`), `R_UNDEF`, `pi`, `wp`; all exist in
   v2.1.4cesm.
2. **`cosp_config.F90`:** copy the "HARP2 simulator ReffLIQ/VeffLIQ joint-histogram"
   block. It reuses `nReffLiq`, `reffLIQ_binBounds`, `reffLIQ_binCenters`,
   `reffLIQ_binEdges`, which exist in v2.1.4cesm (6 bins there, so the HARP2 re axis gets
   6 bins).
3. **`cosp_harp2_interface.F90`:** remove `use mod_cosp_stats, only: compute_orbitmasks,
   cosp_optical_inputs, cosp_column_inputs` (in v2.1.4 the types live in `mod_cosp`, which
   uses the interface modules, so importing them here would be circular). Keep the
   parameters, `COSP_HARP2_INIT`, `READ_HARP2_LUT`, `uniform_grid`,
   `HARP2_VIEW_GEOMETRY` and the `harp2_IN` type. Delete `COSP_ASSIGN_harp2IN` and its
   `_CLEAN` routine; do that assignment inline in `cosp.F90`, as v2.1.4 does for MODIS
   (v2.1.4cesm `cosp.F90` around line 680: pointer assignments plus `pack` of sunlit
   indices), without any swath logic.
4. **`cosp.F90` (v2.1.4):** add `reffLiq`, `veffLiq` to `cosp_optical_inputs`; add
   `real(wp),allocatable,dimension(:) :: sza` to `cosp_column_inputs`; add the four output
   pointers to `cosp_outputs`; add `Lharp2_subcolumn/column` switches, the error checks,
   the subcolumn loop, the column statistics and the night fill (copy from the HARP2
   branch `cosp.F90`, adapting names); add optional trailing arguments
   `Lharp2, harp2_lut_file` to `COSP_INIT` (v2.1.4's `COSP_INIT` has no optional
   arguments yet; appending optional ones keeps CAM's positional call valid).
5. **CAM release `cospsimulator_intr.F90`:** as in 4.4, plus allocate and fill
   `cospstateIN%sza = acos(coszrs(1:ncol))*180._r8/pi` in `construct_cospstateIN` / after
   the `sunlit` block (around line 2035 to 2077). Levels are `1:pver`. No swath changes.
6. **Build:** same `Makefile.in` additions as 4.3 (the release file uses the same rule
   style).

## 6. Validation checklist inside CAM

1. Build with HARP2 compiled in but switched off; run a short case (e.g. 5 days) and
   compare all existing COSP history fields with an unmodified control build: they must
   be bit-for-bit identical (same compiler and flags).
2. Switch HARP2 on (with MODIS on): the log has no `HARP2` error lines; fields exist.
3. Sanity: `CLWHARP2` is 0 to 100 and missing at night; `REFFCLWHARP2/CLWHARP2` is within
   4 to 30 microns; `VEFFCLWHARP2/CLWHARP2` is around 0.2 or larger with MG2/PUMAS `mu`
   (with a constant ve it should sit near that constant); the joint histogram summed
   over bins equals `CLWHARP2`.
4. Physics check: compare HARP2 re with `REFFCLWMODIS/CLWMODIS`. In the COSP offline test
   (UKMO data) HARP2 re was on average 1 micron smaller, because the model re peaks in
   mid-cloud and HARP2 samples only the top few tenths of optical depth. The sign of the
   difference in CAM should follow the model's re profile near cloud top.
5. Timing: check the COSP timers against section 3.6.

## 7. Open items on the COSP side (check the branch log before starting)

`git log --oneline FETCH_HEAD | head` (or on GitHub) shows whether these were done after
this document was written:

1. Make the swath change non-breaking (give HARP2 its own swath component and keep
   `cospswathsIN` at `dimension(6)`), removing the need for 4.2.
2. Extend the LUT ve range to about 0.35 (section 3.2).
