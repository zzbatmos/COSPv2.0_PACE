! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Unit tests for the HARP2 polarimetric cloudbow simulator (mod_harp2_sim)
!
! Usage: ./test_harp2_simulator [path to harp2_lut_670nm.txt]
! Exits with a non-zero status if any test fails.
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program test_harp2_simulator
  USE COSP_KINDS,      ONLY: wp
  USE MOD_COSP_CONFIG, ONLY: R_UNDEF,numHARP2ReffBins,numHARP2VeffBins,numHARP2Flags
  use mod_harp2_sim
  use mod_cosp_harp2_interface, ONLY: cosp_harp2_init, harp2_view_geometry, HARP2_NVIEW,    &
                                      harp2_lut_loaded
  implicit none

  integer, parameter :: nLevels = 20, nCases = 15
  real(wp), dimension(nLevels+1)       :: pLev
  real(wp), dimension(nCases, nLevels) :: tau, fLiq, reff, veff
  integer,  dimension(nCases)          :: flag
  real(wp), dimension(nCases)          :: reRet, veRet
  logical,  dimension(nCases)          :: clamp, atLim
  real(wp), dimension(nCases, HARP2_NVIEW) :: Rp
  real(wp), dimension(HARP2_NVIEW)     :: muView, scatAngle
  real(wp) :: mu0
  character(len=512) :: lutFile
  character(len=40), dimension(nCases) :: label
  integer :: k, nFail

  ! Column statistics
  real(wp), dimension(1) :: cfLiq, reMean, veMean, limFrac, clampFrac
  real(wp), dimension(1, numHARP2ReffBins, numHARP2VeffBins) :: hist
  real(wp), dimension(1, numHARP2Flags) :: flagFrac

  nFail = 0
  lutFile = '../../src/simulator/HARP2_simulator/harp2_lut_670nm.txt'
  if (command_argument_count() > 0) call get_command_argument(1, lutFile)
  call cosp_harp2_init(trim(lutFile))
  if (.not. harp2_lut_loaded()) then
     print*, 'Could not load the look-up table ', trim(lutFile)
     stop 1
  endif
  print '(a,f6.3,a,3i5)', ' LUT: wavelength ', LUT_wavelength, ' um, dims ', nLUT_theta, nLUT_re, nLUT_ve
  print '(a,f7.4)', ' Rayleigh optical depth: ', rayleigh_OpticalDepth

  ! 20 layers of 50 hPa from the top of the atmosphere (level 1) to the surface
  pLev(1:nLevels+1) = (/ (5000._wp*k, k = 0, nLevels) /)
  call harp2_view_geometry(30._wp, mu0, muView, scatAngle)
  print '(a,i3)', ' SZA = 30 deg, views in cloudbow window: ',                              &
       count(scatAngle >= cloudbow_ThetaMin .and. scatAngle <= cloudbow_ThetaMax)

  tau = 0._wp; fLiq = 0._wp; reff = 0._wp; veff = 0._wp

  ! 1-4) Homogeneous liquid clouds (layers 15-17, 700-850 hPa), tau = 10
  call liquid_cloud(1, 15, 17, 10._wp, 10._wp,  0.05_wp); label(1) = 'liquid re=10 ve=0.05'
  call liquid_cloud(2, 15, 17, 10._wp, 12.3_wp, 0.07_wp); label(2) = 'liquid re=12.3 ve=0.07'
  call liquid_cloud(3, 15, 17, 10._wp, 20._wp,  0.02_wp); label(3) = 'liquid re=20 ve=0.02'
  call liquid_cloud(4, 15, 17, 10._wp, 6._wp,   0.12_wp); label(4) = 'liquid re=6 ve=0.12'
  ! 5) Ice-only cloud
  tau(5, 5:7) = 5._wp/3._wp;                        label(5) = 'ice only'
  ! 6) Thick ice (tau=3) above liquid
  call liquid_cloud(6, 15, 17, 10._wp, 10._wp, 0.05_wp)
  tau(6, 5:7) = 1._wp;                              label(6) = 'ice tau=3 over liquid'
  ! 7) Thin ice (tau=0.1) above liquid
  call liquid_cloud(7, 15, 17, 10._wp, 10._wp, 0.05_wp)
  tau(7, 5) = 0.1_wp;                               label(7) = 'ice tau=0.1 over liquid'
  ! 8) Adiabatic-like cloud: re increases from 7 (base) to 14 microns (top) in optically
  !    thin layers (tau = 0.25 each) so that several layers contribute to the signal
  do k = 10, 17
     call liquid_cloud(8, k, k, 0.25_wp, 14._wp - 1._wp*(k-10), 0.05_wp)
  enddo
  call liquid_cloud(8, 18, 18, 8._wp, 7._wp, 0.05_wp)
  label(8) = 're 14 (top) to 7 (base)'
  ! 9) Clear sky
  label(9) = 'clear'
  ! 10) Thin liquid cloud
  call liquid_cloud(10, 16, 16, 0.5_wp, 10._wp, 0.05_wp); label(10) = 'thin liquid tau=0.5'
  ! 11) Droplets larger than the LUT
  call liquid_cloud(11, 15, 17, 10._wp, 35._wp, 0.05_wp); label(11) = 'liquid re=35 (beyond LUT)'
  ! 12) Two-layer liquid: thin small-droplet layer over large droplets
  call liquid_cloud(12, 12, 12, 0.3_wp, 6._wp, 0.02_wp)
  call liquid_cloud(12, 15, 17, 10._wp, 20._wp, 0.02_wp); label(12) = 'liquid re=6 (tau .3) over re=20'
  ! 13) Broad distribution inside the table (CAM6/MG2-like ve)
  call liquid_cloud(13, 15, 17, 10._wp, 10._wp, 0.35_wp); label(13) = 'liquid re=10 ve=0.35'
  ! 14) ve beyond the table: clamped in the forward model, retrieval at the table limit
  call liquid_cloud(14, 15, 17, 10._wp, 10._wp, 0.45_wp); label(14) = 'liquid re=10 ve=0.45 (beyond LUT)'
  ! 15) Host forgot to set ve (zero): clamped
  call liquid_cloud(15, 15, 17, 10._wp, 10._wp, 0._wp);   label(15) = 'liquid re=10 ve=0 (unset)'

  call harp2_subcolumn(nCases, nLevels, HARP2_NVIEW, mu0, muView, scatAngle, pLev, tau, fLiq, &
                       reff, veff, flag, reRet, veRet, clamp, atLim, polarizedReflectance=Rp)
  call report('With Rayleigh scattering')

  call check(flag(1) == harp2_flagCloudbow .and. abs(reRet(1)*1e6 - 10._wp) < 0.5_wp .and.  &
             abs(veRet(1) - 0.05_wp) < 0.015_wp, 'homogeneous re=10, ve=0.05')
  call check(flag(2) == harp2_flagCloudbow .and. abs(reRet(2)*1e6 - 12.3_wp) < 0.5_wp .and. &
             abs(veRet(2) - 0.07_wp) < 0.02_wp, 'homogeneous re=12.3, ve=0.07')
  call check(flag(3) == harp2_flagCloudbow .and. abs(reRet(3)*1e6 - 20._wp) < 0.5_wp .and.  &
             abs(veRet(3) - 0.02_wp) < 0.01_wp, 'homogeneous re=20, ve=0.02')
  call check(flag(4) == harp2_flagCloudbow .and. abs(reRet(4)*1e6 - 6._wp) < 0.5_wp .and.   &
             abs(veRet(4) - 0.12_wp) < 0.03_wp, 'homogeneous re=6, ve=0.12')
  call check(flag(5) == harp2_flagNoCloudbow, 'ice-only cloud has no cloudbow')
  call check(flag(6) == harp2_flagNoCloudbow, 'thick ice hides the cloudbow')
  call check(flag(7) == harp2_flagCloudbow .and. abs(reRet(7)*1e6 - 10._wp) < 0.5_wp,         &
             'thin ice: cloudbow retrieval still works')
  call check(flag(8) == harp2_flagCloudbow .and. reRet(8)*1e6 > 12._wp .and.                  &
             reRet(8)*1e6 < 14._wp .and. veRet(8) > 0.06_wp,                                  &
             'cloud-top weighting of re and broadened ve')
  call check(flag(9) == harp2_flagClear, 'clear sky')
  call check(flag(10) == harp2_flagCloudbow .and. abs(reRet(10)*1e6 - 10._wp) < 0.5_wp,       &
             'thin liquid cloud')
  call check(flag(11) == harp2_flagFitFailed, 're beyond the LUT is not retrieved')
  call check(flag(12) == harp2_flagFitFailed .or. (flag(12) == harp2_flagCloudbow .and.        &
             reRet(12)*1e6 < 20._wp), 'two-layer liquid')
  call check(all(Rp(9,:) == R_UNDEF .or. abs(Rp(9,:)) < 0.02_wp), 'clear-sky Rayleigh signal')
  call check(flag(13) == harp2_flagCloudbow .and. abs(reRet(13)*1e6 - 10._wp) < 0.5_wp .and. &
             abs(veRet(13) - 0.35_wp) < 0.03_wp .and. .not. atLim(13) .and. .not. clamp(13), &
             'broad distribution ve=0.35 recovered, not at the table limit')
  call check(flag(14) == harp2_flagCloudbow .and. clamp(14) .and. atLim(14) .and.            &
             veRet(14) >= LUT_ve(nLUT_ve), 've beyond the table: clamped and at the limit')
  call check(abs(veRet(14) - veRet(13)) > 0.02_wp, 've 0.35 and 0.45 are distinguished')
  call check(clamp(15) .and. clamp(11) .and. .not. any(clamp(1:10)),                          &
             'input clamping flagged only for out-of-table inputs')
  call check(.not. any(atLim(1:4)), 'no table-limit flag for in-range retrievals')

  ! Without Rayleigh scattering the fit should be nearly exact on LUT nodes
  include_Rayleigh = .false.
  call harp2_subcolumn(nCases, nLevels, HARP2_NVIEW, mu0, muView, scatAngle, pLev, tau, fLiq, &
                       reff, veff, flag, reRet, veRet, clamp, atLim)
  call report('Without Rayleigh scattering')
  call check(abs(reRet(1)*1e6 - 10._wp) < 0.1_wp .and. abs(veRet(1) - 0.05_wp) < 0.005_wp,    &
             'no Rayleigh: re=10, ve=0.05 recovered')
  call check(abs(reRet(3)*1e6 - 20._wp) < 0.1_wp .and. abs(veRet(3) - 0.02_wp) < 0.005_wp,    &
             'no Rayleigh: re=20, ve=0.02 recovered')
  include_Rayleigh = .true.

  ! Geometry in which the cloudbow is not observed (all scattering angles below 120 deg)
  scatAngle(:) = 100._wp + 20._wp*(/ (real(k,wp)/HARP2_NVIEW, k = 1, HARP2_NVIEW) /)
  call harp2_subcolumn(nCases, nLevels, HARP2_NVIEW, mu0, muView, scatAngle, pLev, tau, fLiq, &
                       reff, veff, flag, reRet, veRet, clamp, atLim)
  call check(flag(1) == harp2_flagGeometry .and. flag(9) == harp2_flagClear .and.             &
             reRet(1) == R_UNDEF, 'cloudbow not sampled')

  ! Retrievals at several solar zenith angles
  do k = 0, 70, 10
     call harp2_view_geometry(real(k,wp), mu0, muView, scatAngle)
     call harp2_subcolumn(nCases, nLevels, HARP2_NVIEW, mu0, muView, scatAngle, pLev, tau,   &
                          fLiq, reff, veff, flag, reRet, veRet, clamp, atLim)
     print '(a,i3,a,i3,a,f7.2,a,f6.3)', ' SZA ', k, ': views in window ',                      &
          count(scatAngle >= cloudbow_ThetaMin .and. scatAngle <= cloudbow_ThetaMax),         &
          ', re(case 1) = ', reRet(1)*1e6, ', ve = ', veRet(1)
     call check(flag(1) == harp2_flagCloudbow .and. abs(reRet(1)*1e6 - 10._wp) < 0.5_wp,     &
                'retrieval across solar zenith angles')
  enddo

  ! Column statistics: 10 subcolumns, 4 retrievals
  call harp2_view_geometry(30._wp, mu0, muView, scatAngle)
  call harp2_subcolumn(nCases, nLevels, HARP2_NVIEW, mu0, muView, scatAngle, pLev, tau, fLiq, &
                       reff, veff, flag, reRet, veRet, clamp, atLim)
  call harp2_column(1, 4, reshape(flag(1:4),(/1,4/)), reshape(reRet(1:4),(/1,4/)),            &
                    reshape(veRet(1:4),(/1,4/)), reshape(clamp(1:4),(/1,4/)),                 &
                    reshape(atLim(1:4),(/1,4/)), cfLiq, reMean, veMean, hist, flagFrac,       &
                    limFrac, clampFrac)
  print '(a,f6.1,a,f7.2,a,f6.3,a,f6.1)', ' Column: cloudbow fraction (%) ', cfLiq(1),          &
       ', mean re (um) ', reMean(1)*1e6, ', mean ve ', veMean(1), ', sum(hist) ', sum(hist)
  call check(abs(cfLiq(1) - 100._wp) < 1e-3_wp .and. abs(sum(hist) - 100._wp) < 1e-3_wp .and.  &
             abs(reMean(1) - sum(reRet(1:4))/4._wp) < 1e-9_wp, 'column statistics')
  call harp2_column(1, 2, reshape(flag(5:6),(/1,2/)), reshape(reRet(5:6),(/1,2/)),            &
                    reshape(veRet(5:6),(/1,2/)), reshape(clamp(5:6),(/1,2/)),                 &
                    reshape(atLim(5:6),(/1,2/)), cfLiq, reMean, veMean, hist, flagFrac,       &
                    limFrac, clampFrac)
  call check(cfLiq(1) == 0._wp .and. reMean(1) == R_UNDEF .and. sum(hist) == 0._wp,          &
             'column statistics without retrievals')
  call check(abs(flagFrac(1,harp2_flagNoCloudbow+1) - 100._wp) < 1e-3_wp,                    &
             'column diagnostics: no-cloudbow fraction')
  ! Diagnostics over all cases: outcome fractions sum to 100%, limit and clamping fractions
  call harp2_column(1, nCases, reshape(flag,(/1,nCases/)), reshape(reRet,(/1,nCases/)),       &
                    reshape(veRet,(/1,nCases/)), reshape(clamp,(/1,nCases/)),                 &
                    reshape(atLim,(/1,nCases/)), cfLiq, reMean, veMean, hist, flagFrac,       &
                    limFrac, clampFrac)
  print '(a,5f7.1)', ' Outcome fractions (%) for flags 0-4: ', flagFrac(1,:)
  print '(a,f6.1,a,f6.1)', ' ve-at-limit (%): ', limFrac(1), '   input clamped (%): ', clampFrac(1)
  call check(abs(sum(flagFrac) - 100._wp) < 1e-3_wp .and.                                   &
             abs(limFrac(1) - 100._wp*count(atLim .and. flag == harp2_flagCloudbow)/nCases) < 1e-3_wp .and. &
             abs(clampFrac(1) - 100._wp*count(clamp)/nCases) < 1e-3_wp, 'column diagnostics')

  print*
  if (nFail > 0) then
     print '(i3,a)', nFail, ' test(s) FAILED'
     stop 1
  endif
  print*, 'All HARP2 simulator tests passed'

contains

  subroutine liquid_cloud(i, kTop, kBase, tauCloud, re, ve)
    integer, intent(in)  :: i, kTop, kBase
    real(wp),intent(in)  :: tauCloud, re, ve
    tau (i, kTop:kBase) = tauCloud/(kBase - kTop + 1)
    fLiq(i, kTop:kBase) = 1._wp
    reff(i, kTop:kBase) = re
    veff(i, kTop:kBase) = ve
  end subroutine liquid_cloud

  subroutine report(title)
    character(len=*),intent(in) :: title
    integer :: i
    print*
    print*, title
    print '(a)', '  case  description                          flag    re (um)    ve'
    do i = 1, nCases
       if (flag(i) == harp2_flagCloudbow) then
          print '(i5,2x,a38,i4,f10.2,f9.3)', i, label(i), flag(i), reRet(i)*1e6, veRet(i)
       else
          print '(i5,2x,a38,i4,a)', i, label(i), flag(i), '         --       --'
       endif
    enddo
  end subroutine report

  subroutine check(ok, name)
    logical,intent(in)          :: ok
    character(len=*),intent(in) :: name
    if (ok) then
       print '(a,a)', '   PASS: ', name
    else
       print '(a,a)', '   FAIL: ', name
       nFail = nFail + 1
    endif
  end subroutine check

end program test_harp2_simulator
