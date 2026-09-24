! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2026, University of Maryland Baltimore County
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History
! Sep 2026 - Original version
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_HARP2_INTERFACE
  USE COSP_KINDS,          ONLY: wp
  USE COSP_MATH_CONSTANTS, ONLY: pi
  USE MOD_COSP_ERROR,      ONLY: errorMessage
  use mod_harp2_sim,       ONLY: min_OpticalThickness,cloudbow_ThetaMin,cloudbow_ThetaMax,  &
                                 cloudbow_AmplitudeMin,fit_QualityMin,                     &
                                 rayleigh_OpticalDepth,rayleigh_Depolarization,            &
                                 min_NumAnglesInWindow,include_Rayleigh,                   &
                                 nLUT_theta,nLUT_re,nLUT_ve,LUT_wavelength,LUT_theta,      &
                                 LUT_re,LUT_ve,LUT_ssa,LUT_mP12,scattering_angle
  use mod_cosp_stats,      ONLY: compute_orbitmasks,cosp_optical_inputs,cosp_column_inputs
  implicit none

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Viewing geometry of the HARP2 hyper-angular (670 nm) band
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  integer,parameter :: &
       HARP2_NVIEW = 60              ! Number of along-track view angles
  real(wp) :: &
       harp2_viewZenithMax,        & ! Largest along-track view zenith angle (deg)
       harp2_default_sza,          & ! Solar zenith angle (deg) used when cospgridIN%sza is
                                     ! not provided
       harp2_max_sza                 ! No retrievals for larger solar zenith angles (deg)
  real(wp),dimension(HARP2_NVIEW) :: &
       harp2_viewZenith,           & ! View zenith angles (deg)
       harp2_relAzimuth              ! Relative azimuth between sun and sensor (deg)

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  TYPE harp2_in
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type harp2_IN
     integer,pointer :: &
          Npoints,        & ! Number of horizontal gridpoints
          Ncolumns,       & ! Number of subcolumns
          Nlevels           ! Number of vertical levels
     integer :: &
          Nsunlit           ! Number of sunlit (and swathed) gridpoints
     integer,allocatable,dimension(:) :: &
          sunlit,         & ! Indices of sunlit scenes
          notSunlit         ! Indices of dark (or not swathed) scenes
     real(wp),allocatable,dimension(:) :: &
          sza               ! Solar zenith angle (deg)
     real(wp),pointer ::  &
          pres(:,:),      & ! Gridmean pressure at layer edges (Pa)
          tau(:,:,:),     & ! Subcolumn optical thickness @ 0.67 microns
          liqFrac(:,:,:), & ! Liquid fraction of the optical thickness
          reffLiq(:,:,:), & ! Subcolumn liquid effective radius (microns)
          veffLiq(:,:,:)    ! Subcolumn liquid effective variance
  end type harp2_IN

contains
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_harp2_init
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_HARP2_INIT(lut_file)
    character(len=*),intent(in) :: &
         lut_file                    ! Path to the -P12 look-up table (harp2_lut_generator.py)
    integer  :: i
    real(wp) :: rhoN

    ! Retrieval parameters
    min_OpticalThickness  = 0.3_wp     ! Minimum column optical thickness (as for MODIS)
    cloudbow_ThetaMin     = 135._wp    ! Cloudbow fitting window (scattering angle, deg)
    cloudbow_ThetaMax     = 165._wp    !
    min_NumAnglesInWindow = 10         ! Minimum number of view angles inside the window
    cloudbow_AmplitudeMin = 0.2_wp     ! Minimum cloudbow amplitude, relative to that of an
                                       ! optically thick liquid cloud without anything above
    fit_QualityMin        = 0.8_wp     ! Minimum fraction of the non-smooth signal that the
                                       ! cloudbow term must explain

    ! Rayleigh scattering (Hansen and Travis 1974, eq. 2.29, and King factor rho_n = 0.0279)
    include_Rayleigh      = .true.
    rhoN                  = 0.0279_wp
    rayleigh_Depolarization = (1._wp - rhoN)/(1._wp + rhoN/2._wp)

    ! Viewing geometry. Along-track views in the solar principal plane, which is close to
    ! the geometry near the sub-satellite track of an early-afternoon sun-synchronous orbit.
    ! Positive view zenith angles look toward the sun (backscattering side).
    harp2_viewZenithMax = 57._wp
    harp2_default_sza   = 30._wp
    harp2_max_sza       = 75._wp
    do i = 1, HARP2_NVIEW
       harp2_viewZenith(i) = -harp2_viewZenithMax + 2._wp*harp2_viewZenithMax*(i-1)/(HARP2_NVIEW-1)
    enddo
    harp2_relAzimuth(:) = merge(0._wp, 180._wp, harp2_viewZenith(:) >= 0._wp)
    harp2_viewZenith(:) = abs(harp2_viewZenith(:))

    ! Polarized phase function look-up table
    call read_harp2_lut(lut_file)
    if (nLUT_re > 0) rayleigh_OpticalDepth = 0.008569_wp/LUT_wavelength**4 *              &
         (1._wp + 0.0113_wp/LUT_wavelength**2 + 0.00013_wp/LUT_wavelength**4)

  END SUBROUTINE COSP_HARP2_INIT

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE read_harp2_lut
  ! Reads the text file written by harp2_lut_generator.py. On failure an error message is
  ! issued and the LUT dimensions are left at zero, which disables the HARP2 simulator in
  ! cosp_errorCheck.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE READ_HARP2_LUT(lut_file)
    character(len=*),intent(in) :: lut_file
    integer :: iu, ios, nT, nR, nV
    character(len=512) :: line
    real(wp),allocatable,dimension(:,:) :: qext

    if (allocated(LUT_theta)) deallocate(LUT_theta, LUT_re, LUT_ve, LUT_ssa, LUT_mP12)
    nLUT_theta = 0; nLUT_re = 0; nLUT_ve = 0

    open(newunit=iu, file=trim(lut_file), status='old', action='read', form='formatted', &
         iostat=ios)
    if (ios /= 0) then
       call errorMessage('ERROR (HARP2 simulator): cannot open look-up table '//trim(lut_file))
       return
    endif

    ! Skip the comment lines
    do
       read(iu,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0) cycle
       if (line(1:1) /= '#') exit
    enddo
    if (ios == 0) read(line,*,iostat=ios) nT, nR, nV
    if (ios == 0) then
       if (nT < 2 .or. nR < 3 .or. nV < 2) ios = 1
    endif
    if (ios == 0) then
       allocate(LUT_theta(nT), LUT_re(nR), LUT_ve(nV), LUT_ssa(nR,nV), LUT_mP12(nT,nR,nV),  &
                qext(nR,nV))
       read(iu,*,iostat=ios) LUT_wavelength
       if (ios == 0) read(iu,*,iostat=ios) LUT_theta
       if (ios == 0) read(iu,*,iostat=ios) LUT_re
       if (ios == 0) read(iu,*,iostat=ios) LUT_ve
       if (ios == 0) read(iu,*,iostat=ios) LUT_ssa
       if (ios == 0) read(iu,*,iostat=ios) qext       ! Extinction efficiency (not used yet)
       if (ios == 0) read(iu,*,iostat=ios) LUT_mP12
       deallocate(qext)
    endif
    close(iu)

    if (ios /= 0) then
       call errorMessage('ERROR (HARP2 simulator): cannot read look-up table '//trim(lut_file))
       if (allocated(LUT_theta)) deallocate(LUT_theta, LUT_re, LUT_ve, LUT_ssa, LUT_mP12)
       return
    endif

    ! The interpolation assumes uniform theta and re grids and increasing ve
    if (.not. uniform_grid(LUT_theta) .or. .not. uniform_grid(LUT_re) .or.               &
        any(LUT_ve(2:nV) <= LUT_ve(1:nV-1))) then
       call errorMessage('ERROR (HARP2 simulator): theta and re must be uniformly spaced and ve increasing in '// &
                         trim(lut_file))
       deallocate(LUT_theta, LUT_re, LUT_ve, LUT_ssa, LUT_mP12)
       return
    endif
    nLUT_theta = nT; nLUT_re = nR; nLUT_ve = nV

  END SUBROUTINE READ_HARP2_LUT

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! FUNCTION uniform_grid: increasing and uniformly spaced (to 0.1% of the spacing)
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  pure function uniform_grid(x)
    real(wp),dimension(:),intent(in) :: x
    logical :: uniform_grid
    real(wp) :: dx
    integer  :: n

    n  = size(x)
    dx = (x(n) - x(1))/(n - 1)
    uniform_grid = dx > 0._wp .and. all(abs(x(2:n) - x(1:n-1) - dx) <= 1.e-3_wp*dx)
  end function uniform_grid

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE harp2_view_geometry
  ! Cosines of the view zenith angles and scattering angles of the HARP2 views for a given
  ! solar zenith angle (deg).
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE HARP2_VIEW_GEOMETRY(sza, mu0, muView, scatAngle)
    real(wp),intent(in)  :: sza
    real(wp),intent(out) :: mu0
    real(wp),intent(out),dimension(HARP2_NVIEW) :: muView, scatAngle

    mu0       = cos(sza*pi/180._wp)
    muView    = cos(harp2_viewZenith*pi/180._wp)
    scatAngle = scattering_angle(mu0, muView, harp2_relAzimuth)
  END SUBROUTINE HARP2_VIEW_GEOMETRY

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE COSP_ASSIGN_harp2IN
  ! As for MODIS, harp2IN covers all gridpoints; the retrievals are only done for the
  ! sunlit (and swathed) gridpoints listed in harp2IN%sunlit.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_harp2IN(cospIN,cospgridIN,Npoints,harp2IN)
     type(cosp_optical_inputs),intent(in),target :: cospIN     ! Optical inputs to COSP simulator
     type(cosp_column_inputs), intent(in),target :: cospgridIN ! Host model inputs to COSP
     integer,intent(in),target :: &
         Npoints
     type(harp2_IN),intent(inout) :: &
         harp2IN
     ! Local variables
     logical,dimension(Npoints) :: &
         HARP2_MASK                  ! Gridpoints observed in daylight
     logical,dimension(:),allocatable :: &
         HARP2_SWATH_MASK
     integer :: &
         N_HARP2_SWATHED,  &
         i

     harp2IN%Ncolumns => cospIN%Ncolumns
     harp2IN%Nlevels  => cospIN%Nlevels
     harp2IN%Npoints  => Npoints
     harp2IN%tau      => cospIN%tau_067
     harp2IN%liqFrac  => cospIN%fracLiq
     harp2IN%reffLiq  => cospIN%reffLiq
     harp2IN%veffLiq  => cospIN%veffLiq
     harp2IN%pres     => cospgridIN%phalf

     allocate(harp2IN%sza(Npoints))
     if (allocated(cospgridIN%sza)) then
        harp2IN%sza(1:Npoints) = cospgridIN%sza(1:Npoints)
     else
        harp2IN%sza(1:Npoints) = harp2_default_sza
     endif

     HARP2_MASK(1:Npoints) = (cospgridIN%sunlit(1:Npoints) > 0) .and.                     &
                             (harp2IN%sza(1:Npoints) <= harp2_max_sza)
     if (cospIN % cospswathsIN(7) % N_inst_swaths .gt. 0) then
         allocate(HARP2_SWATH_MASK(Npoints))
         ! Do swathing to figure out which cells to simulate on
         call compute_orbitmasks(Npoints,                                                &
                                 cospIN % cospswathsIN(7) % N_inst_swaths,               &
                                 cospIN % cospswathsIN(7) % inst_localtimes,             &
                                 cospIN % cospswathsIN(7) % inst_localtime_widths,       &
                                 cospgridIN%lat, cospgridIN%lon,                         &
                                 cospgridIN%rttov_date(:,2), cospgridIN%rttov_date(:,3), & ! Time fields: month, dayofmonth
                                 cospgridIN%rttov_time(:,1), cospgridIN%rttov_time(:,2), & ! Time fields: hour, minute
                                 HARP2_SWATH_MASK,N_HARP2_SWATHED) ! Output: logical mask array
         HARP2_MASK(1:Npoints) = HARP2_MASK(1:Npoints) .and. HARP2_SWATH_MASK(1:Npoints)
         deallocate(HARP2_SWATH_MASK)
     endif

     harp2IN%Nsunlit = count(HARP2_MASK)
     allocate(harp2IN%sunlit(harp2IN%Nsunlit), harp2IN%notSunlit(Npoints - harp2IN%Nsunlit))
     harp2IN%sunlit    = pack((/ (i, i = 1, Npoints) /), mask = HARP2_MASK)
     harp2IN%notSunlit = pack((/ (i, i = 1, Npoints) /), mask = .not. HARP2_MASK)

  END SUBROUTINE COSP_ASSIGN_harp2IN

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE COSP_ASSIGN_harp2IN_CLEAN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_harp2IN_CLEAN(harp2IN)
     type(harp2_IN),intent(inout) :: harp2IN

     if (allocated(harp2IN%sza))       deallocate(harp2IN%sza)
     if (allocated(harp2IN%sunlit))    deallocate(harp2IN%sunlit)
     if (allocated(harp2IN%notSunlit)) deallocate(harp2IN%notSunlit)
     nullify(harp2IN%Npoints, harp2IN%Ncolumns, harp2IN%Nlevels, harp2IN%tau,             &
             harp2IN%liqFrac, harp2IN%reffLiq, harp2IN%veffLiq, harp2IN%pres)
  END SUBROUTINE COSP_ASSIGN_harp2IN_CLEAN

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE MOD_COSP_HARP2_INTERFACE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_HARP2_INTERFACE
