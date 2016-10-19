#include "cppdefs.h"
      MODULE ablm_mod
#ifdef ABLM
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This file contains the code for the atmospheric boundary layer      !
!  model implemented in ROMS. ablm_advdif solves the                   !
!  advection/diffusion of the air temperature and humidity fields.     !
!  ablm_diabatic alters the air temperature and humidity fields        !
!  using the surface heat and moisture fluxes calculated in            !
!  bulk_flux.F. The algorithms used here are adapted from the          !
!  cheapAML package of the MITgcm.                                     !
! 
!  References:                                                         !
!                                                                      !
!    Deremble, B., Wienders, N. and Dewar, W.K., 2013: CheapAML: A     !
!      Simple, Atmospheric Boundary Layer Model for Use in Ocean-Only  !
!      Model Calculations, MWR, 141, 2, 809-821                        !
!                                                                      !
!  Adapted for ROMS by Ryan Holmes, University of New South Wales      !
!      ryan.holmes@unsw.edu.au                                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: ablm_advdif !, ablm_diabatic
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ablm_advdif (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 17)
# endif
      CALL ablm_advdif_tile (ng, tile,                                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     GRID(ng) % on_u,                             &
     &                     GRID(ng) % om_v,                             &
     &                     GRID(ng) % pmon_u,                           &
     &                     GRID(ng) % pnom_v,                           &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
# ifdef WET_DRY
     &                     GRID(ng) % rmask_wet,                        &
     &                     GRID(ng) % umask_wet,                        &
     &                     GRID(ng) % vmask_wet,                        &
# endif
!     &                     FORCES(ng) % Hair,                           &
     &                     FORCES(ng) % Tair,                           &
     &                     FORCES(ng) % Uwind,                          &
     &                     FORCES(ng) % Vwind,                          &
     &                     FORCES(ng) % TairB)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 17)
# endif

      RETURN
      END SUBROUTINE ablm_advdif
!
!***********************************************************************
      SUBROUTINE ablm_advdif_tile (ng, tile,                            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs, pm, pn, on_u, om_v,              &
     &                           pmon_u, pnom_v,                        &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef WET_DRY
     &                           rmask_wet, umask_wet, vmask_wet,       &
# endif
     &                           Tair, Uwind, Vwind, TairB)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: Tair(LBi:,LBj:)
!      real(r8), intent(inout) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
      real(r8), intent(in) :: TairB(LBi:,LBj:)
!      real(r8), intent(in) :: HairB(LBi:,LBj:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: Tair(LBi:UBi,LBj:UBj)
!      real(r8), intent(inout) :: Hair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: TairB(LBi:UBi,LBj:UBj)
!      real(r8), intent(in) :: HairB(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: cff
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: TairTen !Tair tendency current
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: TairTenBef !Tair tendency previous
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wVel ! Vertical velocity on rho-points
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uVel ! Zonal wind on u-points
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vVel ! Meridional wind on v-points
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: fZon ! Zonal flux on u-points
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: fMer ! Meridional flux on v-points
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: fVer ! Vertical flux on rho-points
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ABinc ! Adams-Bashforth Increment
     
!     Local variables to add as parameters:
      integer, parameter :: ntatm = 5 ! Number of atmospheric time steps
      real(r8) :: dtatm               ! Atmospheric time step
      real(r8), parameter :: zWind = 10.0_r8 ! Height of wind velocities
      real(r8), parameter :: ablmK = 4.0E01_r8 ! Diffusion coefficient
      real(r8), parameter :: abEps = 0.01_r8
      real(r8) :: abFac
      
# include "set_bounds.h"
!
!=======================================================================
!  Atmospheric boundary layer model Tair and Qair advection/diffusion
!=======================================================================
!

!
!-----------------------------------------------------------------------
!  Calculate RHS of Tair advection/diffusion equation
!-----------------------------------------------------------------------
!
! Calculate vertical velocity from continuity, assuming w(0) = 0
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          wVel(i,j) =  -(pm(i,j)*(uWind(i+1,j)-uWind(i-1,j))            &
     &                + pn(i,j)*(vWind(i,j+1)-vWind(i,j-1)))*zWind
        END DO
      END DO

! Start atmospheric time loop
      dtatm = dt(ng)/REAL(ntatm)

      DO k=1,ntatm

! Compute zonal advective flux using second-order centered advection
        DO j=JstrV-2,Jendp2
          DO i=IstrU-1,Iendp2
             cff=0.25_r8*on_u(i,j)
             fZon(i,j) = cff*(uWind(i-1,j)+uWind(i,j))*                 & 
     &                       (Tair(i-1,j)+Tair(i,j))
          END DO
        END DO

! Compute meridional advective flux using second-order centered advection
        DO j=JstrV-1,Jendp2
          DO i=IstrU-2,Iendp2
             cff=0.25_r8*om_v(i,j)
             fMer(i,j) = cff*(vWind(i,j-1)+vWind(i,j))*                 & 
     &                       (Tair(i,j-1)+Tair(i,j))
          END DO
        END DO

! Compute vertical advective flux
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            fVer(i,j) = wVel(i,j)*Tair(i,j)  ! (Tair(i,j)+Tair(i,j))*0.5_r8
          END DO
        END DO

! Compute zonal diffusive flux
      DO j=JstrV-2,Jendp2
        DO i=IstrU-1,Iendp2
          fZon(i,j) = fZon(i,j) - ablmK*pmon_u(i,j)                     &
     &                            *(Tair(i,j)-Tair(i-1,j))
        END DO
      END DO

! Compute meridional diffusive flux
      DO j=JstrV-1,Jendp2
        DO i=IstrU-2,Iendp2
          fMer(i,j) = fMer(i,j) - ablmK*pnom_v(i,j)                     &
     &                            *(Tair(i,j)-Tair(i,j-1))
        END DO
      END DO

! Exchange boundaries on zonal and meridional fluxes:
!      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!        CALL exchange_u2d_tile (ng, tile,                               &
!     &                          IminS, ImaxS, JminS, JmaxS,             &
!     &                          fZon)
!        CALL exchange_v2d_tile (ng, tile,                               &
!     &                          IminS, ImaxS, JminS, JmaxS,             &
!     &                          fMer)
!      END IF
!      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
!     &                    IminS, ImaxS, JminS, JmaxS,                   &
!     &                    NghostPoints,                                 &
!     &                    EWperiodic(ng), NSperiodic(ng),               &
!     &                    fZon, fMer)

! Compute divergence of the fluxes
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            TairTen(i,j) = ((fZon(i,j)-fZon(i+1,j)) +                   &
     &                     (fMer(i,j)-fMer(i,j+1)))*pm(i,j)*pn(i,j)+    &
     &                     0.5_r8/zWind*fVer(i,j)
          END DO
        END DO

!
!-----------------------------------------------------------------------
!  Extrapolate to new time tendency using adams-bashforth2
!-----------------------------------------------------------------------
! 
! Determine Adams-Bashforth extrpolation factor:
!      IF (<<<First time, have TairTenBef>>>) THEN
        IF (k.gt.1) THEN ! Just start again for now
          abFac = 0.5_r8 + abEps
        ELSE
          abFac = 0.0_r8
        END IF

! Extrapolate forward tendency in time, with AB weights:
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            ABinc(i,j) = abFac*(TairTen(i,j) - TairTenBef(i,j))
            TairTenBef(i,j) = TairTen(i,j)
            TairTen(i,j) = TairTen(i,j) + ABinc(i,j)
          END DO
        END DO

!
!-----------------------------------------------------------------------
!  Step Tair forward in time
!-----------------------------------------------------------------------
!
! Step forward in time
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            Tair(i,j) = Tair(i,j) + dtatm*TairTen(i,j)
          END DO
        END DO

!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Tair)
        END IF

# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 4,                          &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Tair)
# endif

!
!-----------------------------------------------------------------------
!  Reset boundary values to open boundary conditions
!-----------------------------------------------------------------------
!
!        CALL ABLM_COPY_EDGES()

      END DO ! End atmospheric time loop

      RETURN
      END SUBROUTINE ablm_advdif_tile
#endif
      END MODULE ablm_mod
