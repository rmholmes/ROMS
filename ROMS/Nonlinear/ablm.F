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
     &                           nrhs,                                  &
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

      integer :: i, j
!      integer, parameter :: IterMax = 3
!      real(r8), parameter :: eps = 1.0E-20_r8
!      real(r8) :: Wspeed, ZQoL, ZToL
!      real(r8), dimension(IminS:ImaxS) :: CC
!      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hlv

# include "set_bounds.h"
!
!=======================================================================
!  Atmospheric boundary layer model Tair and Qair advection/diffusion
!=======================================================================
!
! Calculate Tair due to some random function:
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          Tair(i,j) = TairB(i,j) + 1.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Tair)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Tair)
# endif

      RETURN
      END SUBROUTINE ablm_advdif_tile
#endif
      END MODULE ablm_mod
