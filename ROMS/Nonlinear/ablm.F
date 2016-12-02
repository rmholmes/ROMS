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
!  bulk_flux.F. ablm_nudgcof sets the time-scales for nudging of Tair  !
!  and Qair to data. The algorithms used here are adapted from the     !
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
      PUBLIC  :: ablm_advdif, ablm_diabatic, ablm_nudgcof
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
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     GRID(ng) % on_u,                             &
     &                     GRID(ng) % om_v,                             &
     &                     GRID(ng) % pmon_u,                           &
     &                     GRID(ng) % pnom_v,                           &
     &                     FORCES(ng) % Tair,                           &
     &                     FORCES(ng) % TairTenPrev,                    &
# ifdef ABLM_WATER
     &                     FORCES(ng) % Hair,                           &
     &                     FORCES(ng) % HairTenPrev,                    &
     &                     FORCES(ng) % Pair,                           &
# endif
     &                     FORCES(ng) % Uwind,                          &
     &                     FORCES(ng) % Vwind)
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
     &                           pm, pn, on_u, om_v,                    &
     &                           pmon_u, pnom_v,                        &
     &                           Tair, TairTenPrev,                     &
# ifdef ABLM_WATER
     &                           Hair, HairTenPrev, Pair,               &
# endif
     &                           Uwind, Vwind)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
      USE bc_2d_mod, ONLY: bc_r2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(inout) :: Tair(LBi:,LBj:)
      real(r8), intent(inout) :: TairTenPrev(LBi:,LBj:)
#  ifdef ABLM_WATER
      real(r8), intent(inout) :: Hair(LBi:,LBj:)
      real(r8), intent(inout) :: HairTenPrev(LBi:,LBj:)
      real(r8), intent(inout) :: Pair(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: Uwind(LBi:,LBj:)
      real(r8), intent(inout) :: Vwind(LBi:,LBj:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Tair(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: TairTenPrev(LBi:UBi,LBj:UBj)
#  ifdef ABLM_WATER
      real(r8), intent(inout) :: Hair(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: HairTenPrev(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Pair(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: cff
#ifdef ABLM_WATER
      logical :: Qswitch
#endif
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: TairTen !Tair tendency current
#ifdef ABLM_WATER
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: HairTen !Hair tendency current
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Qair    ! Specific humidity (kg/kg)
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wVel ! Vertical velocity on rho-points
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uVel ! Zonal wind on u-points
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vVel ! Meridional wind on v-points
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: fZon ! Zonal flux on u-points
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: fMer ! Meridional flux on v-points
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: fVer ! Vertical flux on rho-points
#ifdef ABLM_WATER
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: fZonQ ! Zonal flux on u-points, water vapor
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: fMerQ ! Meridional flux on v-points, water vapor
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: fVerQ ! Vertical flux on rho-points, water vapor
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ABinc ! Adams-Bashforth Increment
#ifdef ABLM_WATER
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ABincQ ! Adams-Bashforth Increment Qair
#endif     
      

!     Local variables to add as parameters:
      real(r8) :: dtatm               ! Atmospheric time step
      real(r8) :: abFac
      
# include "set_bounds.h"
!
!=======================================================================
!  Atmospheric boundary layer model Tair and Qair advection/diffusion
!=======================================================================
!
!      
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions
!-----------------------------------------------------------------------
!
      CALL bc_r2d_tile(ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Tair)
#ifdef ABLM_WATER
      CALL bc_r2d_tile(ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Hair)
#endif
      CALL bc_r2d_tile(ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Uwind)
      CALL bc_r2d_tile(ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Vwind)
!
!-----------------------------------------------------------------------
!  Calculate RHS of Tair advection/diffusion equation
!-----------------------------------------------------------------------
!
! Calculate vertical velocity from continuity, assuming w(0) = 0
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          wVel(i,j) =  -(pm(i,j)*(Uwind(i+1,j)-Uwind(i-1,j))            &
     &                + pn(i,j)*(Vwind(i,j+1)-Vwind(i,j-1)))*zWind
        END DO
      END DO

! Start atmospheric time loop
      dtatm = dt(ng)/REAL(ntatm)

      DO k=1,ntatm

#ifdef ABLM_WATER
!-----------------------------------------------------------------------
!  Convert humidity input to specific humidity (kg/kg)
!-----------------------------------------------------------------------
!
      DO j=JstrV-2,Jendp2
        DO i=IstrU-2,Iendp2
          cff=(1.0007_r8+3.46E-6_r8*Pair(i,j))*6.1121_r8*               &
     &        EXP(17.502_r8*Tair(i,j)/(240.97_r8+Tair(i,j)))
          IF (Hair(i,j).lt.2.0_r8) THEN                !RH fraction
            Qswitch = .TRUE.
            cff=cff*Hair(i,j)                          !Vapor pres (mb)
            Qair(i,j)=0.62197_r8*(cff/(Pair(i,j)-0.378_r8*cff)) !Spec hum (kg/kg)
          ELSE          !RH input was actually specific humidity in g/kg
            Qswitch = .FALSE.
            Qair(i,j)=Hair(i,j)/1000.0_r8                       !Spec Hum (kg/kg)
          END IF
        END DO
      END DO
#endif

! Compute zonal advective flux using second-order centered advection
        DO j=JstrV-2,Jendp2
          DO i=IstrU-1,Iendp2
             cff=0.25_r8*on_u(i,j)
             fZon(i,j) = cff*(Uwind(i-1,j)+Uwind(i,j))*                 & 
     &                       (Tair(i-1,j)+Tair(i,j))
#ifdef ABLM_WATER
             fZonQ(i,j) = cff*(Uwind(i-1,j)+Uwind(i,j))*                & 
     &                       (Qair(i-1,j)+Qair(i,j))
#endif
          END DO
        END DO

! Compute meridional advective flux using second-order centered advection
        DO j=JstrV-1,Jendp2
          DO i=IstrU-2,Iendp2
             cff=0.25_r8*om_v(i,j)
             fMer(i,j) = cff*(Vwind(i,j-1)+Vwind(i,j))*                 & 
     &                       (Tair(i,j-1)+Tair(i,j))
#ifdef ABLM_WATER
             fMerQ(i,j) = cff*(Vwind(i,j-1)+Vwind(i,j))*                & 
     &                       (Qair(i,j-1)+Qair(i,j))
#endif
          END DO
        END DO

! Compute vertical advective flux
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            fVer(i,j) = wVel(i,j)*Tair(i,j)  ! (Tair(i,j)+Tair(i,j))*0.5_r8
#ifdef ABLM_WATER
            fVerQ(i,j) = wVel(i,j)*Qair(i,j)
#endif
          END DO
        END DO

! Compute zonal diffusive flux
      DO j=JstrV-2,Jendp2
        DO i=IstrU-1,Iendp2
          fZon(i,j) = fZon(i,j) - ablmK*pmon_u(i,j)                     &
     &                            *(Tair(i,j)-Tair(i-1,j))
#ifdef ABLM_WATER
          fZonQ(i,j) = fZonQ(i,j) - ablmK*pmon_u(i,j)                   &
     &                            *(Qair(i,j)-Qair(i-1,j))
#endif
        END DO
      END DO

! Compute meridional diffusive flux
      DO j=JstrV-1,Jendp2
        DO i=IstrU-2,Iendp2
          fMer(i,j) = fMer(i,j) - ablmK*pnom_v(i,j)                     &
     &                            *(Tair(i,j)-Tair(i,j-1))
#ifdef ABLM_WATER
          fMerQ(i,j) = fMerQ(i,j) - ablmK*pnom_v(i,j)                   &
     &                            *(Qair(i,j)-Qair(i,j-1))
#endif
        END DO
      END DO

! Compute divergence of the fluxes
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            TairTen(i,j) = ((fZon(i,j)-fZon(i+1,j)) +                   &
     &                     (fMer(i,j)-fMer(i,j+1)))*pm(i,j)*pn(i,j)-    &
     &                     0.5_r8/zWind*fVer(i,j)
#ifdef ABLM_WATER
            HairTen(i,j) = ((fZonQ(i,j)-fZonQ(i+1,j)) +                 &
     &                     (fMerQ(i,j)-fMerQ(i,j+1)))*pm(i,j)*pn(i,j)-  &
     &                     0.5_r8/zWind*fVerQ(i,j)
#endif
          END DO
        END DO

!
!-----------------------------------------------------------------------
!  Extrapolate to new time tendency using adams-bashforth2
!-----------------------------------------------------------------------
! 
! Determine Adams-Bashforth extrpolation factor:
        IF (ABfirst) THEN ! First step, don't have TairTenPrev
          abFac = 0.0_r8
          ABfirst = .FALSE.
        ELSE ! Otherwise, do have TairTenPrev
          abFac = 0.5_r8 + ABeps
        END IF

! Extrapolate forward tendency in time, with AB weights:
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            ABinc(i,j) = abFac*(TairTen(i,j) - TairTenPrev(i,j))
            TairTenPrev(i,j) = TairTen(i,j)
            TairTen(i,j) = TairTen(i,j) + ABinc(i,j)
#ifdef ABLM_WATER
            ABincQ(i,j) = abFac*(HairTen(i,j) - HairTenPrev(i,j))
            HairTenPrev(i,j) = HairTen(i,j)
            HairTen(i,j) = HairTen(i,j) + ABincQ(i,j)
#endif
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
#ifdef ABLM_WATER
            Qair(i,j) = Qair(i,j) + dtatm*HairTen(i,j)
#endif
          END DO
        END DO

#ifdef ABLM_WATER
!-----------------------------------------------------------------------
!  Convert specific humidity back to previous units
!-----------------------------------------------------------------------
!
      DO j=JstrV-2,Jendp2
        DO i=IstrU-2,Iendp2
          cff=(1.0007_r8+3.46E-6_r8*Pair(i,j))*6.1121_r8*               &
     &        EXP(17.502_r8*Tair(i,j)/(240.97_r8+Tair(i,j)))
          IF (Qswitch) THEN                !RH fraction
             Hair(i,j) = Qair(i,j)/0.62197_r8*Pair(i,j)/cff/            &
     &                 (1.0_r8 + 0.378_r8/0.62197_r8*Qair(i,j))
          ELSE          !RH input was actually specific humidity in g/kg
            Hair(i,j) = 1000.0_r8*Qair(i,j)
          END IF
        END DO
      END DO
#endif

!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Tair)
#ifdef ABLM_WATER
          CALL exchange_r2d_tile (ng, tile,                             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Hair)
#endif
        END IF

# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 4,                          &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
#ifdef ABLM_WATER
     &                    Hair,                                         &
#endif
     &                    Tair)
# endif

      END DO ! End atmospheric time loop


      RETURN
      END SUBROUTINE ablm_advdif_tile

!
!***********************************************************************
      SUBROUTINE ablm_diabatic (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_stepping
      USE mod_ocean
      USE mod_clima
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
      CALL ablm_diabatic_tile (ng, tile,                                &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     FORCES(ng) % Tair,                           &
     &                     FORCES(ng) % TairB,                          &
#ifdef ABLM_WATER
     &                     FORCES(ng) % Hair,                           &
     &                     FORCES(ng) % HairB,                          &
     &                     FORCES(ng) % Pair,                           &
     &                     FORCES(ng) % Q100,                           &
#endif
     &                     FORCES(ng) % lhflx,                          &
     &                     FORCES(ng) % shflx,                          &
     &                     CLIMA(ng) % ABLMnudgcof,                     &
     &                     OCEAN(ng) % t)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 17)
# endif

      RETURN
      END SUBROUTINE ablm_diabatic
!
!***********************************************************************
      SUBROUTINE ablm_diabatic_tile (ng, tile,                          &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
# ifdef MASKING
     &                           rmask,                                 &
# endif
     &                           Tair, TairB,                           &
# ifdef ABLM_WATER
     &                           Hair, HairB, Pair, Q100,               &
# endif
     &                           lhflx, shflx,                          &
     &                           ABLMnudgcof, t)
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
#  endif
      real(r8), intent(inout) :: Tair(LBi:,LBj:)
      real(r8), intent(in) :: TairB(LBi:,LBj:)
# ifdef ABLM_WATER
      real(r8), intent(inout) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: HairB(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
      real(r8), intent(in) :: Q100(LBi:,LBj:)
# endif
      real(r8), intent(in) :: lhflx(LBi:,LBj:)
      real(r8), intent(in) :: shflx(LBi:,LBj:)
      real(r8), intent(in) :: ABLMnudgcof(LBi:,LBj:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: Tair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: TairB(LBi:UBi,LBj:UBj)
#  ifdef ABLM_WATER
      real(r8), intent(inout) :: Hair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: HairB(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Q100(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ABLMnudgcof(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: cff, HscaleA, TairK, TseaC, TseaK

      real(r8) :: Ttop ! temperature for outgoing longwave calculation
      real(r8) :: TairTen ! temperature for outgoing longwave calculation
      real(r8) :: LWup ! atmospheric upwelled lw radiation
      real(r8) :: OCup ! oceanic upwelled lw radiation
      real(r8) :: LWdn ! atmospheric downwelled lw radiation
      real(r8) :: NETup ! atmospheric upwelled lw radiation
      real(r8) :: NETdn ! atmospheric downwelled lw radiation
#ifdef ABLM_WATER
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Qair    ! Specific humidity (kg/kg)
      real(r8) :: HairTen, Qsea, HscaleW, evap, entrain, lhvap
      logical :: Qswitch
#endif      
!
           
# include "set_bounds.h"
#ifdef ABLM_WATER
!-----------------------------------------------------------------------
!  Convert humidity input to specific humidity (kg/kg)
!-----------------------------------------------------------------------
!
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          cff=(1.0007_r8+3.46E-6_r8*Pair(i,j))*6.1121_r8*               &
     &        EXP(17.502_r8*Tair(i,j)/(240.97_r8+Tair(i,j)))
          IF (Hair(i,j).lt.2.0_r8) THEN                !RH fraction
            Qswitch = .TRUE.
            cff=cff*Hair(i,j)                          !Vapor pres (mb)
            Qair(i,j)=0.62197_r8*(cff/(Pair(i,j)-0.378_r8*cff)) !Spec hum (kg/kg)
          ELSE          !RH input was actually specific humidity in g/kg
            Qswitch = .FALSE.
            Qair(i,j)=Hair(i,j)/1000.0_r8                       !Spec Hum (kg/kg)
          END IF
        END DO
      END DO
#endif
!
!=======================================================================
!  Atmospheric boundary layer model diabatic steps
!=======================================================================
!
!
      HscaleA=ABLM_rhoa*blk_Cpa
# ifdef ABLM_WATER
      HscaleW=1.0_r8/(rho0*Cp)
# endif
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
!
!-----------------------------------------------------------------------
!  Air temperature diabatic changes
!-----------------------------------------------------------------------
! 
          TseaC = t(i,j,N(ng),nrhs,itemp)
          TseaK = TseaC+273.16_r8
          TairK = Tair(i,j)+273.16_r8
          Ttop = TairK - ABLM_lapse * ( hABLM - zWind )

          LWup = 0.5_r8*StefBo*(0.5_r8*TairK+0.5_r8*Ttop)**4
          LWdn = 0.5_r8*StefBo*TairK**4
          OCup = emmiss*StefBo*TseaK**4

          NETup = LWup/HscaleA
          NETdn = (OCup-LWdn-shflx(i,j))/HscaleA

          TairTen = (NETdn - NETup) / hABLM
#ifdef MASKING
          TairTen = TairTen*rmask(i,j)
#endif

          Tair(i,j) = Tair(i,j) + dt(ng)*TairTen

#ifdef ABLM_WATER
!
!-----------------------------------------------------------------------
!  Water vapor diabatic changes
!-----------------------------------------------------------------------
! 
          cff=(1.0007_r8+3.46E-6_r8*Pair(i,j))*6.1121_r8*               &
     &        EXP(17.502_r8*TseaC/(240.97_r8+TseaC))
          cff=cff*0.98_r8
          Qsea=0.62197_r8*(cff/(Pair(i,j)-0.378_r8*cff))

          lhvap = (2.501_r8-0.00237_r8*TseaC)*1.0E+6_r8
          evap = -lhflx(i,j) / HscaleW / lhvap
          entrain = 0.25*(evap / (Qsea - Qair(i,j)))*Q100(i,j)

          HairTen = (evap - entrain) / hABLM / ABLM_rhoa
#ifdef MASKING
          HairTen = HairTen*rmask(i,j)
#endif

          Qair(i,j) = Qair(i,j) + dt(ng)*HairTen
#endif
        END DO
      END DO

#ifdef ABLM_WATER
!-----------------------------------------------------------------------
!  Convert specific humidity back to previous units
!-----------------------------------------------------------------------
!
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          cff=(1.0007_r8+3.46E-6_r8*Pair(i,j))*6.1121_r8*               &
     &        EXP(17.502_r8*Tair(i,j)/(240.97_r8+Tair(i,j)))
          IF (Qswitch) THEN                !RH fraction
            Hair(i,j) = Qair(i,j)/0.62197_r8*Pair(i,j)/cff/             &
     &                 (1.0_r8 + 0.378_r8/0.62197_r8*Qair(i,j))
          ELSE          !RH input was actually specific humidity in g/kg
            Hair(i,j) = 1000.0_r8*Qair(i,j)
          END IF
        END DO
      END DO
#endif

!
!-----------------------------------------------------------------------
!  Nudging (implicit) over land and near boundaries
!-----------------------------------------------------------------------
! 
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          cff = ABLMnudgcof(i,j)*dt(ng)
          Tair(i,j)=Tair(i,j) + (TairB(i,j)-Tair(i,j))*                 &   
     &                           cff/(1.0_r8+cff)
#ifdef ABLM_WATER
          Hair(i,j)=Hair(i,j) + (HairB(i,j)-Hair(i,j))*                 &   
     &                           cff/(1.0_r8+cff)
#endif
# ifdef ABLM_CONVTHR
!
!-----------------------------------------------------------------------
!  Apply a convective threshold limit to Tair
!-----------------------------------------------------------------------
! 
          Tair(i,j) = MIN(ABLM_cvthres,Tair(i,j))
# endif
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
#ifdef ABLM_WATER
          CALL exchange_r2d_tile (ng, tile,                             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Hair)
#endif
        END IF

# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 4,                          &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
#ifdef ABLM_WATER
     &                    Hair,                                         &
#endif
     &                    Tair)
# endif


      RETURN
      END SUBROUTINE ablm_diabatic_tile

!
!***********************************************************************
      SUBROUTINE ablm_nudgcof (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ablm_nudgcof_tile (ng, tile, model,                          &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS)

      RETURN
      END SUBROUTINE ablm_nudgcof
!
!***********************************************************************
      SUBROUTINE ablm_nudgcof_tile (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_clima
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer ::  i, j
      real(r8) :: mw, rmw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set up nudging towards data time-scale coefficients (1/s).
!-----------------------------------------------------------------------
      mw = ABLM_maskwid
      rmw = (REAL(mw,r8)-1.0_r8)
      IF (mw.ne.1) THEN
        rmw = 1.0_r8 / rmw
      END IF

      DO j=JstrT,JendT
        DO i=IstrT,IendT
          wrk(i,j)=0.0_r8
          IF (j.gt.mw) THEN
            IF (j.lt.Mm(ng)-mw+1) THEN
              IF (i.le.mw) THEN
                wrk(i,j) = 1.0_r8 - REAL(i-1,r8)*rmw
              END IF
              IF (i.ge.Lm(ng)-mw+1) THEN
                wrk(i,j) = 1.0_r8 - REAL(Lm(ng)-i,r8)*rmw
              END IF
            ELSE
              wrk(i,j) = 1.0_r8 - REAL(Mm(ng)-j,r8)*rmw
              IF (i.le.mw) THEN
                wrk(i,j) = 1.0_r8 - REAL(i-1,r8)*rmw*                   &
     &                              REAL(Mm(ng)-j,r8)*rmw
              ELSE IF (i.ge.Lm(ng)-mw+1) THEN
                wrk(i,j) = 1.0_r8 - REAL(Lm(ng)-i,r8)*rmw*              &
     &                              REAL(Mm(ng)-j,r8)*rmw
              END IF
            END IF
          ELSE
            wrk(i,j) = 1.0_r8 - REAL(j-1,r8)*rmw
            IF (i.le.mw) THEN
               wrk(i,j) = 1.0_r8 - REAL(i-1,r8)*rmw*                    &
     &                             REAL(j-1,r8)*rmw
            ELSE IF (i.ge.Lm(ng)-mw+1) THEN
               wrk(i,j) = 1.0_r8 - REAL(Lm(ng)-i,r8)*rmw*               &
     &                             REAL(j-1,r8)*rmw
            END IF
          END IF
#ifdef MASKING
          IF (GRID(ng)%rmask(i,j).eq.0.0_r8) THEN
            wrk(i,j) = 1.0_r8
          END IF
#endif
        END DO
      END DO

      DO j=JstrT,JendT
        DO i=IstrT,IendT
          IF (wrk(i,j).eq.0.0_r8) THEN
            wrk(i,j) = ABLM_orelax
          ELSE
            wrk(i,j) = ABLM_lrelax / wrk(i,j)
          END IF
        END DO
      END DO

      DO j=JstrT,JendT
        DO i=IstrT,IendT
          IF (wrk(i,j).ne.0.0_r8) THEN
            CLIMA(ng)%ABLMnudgcof(i,j)=1.0_r8 / wrk(i,j) / 86400.0_r8
          ELSE
            CLIMA(ng)%ABLMnudgcof(i,j)=0.0_r8
          END IF
        END DO
      END DO

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchage nudging coefficients information.
!-----------------------------------------------------------------------
!
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      CLIMA(ng)%ABLMnudgcof)
#endif
      RETURN
      END SUBROUTINE ablm_nudgcof_tile
#endif
      END MODULE ablm_mod

