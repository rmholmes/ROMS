#include "cppdefs.h"
      SUBROUTINE uv3dmix2 (ng, tile)
!
!***********************************************************************
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!************************************************** Hernan G. Arango ***
!                                                                      !
!  This routine computes harmonic mixing of momentum, along constant   !
!  S-surfaces,  from the horizontal divergence of the stress tensor.   !
!  A transverse  isotropy  is  assumed so the stress tensor is split   !
!  into vertical and horizontal subtensors.                            !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!      Wajsowicz, R.C, 1993: A consistent formulation of the           !
!         anisotropic stress tensor for use in models of the           !
!         large-scale ocean circulation, JCP, 105, 333-338.            !
!                                                                      !
!      Sadourny, R. and K. Maynard, 1997: Formulations of              !
!         lateral diffusion in geophysical fluid dynamics              !
!         models, In Numerical Methods of Atmospheric and              !
!         Oceanic Modelling. Lin, Laprise, and Ritchie,                !
!         Eds., NRC Research Press, 547-556.                           !
!                                                                      !
!      Griffies, S.M. and R.W. Hallberg, 2000: Biharmonic              !
!         friction with a Smagorinsky-like viscosity for               !
!         use in large-scale eddy-permitting ocean models,             !
!         Monthly Weather Rev., 128, 8, 2935-2946.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_coupling
#ifdef DIAGNOSTICS_UV
      USE mod_diags
#endif
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 30)
#endif
      CALL uv3dmix2_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    nrhs(ng), nnew(ng),                           &
#ifdef MASKING
     &                    GRID(ng) % pmask,                             &
#endif
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % om_p,                              &
     &                    GRID(ng) % om_r,                              &
     &                    GRID(ng) % on_p,                              &
     &                    GRID(ng) % on_r,                              &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pmon_p,                            &
     &                    GRID(ng) % pmon_r,                            &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % pnom_p,                            &
     &                    GRID(ng) % pnom_r,                            &
     &                    MIXING(ng) % visc2_p,                         &
     &                    MIXING(ng) % visc2_r,                         &
#ifdef DIAGNOSTICS_UV
     &                    DIAGS(ng) % DiaRUfrc,                         &
     &                    DIAGS(ng) % DiaRVfrc,                         &
     &                    DIAGS(ng) % DiaU3wrk,                         &
     &                    DIAGS(ng) % DiaV3wrk,                         &
#endif
     &                    COUPLING(ng) % rufrc,                         &
     &                    COUPLING(ng) % rvfrc,                         &
     &                    OCEAN(ng) % u,                                &
     &                    OCEAN(ng) % v)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 30)
#endif
      RETURN
      END SUBROUTINE uv3dmix2
!
!***********************************************************************
      SUBROUTINE uv3dmix2_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          nrhs, nnew,                             &
#ifdef MASKING
     &                          pmask,                                  &
#endif
     &                          Hz,                                     &
     &                          om_p, om_r, on_p, on_r,                 &
     &                          pm, pmon_p, pmon_r,                     &
     &                          pn, pnom_p, pnom_r,                     &
     &                          visc2_p, visc2_r,                       &
#ifdef DIAGNOSTICS_UV
     &                          DiaRUfrc, DiaRVfrc,                     &
     &                          DiaU3wrk, DiaV3wrk,                     &
#endif
     &                          rufrc, rvfrc, u, v)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_p(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: on_p(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: pnom_p(LBi:,LBj:)
      real(r8), intent(in) :: pnom_r(LBi:,LBj:)
      real(r8), intent(in) :: visc2_p(LBi:,LBj:)
      real(r8), intent(in) :: visc2_r(LBi:,LBj:)

# ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRUfrc(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaRVfrc(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaU3wrk(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaV3wrk(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: rufrc(LBi:,LBj:)
      real(r8), intent(inout) :: rvfrc(LBi:,LBj:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: visc2_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: visc2_r(LBi:UBi,LBj:UBj)

# ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRUfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
      real(r8), intent(inout) :: DiaRVfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
      real(r8), intent(inout) :: DiaU3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
      real(r8), intent(inout) :: DiaV3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
# endif
      real(r8), intent(inout) :: rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rvfrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8) :: cff, cff1, cff2

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFe
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: VFe
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFx
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: VFx

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal harmonic viscosity along constant S-surfaces.
!-----------------------------------------------------------------------
!
      K_LOOP : DO k=1,N(ng)
!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.
!
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff=visc2_r(i,j)*Hz(i,j,k)*0.5_r8*                          &
     &          (pmon_r(i,j)*                                           &
     &           ((pn(i  ,j)+pn(i+1,j))*u(i+1,j,k,nrhs)-                &
     &            (pn(i-1,j)+pn(i  ,j))*u(i  ,j,k,nrhs))-               &
     &           pnom_r(i,j)*                                           &
     &           ((pm(i,j  )+pm(i,j+1))*v(i,j+1,k,nrhs)-                &
     &            (pm(i,j-1)+pm(i,j  ))*v(i,j  ,k,nrhs)))
            UFx(i,j)=on_r(i,j)*on_r(i,j)*cff
            VFe(i,j)=om_r(i,j)*om_r(i,j)*cff
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend+1
            cff=visc2_p(i,j)*0.125_r8*(Hz(i-1,j  ,k)+Hz(i,j  ,k)+       &
     &                                 Hz(i-1,j-1,k)+Hz(i,j-1,k))*      &
     &          (pmon_p(i,j)*                                           &
     &           ((pn(i  ,j-1)+pn(i  ,j))*v(i  ,j,k,nrhs)-              &
     &            (pn(i-1,j-1)+pn(i-1,j))*v(i-1,j,k,nrhs))+             &
     &           pnom_p(i,j)*                                           &
     &           ((pm(i-1,j  )+pm(i,j  ))*u(i,j  ,k,nrhs)-              &
     &            (pm(i-1,j-1)+pm(i,j-1))*u(i,j-1,k,nrhs)))
#ifdef MASKING
            cff=cff*pmask(i,j)
#endif
            UFe(i,j)=om_p(i,j)*om_p(i,j)*cff
            VFx(i,j)=on_p(i,j)*on_p(i,j)*cff
          END DO
        END DO
!
! Time-step harmonic, S-surfaces viscosity term. Notice that momentum
! at this stage is HzU and HzV and has m2/s units. Add contribution for
! barotropic forcing terms.
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.25_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
            cff1=0.5_r8*((pn(i-1,j)+pn(i,j))*                           &
     &                   (UFx(i,j  )-UFx(i-1,j))+                       &
     &                   (pm(i-1,j)+pm(i,j))*                           &
     &                   (UFe(i,j+1)-UFe(i  ,j)))
            cff2=dt(ng)*cff*cff1
            rufrc(i,j)=rufrc(i,j)+cff1
            u(i,j,k,nnew)=u(i,j,k,nnew)+cff2
#ifdef DIAGNOSTICS_UV
            DiaRUfrc(i,j,3,M2hvis)=DiaRUfrc(i,j,3,M2hvis)+cff1
            DiaU3wrk(i,j,k,M3hvis)=cff2
#endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.25_r8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
            cff1=0.5_r8*((pn(i,j-1)+pn(i,j))*                           &
     &                   (VFx(i+1,j)-VFx(i,j  ))-                       &
     &                   (pm(i,j-1)+pm(i,j))*                           &
     &                   (VFe(i  ,j)-VFe(i,j-1)))
            cff2=dt(ng)*cff*cff1
            rvfrc(i,j)=rvfrc(i,j)+cff1
            v(i,j,k,nnew)=v(i,j,k,nnew)+cff2
#ifdef DIAGNOSTICS_UV
            DiaRVfrc(i,j,3,M2hvis)=DiaRVfrc(i,j,3,M2hvis)+cff1
            DiaV3wrk(i,j,k,M3hvis)=cff2
#endif
          END DO
        END DO
      END DO K_LOOP
      RETURN
      END SUBROUTINE uv3dmix2_tile