      SUBROUTINE t3dmix2 (ng, tile)
!
!***********************************************************************
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!************************************************** Hernan G. Arango ***
!                                                                      !
!  This subroutine computes horizontal harmonic mixing of tracers      !
!  along geopotential surfaces.                                        !
!                                                                      !
!***********************************************************************
!
      USE mod_param
#ifdef DIAGNOSTICS_TS
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
      CALL wclock_on (ng, iNLM, 25)
#endif
      CALL t3dmix2_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   nrhs(ng), nnew(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
#endif
     &                   GRID(ng) % om_v,                               &
     &                   GRID(ng) % on_u,                               &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   MIXING(ng) % diff2,                            &
#ifdef DIAGNOSTICS_TS
     &                   DIAGS(ng) % DiaTwrk,                           &
#endif
     &                   OCEAN(ng) % t)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 25)
#endif
      RETURN
      END SUBROUTINE t3dmix2
!
!***********************************************************************
      SUBROUTINE t3dmix2_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         nrhs, nnew,                              &
#ifdef MASKING
     &                         umask, vmask,                            &
#endif
     &                         om_v, on_u, pm, pn,                      &
     &                         Hz, z_r,                                 &
     &                         diff2,                                   &
#ifdef DIAGNOSTICS_TS
     &                         DiaTwrk,                                 &
#endif
     &                         t)
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
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: diff2(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
# ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: diff2(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
# ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
     &                                   NDT)
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, itrc, j, k, k1, k2

      real(r8) :: cff, cff1, cff2, cff3, cff4

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FE
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: FS
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dTdz
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dTdx
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dTde
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dZdx
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: dZde

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal harmonic diffusion along geopotential surfaces.
!-----------------------------------------------------------------------
!
!  Compute horizontal and vertical gradients.  Notice the recursive
!  blocking sequence.  The vertical placement of the gradients is:
!
!        dTdx,dTde(:,:,k1) k     rho-points
!        dTdx,dTde(:,:,k2) k+1   rho-points
!          FS,dTdz(:,:,k1) k-1/2   W-points
!          FS,dTdz(:,:,k2) k+1/2   W-points
!
      T_LOOP : DO itrc=1,NT(ng)
        k2=1
        K_LOOP : DO k=0,N(ng)
          k1=k2
          k2=3-k1
          IF (k.lt.N(ng)) THEN
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
                cff=0.5_r8*(pm(i,j)+pm(i-1,j))
#ifdef MASKING
                cff=cff*umask(i,j)
#endif
                dZdx(i,j,k2)=cff*(z_r(i  ,j,k+1)-                       &
     &                            z_r(i-1,j,k+1))
                dTdx(i,j,k2)=cff*(t(i  ,j,k+1,nrhs,itrc)-               &
     &                            t(i-1,j,k+1,nrhs,itrc))
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
                cff=0.5_r8*(pn(i,j)+pn(i,j-1))
#ifdef MASKING
                cff=cff*vmask(i,j)
#endif
                dZde(i,j,k2)=cff*(z_r(i,j  ,k+1)-                       &
     &                            z_r(i,j-1,k+1))
                dTde(i,j,k2)=cff*(t(i,j  ,k+1,nrhs,itrc)-               &
     &                            t(i,j-1,k+1,nrhs,itrc))
              END DO
            END DO
          END IF
          IF ((k.eq.0).or.(k.eq.N(ng))) THEN
            DO j=Jstr-1,Jend+1
              DO i=Istr-1,Iend+1
                dTdz(i,j,k2)=0.0_r8
                FS(i,j,k2)=0.0_r8
              END DO
            END DO
          ELSE
            DO j=Jstr-1,Jend+1
              DO i=Istr-1,Iend+1
                cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
                dTdz(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)-                 &
     &                            t(i,j,k  ,nrhs,itrc))
              END DO
            END DO
          END IF
!
!  Compute components of the rotated tracer flux (T m3/s) along
!  geopotential surfaces.
!
          IF (k.gt.0) THEN
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
                cff=0.25_r8*(diff2(i,j,itrc)+diff2(i-1,j,itrc))*        &
     &              on_u(i,j)
                FX(i,j)=cff*                                            &
     &                  (Hz(i,j,k)+Hz(i-1,j,k))*                        &
     &                  (dTdx(i,j,k1)-                                  &
     &                   0.5_r8*(MIN(dZdx(i,j,k1),0.0_r8)*              &
     &                              (dTdz(i-1,j,k1)+                    &
     &                               dTdz(i  ,j,k2))+                   &
     &                           MAX(dZdx(i,j,k1),0.0_r8)*              &
     &                              (dTdz(i-1,j,k2)+                    &
     &                               dTdz(i  ,j,k1))))
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
                cff=0.25_r8*(diff2(i,j,itrc)+diff2(i,j-1,itrc))*        &
     &              om_v(i,j)
                FE(i,j)=cff*                                            &
     &                  (Hz(i,j,k)+Hz(i,j-1,k))*                        &
     &                  (dTde(i,j,k1)-                                  &
     &                   0.5_r8*(MIN(dZde(i,j,k1),0.0_r8)*              &
     &                              (dTdz(i,j-1,k1)+                    &
     &                               dTdz(i,j  ,k2))+                   &
     &                           MAX(dZde(i,j,k1),0.0_r8)*              &
     &                              (dTdz(i,j-1,k2)+                    &
     &                               dTdz(i,j  ,k1))))
              END DO
            END DO
            IF (k.lt.N(ng)) THEN
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  cff=0.5_r8*diff2(i,j,itrc)
                  cff1=MIN(dZdx(i  ,j,k1),0.0_r8)
                  cff2=MIN(dZdx(i+1,j,k2),0.0_r8)
                  cff3=MAX(dZdx(i  ,j,k2),0.0_r8)
                  cff4=MAX(dZdx(i+1,j,k1),0.0_r8)
                  FS(i,j,k2)=cff*                                       &
     &                       (cff1*(cff1*dTdz(i,j,k2)-dTdx(i  ,j,k1))+  &
     &                        cff2*(cff2*dTdz(i,j,k2)-dTdx(i+1,j,k2))+  &
     &                        cff3*(cff3*dTdz(i,j,k2)-dTdx(i  ,j,k2))+  &
     &                        cff4*(cff4*dTdz(i,j,k2)-dTdx(i+1,j,k1)))
                  cff1=MIN(dZde(i,j  ,k1),0.0_r8)
                  cff2=MIN(dZde(i,j+1,k2),0.0_r8)
                  cff3=MAX(dZde(i,j  ,k2),0.0_r8)
                  cff4=MAX(dZde(i,j+1,k1),0.0_r8)
                  FS(i,j,k2)=FS(i,j,k2)+                                &
     &                       cff*                                       &
     &                       (cff1*(cff1*dTdz(i,j,k2)-dTde(i,j  ,k1))+  &
     &                        cff2*(cff2*dTdz(i,j,k2)-dTde(i,j+1,k2))+  &
     &                        cff3*(cff3*dTdz(i,j,k2)-dTde(i,j  ,k2))+  &
     &                        cff4*(cff4*dTdz(i,j,k2)-dTde(i,j+1,k1)))
                END DO
              END DO
            END IF
!
!  Time-step harmonic, geopotential diffusion term (m Tunits).
!
            DO j=Jstr,Jend
              DO i=Istr,Iend
                cff=dt(ng)*pm(i,j)*pn(i,j)*                             &
     &                     (FX(i+1,j  )-FX(i,j)+                        &
     &                      FE(i  ,j+1)-FE(i,j))+                       &
     &              dt(ng)*(FS(i,j,k2)-FS(i,j,k1))
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+cff
#ifdef TS_MPDATA
                t(i,j,k,3,itrc)=t(i,j,k,nnew,itrc)
#endif
#ifdef DIAGNOSTICS_TS
                DiaTwrk(i,j,k,itrc,iThdif)=cff
#endif
              END DO
            END DO
          END IF
        END DO K_LOOP
      END DO T_LOOP
      RETURN
      END SUBROUTINE t3dmix2_tile