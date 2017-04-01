#include "cppdefs.h"
      MODULE svd_mod
#ifdef SVD_ENSO
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the SVD regression coefficients using the SST !
!  anomaly field and input SST SVD modes.                              !
!  Note: This is called as a distributed/parallel subroutine within    !
!  the serial get_data routine                                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: svd_regcoef
!
      CONTAINS
!
!***********************************************************************      
      SUBROUTINE svd_regcoef (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
      CALL svd_regcoef_tile (ng, tile,                                  &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                GRID(ng) % omn,                                   &
     &                FORCES(ng) % ssta,                                &
     &                FORCES(ng) % sstPG)
      RETURN
      END SUBROUTINE svd_regcoef
!
!***********************************************************************
      SUBROUTINE svd_regcoef_tile (ng, tile,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      omn, ssta, sstPG)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars

!
      USE distribute_mod, ONLY : mp_reduce
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: ssta(LBi:,LBj:,:)
      real(r8), intent(in) :: sstPG(LBi:,LBj:,:,:)
# else
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ssta(LBi:UBi,LBj:UBj,:)
      real(r8), intent(in) :: sstPG(LBi:UBi,LBj:UBj,:,:)
# endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k
      integer, parameter :: Nreduce = 4
      real(r8), parameter :: small = 1.0E-05_r8
      real(r8), dimension(stenNmod,2) :: my_cov, my_var, cov, var
      real(r8), dimension(Nreduce) :: buffer
      character (len=3), dimension(Nreduce) :: op_handle

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: avg_ssta

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute averaged SST anomaly
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          avg_ssta(i,j) = 0.0_r8
          DO k=1,stenNdays
            avg_ssta(i,j) = avg_ssta(i,j) + ssta(i,j,k)
          END DO
          avg_ssta(i,j) = avg_ssta(i,j) / stenNdays
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Global sum to compute variances and covariances
!-----------------------------------------------------------------------
!
      NSUB=1
      DO k=1,stenNmod
         
        my_cov(k,1) = 0.0_r8
        my_cov(k,2) = 0.0_r8
        my_var(k,1) = 0.0_r8
        my_var(k,2) = 0.0_r8
         
        DO j=JstrR,JendR
          DO i=IstrR,IendR
             my_cov(k,1) = my_cov(k,1) + avg_ssta(i,j) *                &
     &                                  sstPG(i,j,k,1) * omn(i,j)
             my_cov(k,2) = my_cov(k,2) + avg_ssta(i,j) *                &
     &                                  sstPG(i,j,k,2) * omn(i,j)
             my_var(k,1) = my_var(k,1) + sstPG(i,j,k,1) *               &
     &                                  sstPG(i,j,k,1) * omn(i,j)
             my_var(k,2) = my_var(k,2) + sstPG(i,j,k,2) *               &
     &                                  sstPG(i,j,k,2) * omn(i,j)
          END DO
        END DO
        IF (tile_count.eq.0) THEN
          cov(k,1) = my_cov(k,1)
          cov(k,2) = my_cov(k,2)
          var(k,1) = my_var(k,1)
          var(k,2) = my_var(k,2)
        ELSE
          cov(k,1) = cov(k,1) + my_cov(k,1)
          cov(k,2) = cov(k,2) + my_cov(k,2)
          var(k,1) = var(k,1) + my_var(k,1)
          var(k,2) = var(k,2) + my_var(k,2)
        END IF
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
          buffer(1)=cov(k,1)
          buffer(2)=cov(k,2)
          buffer(3)=var(k,1)
          buffer(4)=var(k,2)
          op_handle(1)='SUM'
          op_handle(2)='SUM'
          op_handle(3)='SUM'
          op_handle(4)='SUM'
          CALL mp_reduce (ng, iNLM, Nreduce, buffer, op_handle)
          cov(k,1)=buffer(1)
          cov(k,2)=buffer(2)
          var(k,1)=buffer(3)
          var(k,2)=buffer(4)
          stenRcoef(k,1,ng) = cov(k,1)/(var(k,1) + small)
          stenRcoef(k,2,ng) = cov(k,2)/(var(k,2) + small)
        END IF
      END DO

      RETURN
      END SUBROUTINE svd_regcoef_tile
#endif
      END MODULE svd_mod

