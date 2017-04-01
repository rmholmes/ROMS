#include "cppdefs.h"
      SUBROUTINE get_data (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in forcing, climatology and other data from      !
!  NetCDF files.  If there is more than one time-record,  data is      !
!  loaded into global  two-time  record arrays. The interpolation      !
!  is carried elsewhere.                                               !
!                                                                      !
!  Currently, this routine is only executed in serial mode by the      !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
#if defined RED_TIDE
      USE mod_biology
#endif
      USE mod_boundary
      USE mod_clima
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
#if defined NLM_OUTER || defined RED_TIDE            || \
    defined TLM_CHECK || defined TL_W4DPSAS          || \
    defined W4DPSAS   || defined W4DPSAS_SENSITIVITY
      USE mod_ocean
#endif
      USE mod_scalars
      USE mod_sources
      USE mod_stepping
#ifdef SVD_ENSO
      USE mod_parallel
      USE svd_mod, ONLY : svd_regcoef
#endif

!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)

      integer :: ILB, IUB, JLB, JUB
      integer :: LBi, UBi, LBj, UBj
      integer :: i, ic, my_tile
#ifdef STATS_ENSO
      integer :: j
# ifdef SVD_ENSO
      integer :: k, tile
# endif
#endif
!
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      my_tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(my_tile)
      IUB=BOUNDS(ng)%UBi(my_tile)
      JLB=BOUNDS(ng)%LBj(my_tile)
      JUB=BOUNDS(ng)%UBj(my_tile)
!
!  Lower and upper bounds for tiled arrays.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
#endif
!
!=======================================================================
!  Read in forcing data from FORCING NetCDF file.
!=======================================================================
#if defined STATS_ENSO
# ifdef SVD_ENSO
!
!-----------------------------------------------------------------------
!  Read-in SST SVD patterns
!-----------------------------------------------------------------------
!
      CALL get_svdfld (ng, iNLM, idSSTP, ncFRCid(idSSTP,ng),            &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, stenNmod,                  &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % sstPG)
      IF (exit_flag.ne.NoError) RETURN
      IF (update(1)) stenUpdte(ng) = .TRUE.
!
!-----------------------------------------------------------------------
!  Calculate regression coefficients
!-----------------------------------------------------------------------
!
      IF (stenUpdte(ng)) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
           CALL svd_regcoef (ng, tile)
        END DO
        IF (exit_flag.ne.NoError) RETURN
      END IF
# else
!
!-----------------------------------------------------------------------
!  Mean SST index values
!-----------------------------------------------------------------------
!
        CALL get_ngfld (ng, iNLM, idSsti, ncFRCid(idSsti,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  1, 1, 1, 2, 1, 1, 1,                            &
     &                  FORCES(ng) % SSTiClimG )
        IF (exit_flag.ne.NoError) RETURN

# endif
#endif
#ifdef ADD_WWB
!
!-----------------------------------------------------------------------
!  WWB amplitude
!-----------------------------------------------------------------------
!
        CALL get_ngfld (ng, iNLM, idWwb, ncFRCid(idWwb,ng),             &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  1, 1, 1, 2, 1, 1, 1,                            &
     &                  FORCES(ng) % WWBampG )
        IF (exit_flag.ne.NoError) RETURN
#endif
#ifndef ANA_PSOURCE
!
!-----------------------------------------------------------------------
!  Point Sources/Sinks time dependent data.
!-----------------------------------------------------------------------
!
!  Point Source/Sink vertically integrated mass transport.
!
      IF (LuvSrc(ng).or.LwSrc(ng)) THEN
        CALL get_ngfld (ng, iNLM, idRtra, SSF(ng)%ncid,                 &
     &                  1, SSF(ng), update(1),                          &
     &                  1, Nsrc(ng), 1, 2, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % QbarG)
        IF (exit_flag.ne.NoError) RETURN
      END IF

# ifdef SOLVE3D
!
!  Tracer Sources/Sinks.
!
      DO i=1,NT(ng)
        IF (LtracerSrc(i,ng)) THEN
          CALL get_ngfld (ng, iNLM, idRtrc(i), SSF(ng)%ncid,            &
     &                    1, SSF(ng), update(1),                        &
     &                    1, Nsrc(ng), N(ng), 2, 1, Nsrc(ng), N(ng),    &
     &                    SOURCES(ng) % TsrcG(:,:,:,i))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
# endif
#endif

#if !defined ANA_WINDS && (defined BULK_FLUXES || defined ECOSIM)
!
!-----------------------------------------------------------------------
!  Surface wind components.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idUair, ncFRCid(idUair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
# ifdef STATS_ENSO
     &                FORCES(ng) % UwindOG)
# else
     &                FORCES(ng) % UwindG)
# endif
      IF (exit_flag.ne.NoError) RETURN
# ifdef STATS_ENSO
      IF (update(1)) stenUpdte(ng) = .TRUE.
#  ifdef SVD_ENSO
      CALL get_svdfld (ng, iNLM, idUairP, ncFRCid(idUairP,ng),          &
#  else
      CALL get_2dfld (ng, iNLM, idUairP, ncFRCid(idUairP,ng),           &
#  endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
#  ifdef SVD_ENSO
     &                LBi, UBi, LBj, UBj, 2, stenNmod,                  &
#  else
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  endif
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % UwindPG)
      IF (exit_flag.ne.NoError) RETURN
      IF (update(1)) stenUpdte(ng) = .TRUE.

      IF (stenUpdte(ng)) THEN
#  ifdef SVD_ENSO
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%UwindG(i,j,1)=FORCES(ng)%UwindOG(i,j,1)          
            FORCES(ng)%UwindG(i,j,2)=FORCES(ng)%UwindOG(i,j,2)          
            DO k=1,stenNmod
              FORCES(ng)%UwindG(i,j,1) = FORCES(ng)%UwindG(i,j,1) +     &
     &          stenRcoef(k,1,ng)*stenRuwnd(ng) *                       &
     &          FORCES(ng)%UwindPG(i,j,k,1)
              FORCES(ng)%UwindG(i,j,2) = FORCES(ng)%UwindG(i,j,2) +     &
     &          stenRcoef(k,2,ng)*stenRuwnd(ng) *                       &
     &          FORCES(ng)%UwindPG(i,j,k,2)
            END DO
          END DO
        END DO
#  else
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%UwindG(i,j,1)=FORCES(ng)%UwindOG(i,j,1)+         &
     &      stenIanom(ng)*stenRuwnd(ng)*FORCES(ng)%UwindPG(i,j,1)
            FORCES(ng)%UwindG(i,j,2)=FORCES(ng)%UwindOG(i,j,2)+         &
     &      stenIanom(ng)*stenRuwnd(ng)*FORCES(ng)%UwindPG(i,j,2)
          END DO
        END DO
#  endif
      END IF
# endif

      CALL get_2dfld (ng , iNLM, idVair, ncFRCid(idVair,ng),            &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
# ifdef STATS_ENSO
     &                FORCES(ng) % VwindOG)
# else
     &                FORCES(ng) % VwindG)
# endif
      IF (exit_flag.ne.NoError) RETURN
# ifdef STATS_ENSO
      IF (update(1)) stenUpdte(ng) = .TRUE.
#  ifdef SVD_ENSO
      CALL get_svdfld (ng, iNLM, idVairP, ncFRCid(idVairP,ng),          &
#  else
      CALL get_2dfld (ng , iNLM, idVairP, ncFRCid(idVairP,ng),          &
#  endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
#  ifdef SVD_ENSO
     &                LBi, UBi, LBj, UBj, 2, stenNmod,                  &
#  else
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  endif
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % VwindPG)
      IF (exit_flag.ne.NoError) RETURN
      IF (update(1)) stenUpdte(ng) = .TRUE.
      
      IF (stenUpdte(ng)) THEN
#  ifdef SVD_ENSO
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%VwindG(i,j,1)=FORCES(ng)%VwindOG(i,j,1)          
            FORCES(ng)%VwindG(i,j,2)=FORCES(ng)%VwindOG(i,j,2)          
            DO k=1,stenNmod
              FORCES(ng)%VwindG(i,j,1) = FORCES(ng)%VwindG(i,j,1) +     &
     &          stenRcoef(k,1,ng)*stenRvwnd(ng) *                       &
     &          FORCES(ng)%VwindPG(i,j,k,1)
              FORCES(ng)%VwindG(i,j,2) = FORCES(ng)%VwindG(i,j,2) +     &
     &          stenRcoef(k,2,ng)*stenRvwnd(ng) *                       &
     &          FORCES(ng)%VwindPG(i,j,k,2)
            END DO
          END DO
        END DO
#  else
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%VwindG(i,j,1)=FORCES(ng)%VwindOG(i,j,1)+         &
     &      stenIanom(ng)*stenRvwnd(ng)*FORCES(ng)%VwindPG(i,j,1)
            FORCES(ng)%VwindG(i,j,2)=FORCES(ng)%VwindOG(i,j,2)+         &
     &      stenIanom(ng)*stenRvwnd(ng)*FORCES(ng)%VwindPG(i,j,2)
          END DO
        END DO
#  endif
      END IF
# endif
# ifdef SET_WINDSPEED
!
!-----------------------------------------------------------------------
!  Surface wind speed magnitude.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idSair, ncFRCid(idSair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
#  ifdef STATS_ENSO
     &                FORCES(ng) % wMagOG)
#  else
     &                FORCES(ng) % wMagG)
#  endif
      IF (exit_flag.ne.NoError) RETURN
#  ifdef STATS_ENSO
      IF (update(1)) stenUpdte(ng) = .TRUE.
#  ifdef SVD_ENSO
      CALL get_svdfld (ng, iNLM, idSairP, ncFRCid(idSairP,ng),          &
#  else
      CALL get_2dfld (ng , iNLM, idSairP, ncFRCid(idSairP,ng),          &
#  endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
#   ifdef SVD_ENSO
     &                LBi, UBi, LBj, UBj, 2, stenNmod,                  &
#   else
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   endif
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % wMagPG)
      IF (exit_flag.ne.NoError) RETURN
      IF (update(1)) stenUpdte(ng) = .TRUE.
      
      IF (stenUpdte(ng)) THEN
#   ifdef SVD_ENSO
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%wMagG(i,j,1)=FORCES(ng)%wMagOG(i,j,1)            
            FORCES(ng)%wMagG(i,j,2)=FORCES(ng)%wMagOG(i,j,2)            
            DO k=1,stenNmod
              FORCES(ng)%wMagG(i,j,1) = FORCES(ng)%wMagG(i,j,1) +       &
     &          stenRcoef(k,1,ng)*stenRwmag(ng) *                       &
     &          FORCES(ng)%wMagPG(i,j,k,1)
              FORCES(ng)%wMagG(i,j,2) = FORCES(ng)%wMagG(i,j,2) +       &
     &          stenRcoef(k,2,ng)*stenRwmag(ng) *                       &
     &          FORCES(ng)%wMagPG(i,j,k,2)
            END DO
          END DO
        END DO
#   else
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%wMagG(i,j,1)=FORCES(ng)%wMagOG(i,j,1)+           &
     &      stenIanom(ng)*stenRwmag(ng)*FORCES(ng)%wMagPG(i,j,1)
            FORCES(ng)%wMagG(i,j,2)=FORCES(ng)%wMagOG(i,j,2)+           &
     &      stenIanom(ng)*stenRwmag(ng)*FORCES(ng)%wMagPG(i,j,2)
          END DO
        END DO
#   endif
      END IF
#  endif
# endif
#endif

#if defined ABLM && defined ABLM_BLHVAR
!
!-----------------------------------------------------------------------
!  Atmospheric boundary layer depth.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idHbly, ncFRCid(idHbly,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
# ifdef STATS_ENSO
     &                FORCES(ng) % blhOG)
# else
     &                FORCES(ng) % blhG)
# endif

      IF (exit_flag.ne.NoError) RETURN
# ifdef STATS_ENSO
      IF (update(1)) stenUpdte(ng) = .TRUE.
#  ifdef SVD_ENSO
      CALL get_svdfld (ng, iNLM, idHblyP, ncFRCid(idHblyP,ng),          &
#  else
      CALL get_2dfld (ng , iNLM, idHblyP, ncFRCid(idHblyP,ng),          &
#  endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
#  ifdef SVD_ENSO
     &                LBi, UBi, LBj, UBj, 2, stenNmod,                  &
#  else
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  endif
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % blhPG)
      IF (exit_flag.ne.NoError) RETURN
      IF (update(1)) stenUpdte(ng) = .TRUE.
      
      IF (stenUpdte(ng)) THEN
#   ifdef SVD_ENSO
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%blhG(i,j,1)=FORCES(ng)%blhOG(i,j,1)              
            FORCES(ng)%blhG(i,j,2)=FORCES(ng)%blhOG(i,j,2)              
            DO k=1,stenNmod
              FORCES(ng)%blhG(i,j,1) = FORCES(ng)%blhG(i,j,1) +         &
     &        stenRcoef(k,1,ng)*stenRblh(ng)*FORCES(ng)%blhPG(i,j,k,1)
              FORCES(ng)%blhG(i,j,2) = FORCES(ng)%blhG(i,j,2) +         &
     &        stenRcoef(k,2,ng)*stenRblh(ng)*FORCES(ng)%blhPG(i,j,k,2)
            END DO
          END DO
        END DO
#   else
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%blhG(i,j,1)=FORCES(ng)%blhOG(i,j,1)+             &
     &      stenIanom(ng)*stenRblh(ng)*FORCES(ng)%blhPG(i,j,1)
            FORCES(ng)%blhG(i,j,2)=FORCES(ng)%blhOG(i,j,2)+             &
     &      stenIanom(ng)*stenRblh(ng)*FORCES(ng)%blhPG(i,j,2)
          END DO
        END DO
#   endif
      END IF
# endif
#endif

#ifndef AIR_OCEAN
# if (!defined ANA_SMFLUX && !defined BULK_FLUXES) || \
       defined BULK_STRESS
!
!-----------------------------------------------------------------------
!  Surface wind stress components.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idUsms, ncFRCid(idUsms,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % umask,                                 &
#  endif
#  ifdef STATS_ENSO
     &                FORCES(ng) % sustrOG)
#  else
     &                FORCES(ng) % sustrG)
#  endif
      IF (exit_flag.ne.NoError) RETURN
#  ifdef STATS_ENSO
      IF (update(1)) stenUpdte(ng) = .TRUE.
#  ifdef SVD_ENSO
      CALL get_svdfld (ng, iNLM, idUsmsP, ncFRCid(idUsmsP,ng),          &
#  else
      CALL get_2dfld (ng, iNLM, idUsmsP, ncFRCid(idUsmsP,ng),           &
#  endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
#  ifdef SVD_ENSO
     &                LBi, UBi, LBj, UBj, 2, stenNmod,                  &
#  else
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  endif
#   ifdef MASKING
     &                GRID(ng) % umask,                                 &
#   endif
     &                FORCES(ng) % sustrPG)
      IF (exit_flag.ne.NoError) RETURN
      IF (update(1)) stenUpdte(ng) = .TRUE.
      IF (stenUpdte(ng)) THEN
#   ifdef SVD_ENSO
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%sustrG(i,j,1)=FORCES(ng)%sustrOG(i,j,1)          
            FORCES(ng)%sustrG(i,j,2)=FORCES(ng)%sustrOG(i,j,2)          
            DO k=1,stenNmod
              FORCES(ng)%sustrG(i,j,1) = FORCES(ng)%sustrG(i,j,1) +     &
     &          stenRcoef(k,1,ng)*stenRustr(ng) *                       &
     &          FORCES(ng)%sustrPG(i,j,k,1)
              FORCES(ng)%sustrG(i,j,2) = FORCES(ng)%sustrG(i,j,2) +     &
     &          stenRcoef(k,2,ng)*stenRustr(ng) *                       &
     &          FORCES(ng)%sustrPG(i,j,k,2)
            END DO
          END DO
        END DO
#   else
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%sustrG(i,j,1)=FORCES(ng)%sustrOG(i,j,1)+         &
     &      stenIanom(ng)*stenRustr(ng)*FORCES(ng)%sustrPG(i,j,1)
            FORCES(ng)%sustrG(i,j,2)=FORCES(ng)%sustrOG(i,j,2)+         &
     &      stenIanom(ng)*stenRustr(ng)*FORCES(ng)%sustrPG(i,j,2)
          END DO
        END DO
#   endif
      END IF
#  endif

      CALL get_2dfld (ng, iNLM, idVsms, ncFRCid(idVsms,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#  endif
#  ifdef STATS_ENSO
     &                FORCES(ng) % svstrOG)
#  else
     &                FORCES(ng) % svstrG)
#  endif
      IF (exit_flag.ne.NoError) RETURN
# endif
#  ifdef STATS_ENSO
      IF (update(1)) stenUpdte(ng) = .TRUE.
#  ifdef SVD_ENSO
      CALL get_svdfld (ng, iNLM, idVsmsP, ncFRCid(idVsmsP,ng),          &
#  else
      CALL get_2dfld (ng, iNLM, idVsmsP, ncFRCid(idVsmsP,ng),           &
#  endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
#  ifdef SVD_ENSO
     &                LBi, UBi, LBj, UBj, 2, stenNmod,                  &
#  else
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  endif
#   ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#   endif
     &                FORCES(ng) % svstrPG)
      IF (exit_flag.ne.NoError) RETURN
      IF (update(1)) stenUpdte(ng) = .TRUE.
      IF (stenUpdte(ng)) THEN
#   ifdef SVD_ENSO
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%svstrG(i,j,1)=FORCES(ng)%svstrOG(i,j,1)          
            FORCES(ng)%svstrG(i,j,2)=FORCES(ng)%svstrOG(i,j,2)          
            DO k=1,stenNmod
              FORCES(ng)%svstrG(i,j,1) = FORCES(ng)%svstrG(i,j,1) +     &
     &          stenRcoef(k,1,ng)*stenRvstr(ng) *                       &
     &          FORCES(ng)%svstrPG(i,j,k,1)
              FORCES(ng)%svstrG(i,j,2) = FORCES(ng)%svstrG(i,j,2) +     &
     &          stenRcoef(k,2,ng)*stenRvstr(ng) *                       &
     &          FORCES(ng)%svstrPG(i,j,k,2)
            END DO
          END DO
        END DO
#   else
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%svstrG(i,j,1)=FORCES(ng)%svstrOG(i,j,1)+         &
     &      stenIanom(ng)*stenRvstr(ng)*FORCES(ng)%svstrPG(i,j,1)
            FORCES(ng)%svstrG(i,j,2)=FORCES(ng)%svstrOG(i,j,2)+         &
     &      stenIanom(ng)*stenRvstr(ng)*FORCES(ng)%svstrPG(i,j,2)
          END DO
        END DO
#   endif
      END IF
#  endif

# if defined FOUR_DVAR && defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  If not first NLM run, get surface wind stress components from initial
!  NLM run.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
        CALL get_2dfld (ng, iNLM, idUsms, BLK(ng)%ncid,                 &
     &                  1, BLK(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
#  endif
     &                  FORCES(ng) % sustrG)
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iNLM, idVsms, BLK(ng)%ncid,                 &
     &                  1, BLK(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  FORCES(ng) % svstrG)
        IF (exit_flag.ne.NoError) RETURN
      END IF
# endif
#endif

#if !defined ANA_PAIR && (defined BULK_FLUXES || defined ECOSIM || \
                          defined ATM_PRESS)
!
!-----------------------------------------------------------------------
!  Surface air pressure.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idPair, ncFRCid(idPair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
# ifdef STATS_ENSO
     &                FORCES(ng) % PairOG)
# else
     &                FORCES(ng) % PairG)
# endif
      IF (exit_flag.ne.NoError) RETURN
# ifdef STATS_ENSO
      IF (update(1)) stenUpdte(ng) = .TRUE.
#  ifdef SVD_ENSO
      CALL get_svdfld (ng, iNLM, idPairP, ncFRCid(idPairP,ng),          &
#  else
      CALL get_2dfld (ng, iNLM, idPairP, ncFRCid(idPairP,ng),           &
#  endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
#  ifdef SVD_ENSO
     &                LBi, UBi, LBj, UBj, 2, stenNmod,                  &
#  else
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  endif
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % PairPG)
      IF (exit_flag.ne.NoError) RETURN
      IF (update(1)) stenUpdte(ng) = .TRUE.
      
      IF (stenUpdte(ng)) THEN
#   ifdef SVD_ENSO
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%PairG(i,j,1)=FORCES(ng)%PairOG(i,j,1)            
            FORCES(ng)%PairG(i,j,2)=FORCES(ng)%PairOG(i,j,2)            
            DO k=1,stenNmod
              FORCES(ng)%PairG(i,j,1) = FORCES(ng)%PairG(i,j,1) +       &
     &          stenRcoef(k,1,ng)*stenRpair(ng) *                       &
     &          FORCES(ng)%PairPG(i,j,k,1)
              FORCES(ng)%PairG(i,j,2) = FORCES(ng)%PairG(i,j,2) +       &
     &          stenRcoef(k,2,ng)*stenRpair(ng) *                       &
     &          FORCES(ng)%PairPG(i,j,k,2)
            END DO
          END DO
        END DO
#   else
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%PairG(i,j,1)=FORCES(ng)%PairOG(i,j,1)+           &
     &      stenIanom(ng)*stenRpair(ng)*FORCES(ng)%PairPG(i,j,1)
            FORCES(ng)%PairG(i,j,2)=FORCES(ng)%PairOG(i,j,2)+           &
     &      stenIanom(ng)*stenRpair(ng)*FORCES(ng)%PairPG(i,j,2)
          END DO
        END DO
#   endif
      END IF
# endif
#endif

#if !defined ANA_WWAVE && defined WAVE_DATA
!
!-----------------------------------------------------------------------
!  Surface wind induced wave amplitude, direction and period.
!-----------------------------------------------------------------------
!
# ifdef WAVES_DIR
      CALL get_2dfld (ng, iNLM, idWdir, ncFRCid(idWdir,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % DwaveG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef WAVES_HEIGHT
      CALL get_2dfld (ng, iNLM, idWamp, ncFRCid(idWamp,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % HwaveG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef WAVES_LENGTH
      CALL get_2dfld (ng, iNLM, idWlen, ncFRCid(idWlen,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % LwaveG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef WAVES_TOP_PERIOD
      CALL get_2dfld (ng, iNLM, idWptp, ncFRCid(idWptp,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Pwave_topG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef WAVES_BOT_PERIOD
      CALL get_2dfld (ng, iNLM, idWpbt, ncFRCid(idWpbt,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Pwave_botG(:,:,1))
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined WAVES_UB
      CALL get_2dfld (ng, iNLM, idWorb, ncFRCid(idWorb,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Ub_swanG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined TKE_WAVEDISS
      CALL get_2dfld (ng, iNLM, idWdis, ncFRCid(idWdis,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Wave_dissipG(:,:,1))
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined SVENDSEN_ROLLER
      CALL get_2dfld (ng, iNLM, idWbrk, ncFRCid(idWbrk,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % Wave_breakG)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#ifdef SOLVE3D

# if !defined ANA_CLOUD && defined CLOUDS
!
!-----------------------------------------------------------------------
!  Cloud fraction.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idCfra, ncFRCid(idCfra,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % cloudG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if !defined ANA_SRFLUX && defined SHORTWAVE
!
!-----------------------------------------------------------------------
!  Surface solar shortwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idSrad, ncFRCid(idSrad,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
#  ifdef STATS_ENSO
     &                FORCES(ng) % srflxOG)
#  else
     &                FORCES(ng) % srflxG)
#  endif
      IF (exit_flag.ne.NoError) RETURN
#  ifdef STATS_ENSO
      IF (update(1)) stenUpdte(ng) = .TRUE.
#  ifdef SVD_ENSO
      CALL get_svdfld (ng, iNLM, idSradP, ncFRCid(idSradP,ng),          &
#  else
      CALL get_2dfld (ng, iNLM, idSradP, ncFRCid(idSradP,ng),           &
#  endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
#  ifdef SVD_ENSO
     &                LBi, UBi, LBj, UBj, 2, stenNmod,                  &
#  else
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  endif
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % srflxPG)
      IF (exit_flag.ne.NoError) RETURN
      IF (update(1)) stenUpdte(ng) = .TRUE.
      
      IF (stenUpdte(ng)) THEN
#   ifdef SVD_ENSO
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%srflxG(i,j,1)=FORCES(ng)%srflxOG(i,j,1)          
            FORCES(ng)%srflxG(i,j,2)=FORCES(ng)%srflxOG(i,j,2)          
            DO k=1,stenNmod
              FORCES(ng)%srflxG(i,j,1) = FORCES(ng)%srflxG(i,j,1) +     &
     &          stenRcoef(k,1,ng)*stenRsrad(ng) *                       &
     &          FORCES(ng)%srflxPG(i,j,k,1)
              FORCES(ng)%srflxG(i,j,2) = FORCES(ng)%srflxG(i,j,2) +     &
     &          stenRcoef(k,2,ng)*stenRsrad(ng) *                       &
     &          FORCES(ng)%srflxPG(i,j,k,2)
            END DO
          END DO
        END DO
#   else
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%srflxG(i,j,1)=FORCES(ng)%srflxOG(i,j,1)+         &
     &      stenIanom(ng)*stenRsrad(ng)*FORCES(ng)%srflxPG(i,j,1)
            FORCES(ng)%srflxG(i,j,2)=FORCES(ng)%srflxOG(i,j,2)+         &
     &      stenIanom(ng)*stenRsrad(ng)*FORCES(ng)%srflxPG(i,j,2)
          END DO
        END DO
#   endif
      END IF
#  endif
# endif

# if defined RED_TIDE && defined DAILY_SHORTWAVE
!
!-----------------------------------------------------------------------
!  Daily-averaged Surface solar shortwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idAsrf, ncFRCid(idAsrf,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % srflxG_avg)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined BULK_FLUXES && !defined LONGWAVE && !defined LONGWAVE_OUT
!
!-----------------------------------------------------------------------
!  Surface net longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idLrad, ncFRCid(idLrad,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                      &
#  endif
     &                FORCES(ng) % lrflxG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined BULK_FLUXES && defined LONGWAVE_OUT
!
!-----------------------------------------------------------------------
!  Surface downwelling longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idLdwn, ncFRCid(idLdwn,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
#  ifdef STATS_ENSO
     &                FORCES(ng) % lrflxOG)
#  else
     &                FORCES(ng) % lrflxG)
#  endif
      IF (exit_flag.ne.NoError) RETURN
#  ifdef STATS_ENSO
      IF (update(1)) stenUpdte(ng) = .TRUE.
#  ifdef SVD_ENSO
      CALL get_svdfld (ng, iNLM, idLdwnP, ncFRCid(idLdwnP,ng),          &
#  else
      CALL get_2dfld (ng, iNLM, idLdwnP, ncFRCid(idLdwnP,ng),           &
#  endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
#  ifdef SVD_ENSO
     &                LBi, UBi, LBj, UBj, 2, stenNmod,                  &
#  else
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  endif
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % lrflxPG)
      IF (exit_flag.ne.NoError) RETURN
      IF (update(1)) stenUpdte(ng) = .TRUE.
      
      IF (stenUpdte(ng)) THEN
#   ifdef SVD_ENSO
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%lrflxG(i,j,1)=FORCES(ng)%lrflxOG(i,j,1)          
            FORCES(ng)%lrflxG(i,j,2)=FORCES(ng)%lrflxOG(i,j,2)          
            DO k=1,stenNmod
              FORCES(ng)%lrflxG(i,j,1) = FORCES(ng)%lrflxG(i,j,1) +     &
     &          stenRcoef(k,1,ng)*stenRlrad(ng) *                       &
     &          FORCES(ng)%lrflxPG(i,j,k,1)
              FORCES(ng)%lrflxG(i,j,2) = FORCES(ng)%lrflxG(i,j,2) +     &
     &          stenRcoef(k,2,ng)*stenRlrad(ng) *                       &
     &          FORCES(ng)%lrflxPG(i,j,k,2)
            END DO
          END DO
        END DO
#   else
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%lrflxG(i,j,1)=FORCES(ng)%lrflxOG(i,j,1)+         &
     &      stenIanom(ng)*stenRlrad(ng)*FORCES(ng)%lrflxPG(i,j,1)
            FORCES(ng)%lrflxG(i,j,2)=FORCES(ng)%lrflxOG(i,j,2)+         &
     &      stenIanom(ng)*stenRlrad(ng)*FORCES(ng)%lrflxPG(i,j,2)
          END DO
        END DO
#   endif
      END IF
#  endif
# endif

# if !defined ANA_TAIR && \
    ( defined BULK_FLUXES || defined ECOSIM || \
     (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO) )
!
!-----------------------------------------------------------------------
!  Surface air temperature.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTair, ncFRCid(idTair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
#  ifdef STATS_ENSO
     &                FORCES(ng) % TairOG)
#  else
     &                FORCES(ng) % TairG)
#  endif
      IF (exit_flag.ne.NoError) RETURN
#  ifdef STATS_ENSO
      IF (update(1)) stenUpdte(ng) = .TRUE.
#  ifdef SVD_ENSO
      CALL get_svdfld (ng, iNLM, idTairP, ncFRCid(idTairP,ng),          &
#  else
      CALL get_2dfld (ng, iNLM, idTairP, ncFRCid(idTairP,ng),           &
#  endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
#  ifdef SVD_ENSO
     &                LBi, UBi, LBj, UBj, 2, stenNmod,                  &
#  else
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  endif
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % TairPG)
      IF (exit_flag.ne.NoError) RETURN
      IF (update(1)) stenUpdte(ng) = .TRUE.
      
      IF (stenUpdte(ng)) THEN
#   ifdef SVD_ENSO
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%TairG(i,j,1)=FORCES(ng)%TairOG(i,j,1)            
            FORCES(ng)%TairG(i,j,2)=FORCES(ng)%TairOG(i,j,2)            
            DO k=1,stenNmod
              FORCES(ng)%TairG(i,j,1) = FORCES(ng)%TairG(i,j,1) +       &
     &          stenRcoef(k,1,ng)*stenRtair(ng) *                       &
     &          FORCES(ng)%TairPG(i,j,k,1)
              FORCES(ng)%TairG(i,j,2) = FORCES(ng)%TairG(i,j,2) +       &
     &          stenRcoef(k,2,ng)*stenRtair(ng) *                       &
     &          FORCES(ng)%TairPG(i,j,k,2)
            END DO
          END DO
        END DO
#   else
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%TairG(i,j,1)=FORCES(ng)%TairOG(i,j,1)+           &
     &      stenIanom(ng)*stenRtair(ng)*FORCES(ng)%TairPG(i,j,1)
            FORCES(ng)%TairG(i,j,2)=FORCES(ng)%TairOG(i,j,2)+           &
     &      stenIanom(ng)*stenRtair(ng)*FORCES(ng)%TairPG(i,j,2)
          END DO
        END DO
#   endif
      END IF
#  endif
# endif

# if !defined ANA_HUMIDITY && (defined BULK_FLUXES || defined ECOSIM)
!
!-----------------------------------------------------------------------
!  Surface air humidity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idQair, ncFRCid(idQair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
#  ifdef STATS_ENSO
     &                FORCES(ng) % HairOG)
#  else
     &                FORCES(ng) % HairG)
#  endif
      IF (exit_flag.ne.NoError) RETURN
#  ifdef STATS_ENSO
      IF (update(1)) stenUpdte(ng) = .TRUE.
#  ifdef SVD_ENSO
      CALL get_svdfld (ng, iNLM, idQairP, ncFRCid(idQairP,ng),          &
#  else
      CALL get_2dfld (ng, iNLM, idQairP, ncFRCid(idQairP,ng),           &
#  endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
#  ifdef SVD_ENSO
     &                LBi, UBi, LBj, UBj, 2, stenNmod,                  &
#  else
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  endif
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % HairPG)
      IF (exit_flag.ne.NoError) RETURN
      IF (update(1)) stenUpdte(ng) = .TRUE.
      
      IF (stenUpdte(ng)) THEN
#   ifdef SVD_ENSO
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%HairG(i,j,1)=FORCES(ng)%HairOG(i,j,1)            
            FORCES(ng)%HairG(i,j,2)=FORCES(ng)%HairOG(i,j,2)            
            DO k=1,stenNmod
              FORCES(ng)%HairG(i,j,1) = FORCES(ng)%HairG(i,j,1) +       &
     &          stenRcoef(k,1,ng)*stenRqair(ng) *                       &
     &          FORCES(ng)%HairPG(i,j,k,1)
              FORCES(ng)%HairG(i,j,2) = FORCES(ng)%HairG(i,j,2) +       &
     &          stenRcoef(k,2,ng)*stenRqair(ng) *                       &
     &          FORCES(ng)%HairPG(i,j,k,2)
            END DO
          END DO
        END DO
#   else
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%HairG(i,j,1)=FORCES(ng)%HairOG(i,j,1)+           &
     &      stenIanom(ng)*stenRqair(ng)*FORCES(ng)%HairPG(i,j,1)
            FORCES(ng)%HairG(i,j,2)=FORCES(ng)%HairOG(i,j,2)+           &
     &      stenIanom(ng)*stenRqair(ng)*FORCES(ng)%HairPG(i,j,2)
          END DO
        END DO
#   endif
      END IF
#  endif
# endif

# if !defined ANA_RAIN && defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Rain fall rate.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idrain, ncFRCid(idrain,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % rainG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if !defined ANA_STFLUX && !defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net heat flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTsur(itemp),                          &
     &                ncFRCid(idTsur(itemp),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
#  ifdef STATS_ENSO
     &                FORCES(ng) % stflxOG(:,:,:,itemp))
#  else
     &                FORCES(ng) % stflxG(:,:,:,itemp))
#  endif
      IF (exit_flag.ne.NoError) RETURN
#  ifdef STATS_ENSO
      IF (update(1)) stenUpdte(ng) = .TRUE.
#  ifdef SVD_ENSO
      CALL get_svdfld (ng, iNLM, idTsurP(itemp),                        &
#  else
      CALL get_2dfld (ng, iNLM, idTsurP(itemp),                         &
#  endif
     &                ncFRCid(idTsurP(itemp),ng),                       &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
#  ifdef SVD_ENSO
     &                LBi, UBi, LBj, UBj, 2, stenNmod,                  &
#  else
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  endif
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % stflxPG(:,:,:,itemp))
      IF (exit_flag.ne.NoError) RETURN
      IF (update(1)) stenUpdte(ng) = .TRUE.
      
      IF (stenUpdte(ng)) THEN
#   ifdef SVD_ENSO
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%stflxG(i,j,1,itemp)=                             &
     &                FORCES(ng)%stflxOG(i,j,1,itemp)
            FORCES(ng)%stflxG(i,j,2,itemp)=                             &
     &                FORCES(ng)%stflxOG(i,j,2,itemp)
            DO k=1,stenNmod
              FORCES(ng)%stflxG(i,j,1,itemp) =
     &            FORCES(ng)%stflxG(i,j,1,itemp) +                      &
     &            stenRcoef(k,1,ng)*stenRhflx(ng) *                     &
     &            FORCES(ng)%stflxPG(i,j,k,1,itemp)
              FORCES(ng)%stflxG(i,j,2,itemp) =
     &            FORCES(ng)%stflxG(i,j,2,itemp) +                      &
     &            stenRcoef(k,2,ng)*stenRhflx(ng) *                     &
     &            FORCES(ng)%stflxPG(i,j,k,2,itemp)
            END DO
          END DO
        END DO
#   else
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%stflxG(i,j,1,itemp)=                             &
     &          FORCES(ng)%stflxOG(i,j,1,itemp) + stenIanom(ng)*        &
     &          stenRhflx(ng)*FORCES(ng)%stflxPG(i,j,1,itemp)
            FORCES(ng)%stflxG(i,j,2,itemp)=                             &
     &          FORCES(ng)%stflxOG(i,j,2,itemp) + stenIanom(ng)*        &
     &          stenRhflx(ng)*FORCES(ng)%stflxPG(i,j,2,itemp)
          END DO
        END DO
#   endif
      END IF
#  endif
# endif

# if !defined ANA_STFLUX  && defined FOUR_DVAR      && \
      defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  If not first NLM run, get surface net heat flux from initial NLM run.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
        CALL get_2dfld (ng, iNLM, idTsur(itemp), BLK(ng)%ncid,          &
     &                  1, BLK(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#  ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#  endif
     &                  FORCES(ng) % stflxG(:,:,:,itemp))
        IF (exit_flag.ne.NoError) RETURN
      END IF
# endif

# if (!defined ANA_SST && defined QCORRECTION) || defined SVD_ENSO
!
!-----------------------------------------------------------------------
#  ifdef SVD_ENSO
!  Mean Sea Surface Temperature (SST) field for SVD_ENSO
#  else
!  Surface net heat flux correction fields: sea surface temperature
!  (SST).
#  endif
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idSSTc, ncFRCid(idSSTc,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % sstG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if !defined ANA_DQDSST && defined QCORRECTION
!
!-----------------------------------------------------------------------
!  Surface net heat flux correction fields: heat flux sensitivity to
!  SST (dQdSST).
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, iddQdT, ncFRCid(iddQdT,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % dqdtG)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifndef ANA_BTFLUX
!
!-----------------------------------------------------------------------
!  Bottom net heat flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTbot(itemp),                          &
     &                ncFRCid(idTbot(itemp),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FORCES(ng) % btflxG(:,:,:,itemp))
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef SALINITY
#  if defined NL_BULK_FLUXES && !defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux (E-P) from NLM bulk flux computation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iRPM, idEmPf, ncFRCid(idEmPf,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % stflxG(:,:,:,isalt))
      IF (exit_flag.ne.NoError) RETURN

#  elif !(defined ANA_SSFLUX || defined EMINUSP || defined SRELAXATION)
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux: E-P.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idsfwf, ncFRCid(idsfwf,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % stflxG(:,:,:,isalt))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if !defined ANA_STFLUX  && defined FOUR_DVAR      && \
       defined BULK_FLUXES && defined NL_BULK_FLUXES && \
       defined EMINUSP
!
!-----------------------------------------------------------------------
!  If not first NLM run, get freshwater flux (E-P) from initial NLM run.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
        CALL get_2dfld (ng, iNLM, idEmPf, BLK(ng)%ncid,                 &
     &                  1, BLK(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  FORCES(ng) % stflxG(:,:,:,isalt))
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  endif

#  if !defined ANA_SSS && (defined SCORRECTION || defined SRELAXATION)
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux correction field: sea surface salinity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idSSSc, ncFRCid(idSSSc,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
#   ifdef STATS_ENSO
     &                FORCES(ng) % sssOG)
#   else
     &                FORCES(ng) % sssG)
#   endif
      IF (exit_flag.ne.NoError) RETURN
#   ifdef STATS_ENSO
      IF (update(1)) stenUpdte(ng) = .TRUE.
#  ifdef SVD_ENSO
      CALL get_svdfld (ng, iNLM, idSSSP, ncFRCid(idSSSP,ng),            &
#  else
      CALL get_2dfld (ng, iNLM, idSSSP, ncFRCid(idSSSP,ng),             &
#  endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
#  ifdef SVD_ENSO
     &                LBi, UBi, LBj, UBj, 2, stenNmod,                  &
#  else
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  endif
#    ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#    endif
     &                FORCES(ng) % sssPG)
      IF (exit_flag.ne.NoError) RETURN
      IF (update(1)) stenUpdte(ng) = .TRUE.
      
      IF (stenUpdte(ng)) THEN
#   ifdef SVD_ENSO
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%sssG(i,j,1)=FORCES(ng)%sssOG(i,j,1)              
            FORCES(ng)%sssG(i,j,2)=FORCES(ng)%sssOG(i,j,2)              
            DO k=1,stenNmod
              FORCES(ng)%sssG(i,j,1) = FORCES(ng)%sssG(i,j,1) +         &
     &        stenRcoef(k,1,ng)*stenRsss(ng)*FORCES(ng)%sssPG(i,j,k,1)
              FORCES(ng)%sssG(i,j,2) = FORCES(ng)%sssG(i,j,2) +         &
     &        stenRcoef(k,2,ng)*stenRsss(ng)*FORCES(ng)%sssPG(i,j,k,2)
            END DO
          END DO
        END DO
#   else
        DO j=LBj,UBj
          DO i=LBi,UBi
            FORCES(ng)%sssG(i,j,1)=FORCES(ng)%sssOG(i,j,1)+             &
     &      stenIanom(ng)*stenRsss(ng)*FORCES(ng)%sssPG(i,j,1)
            FORCES(ng)%sssG(i,j,2)=FORCES(ng)%sssOG(i,j,2)+             &
     &      stenIanom(ng)*stenRsss(ng)*FORCES(ng)%sssPG(i,j,2)
          END DO
        END DO
#   endif
      END IF
#   endif
#  endif

#  ifndef ANA_BSFLUX
!
!-----------------------------------------------------------------------
!  Bottom net freshwater flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTbot(isalt),                          &
     &                ncFRCid(idTbot(isalt),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FORCES(ng) % btflxG(:,:,:,isalt))
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

# if defined BIOLOGY || defined SEDIMENT || defined T_PASSIVE
#  ifndef ANA_SPFLUX
!
!-----------------------------------------------------------------------
!  Passive tracers surface fluxes.
!-----------------------------------------------------------------------
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfld (ng, iNLM, idTsur(i), ncFRCid(idTsur(i),ng),     &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  FORCES(ng) % stflxG(:,:,:,i))
        IF (exit_flag.ne.NoError) RETURN
      END DO
#  endif

#  ifndef ANA_BPFLUX
!
!-----------------------------------------------------------------------
!  Passive tracers bottom fluxes.
!-----------------------------------------------------------------------
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfld (ng, iNLM, idTbot(i), ncFRCid(idTbot(i),ng),     &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  FORCES(ng) % btflxG(:,:,:,i))
        IF (exit_flag.ne.NoError) RETURN
      END DO
#  endif
# endif
#endif

#ifdef STATS_ENSO
      stenUpdte(ng) = .FALSE.
#endif
#ifdef RED_TIDE
!
!-----------------------------------------------------------------------
!  Red tide Observed Dissolved Inorganic Nutrient.
!-----------------------------------------------------------------------
!
      CALL get_3dfld (ng, iNLM, idODIN, ncFRCid(idODIN,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                OCEAN(ng) % DIN_obsG)
      IF (exit_flag.ne.NoError) RETURN
#endif
!
!=======================================================================
!  Read in open boundary conditions from BOUNDARY NetCDF file.  In
!  grid refinement, only the coarser grid (RefineScale(ng)=0) open
!  boundary conditions data is processed and needed.
!=======================================================================

#ifndef ANA_FSOBC
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(iwest), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(ieast,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(ieast), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(isouth,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(isouth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(inorth,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(inorth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
#endif

#ifndef ANA_M2OBC
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(iwest), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(iwest,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(iwest), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(ieast,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(ieast), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(ieast,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(ieast), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(isouth,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(isouth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(isouth,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(isouth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(inorth,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(inorth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(inorth,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(inorth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
#endif

#ifdef SOLVE3D
# ifndef ANA_M3OBC
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(iwest), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(iwest,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(iwest), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(ieast,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(ieast), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(ieast,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(ieast), BRY(ng)%ncid,        &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(isouth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(isouth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(isouth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(isouth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(inorth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(inorth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (LBC(inorth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(inorth), BRY(ng)%ncid,       &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
# endif

# ifndef ANA_TOBC
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        DO i=1,NT(ng)
          IF (LBC(iwest,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(iwest,i), BRY(ng)%ncid,    &
     &                      1, BRY(ng), update(1),                      &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_west(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO

        DO i=1,NT(ng)
          IF (LBC(ieast,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(ieast,i), BRY(ng)%ncid,    &
     &                      1, BRY(ng), update(1),                      &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_east(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO

        DO i=1,NT(ng)
          IF (LBC(isouth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(isouth,i), BRY(ng)%ncid,   &
     &                      1, BRY(ng), update(1),                      &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_south(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO

        DO i=1,NT(ng)
          IF (LBC(inorth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(inorth,i), BRY(ng)%ncid,   &
     &                      1, BRY(ng), update(1),                      &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_north(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
      END IF
# endif
#endif
!
!=======================================================================
!  Read in data from Climatology NetCDF file.
!=======================================================================

#ifndef ANA_SSH
!
!  Free-surface.
!
      IF (LsshCLM(ng)) THEN
        CALL get_2dfld (ng, iNLM, idSSHc, CLM(ng)%ncid,                 &
     &                  1, CLM(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
     &                  CLIMA(ng) % sshG)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#endif
#ifndef ANA_M2CLIMA
!
!  2D momentum.
!
      IF (Lm2CLM(ng)) THEN
        CALL get_2dfld (ng, iNLM, idUbcl, CLM(ng)%ncid,                   &
     &                  1, CLM(ng), update(1),                            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                  GRID(ng) % umask,                                 &
# endif
     &                  CLIMA(ng) % ubarclmG)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL get_2dfld (ng, iNLM, idVbcl, CLM(ng)%ncid,                   &
     &                  1, CLM(ng), update(1),                            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                         &
# ifdef MASKING
     &                  GRID(ng) % vmask,                                 &
# endif
     &                  CLIMA(ng) % vbarclmG)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#endif
#ifdef SOLVE3D
# ifndef ANA_M3CLIMA
!
!  3D momentum.
!
      IF (Lm3CLM(ng)) THEN
        CALL get_3dfld (ng, iNLM, idUclm, CLM(ng)%ncid,                 &
     &                  1, CLM(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
#  endif
     &                  CLIMA(ng) % uclmG)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL get_3dfld (ng, iNLM, idVclm, CLM(ng)%ncid,                 &
     &                  1, CLM(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  CLIMA(ng) % vclmG)
        IF (exit_flag.ne.NoError) RETURN
      END IF
# endif
# ifndef ANA_TCLIMA
!
!  Tracers.
!
      ic=0
      DO i=1,NT(ng)
        IF (LtracerCLM(i,ng)) THEN
          ic=ic+1
          CALL get_3dfld (ng, iNLM, idTclm(i), CLM(ng)%ncid,            &
     &                    1, CLM(ng), update(1),                        &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    CLIMA(ng) % tclmG(:,:,:,:,ic))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
# endif
#endif

#ifdef TLM_CHECK
!
!-----------------------------------------------------------------------
!  If tangent linear model check, read in nonlinear forward solution
!  to compute dot product with perturbated nonlinear solution. Time
!  interpolation between snapshot is not required (see subroutine
!  "nl_dotproduct").
!-----------------------------------------------------------------------
!
      IF (outer.ge.1) THEN
!
!  Read in free-surface.
!
        CALL get_2dfld (ng, iNLM, idFsur, FWD(ng)%ncid,                 &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
     &                  OCEAN(ng) % zetaG)
        IF (exit_flag.ne.NoError) RETURN
!
!  Read 2D momentum.
!
        CALL get_2dfld (ng, iNLM, idUbar, FWD(ng)%ncid,                 &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % umask,                               &
# endif
     &                  OCEAN(ng) % ubarG)
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iNLM, idVbar, FWD(ng)%ncid,                 &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % vmask,                               &
# endif
     &                  OCEAN(ng) % vbarG)
        IF (exit_flag.ne.NoError) RETURN
# ifdef SOLVE3D
!
!  Read in 3D momentum.
!
        CALL get_3dfld (ng, iNLM, idUvel, FWD(ng)%ncid,                 &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
#  endif
     &                  OCEAN(ng) % uG)
        IF (exit_flag.ne.NoError) RETURN

        CALL get_3dfld (ng, iNLM, idVvel, FWD(ng)%ncid,                 &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  OCEAN(ng) % vG)
        IF (exit_flag.ne.NoError) RETURN
!
!  Read in 3D tracers.
!
        DO i=1,NT(ng)
          CALL get_3dfld (ng, iNLM, idTvar(i), FWD(ng)%ncid,            &
     &                    1, FWD(ng), update(1),                        &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % tG(:,:,:,:,i))
          IF (exit_flag.ne.NoError) RETURN
        END DO
# endif
      END IF
#endif

#if defined NLM_OUTER || defined TL_W4DPSAS           || \
    defined W4DPSAS   || defined W4DPSAS_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Read weak contraint forcing snapshots. Notice that the forward
!  basic state snapshops arrays are reused here.
!-----------------------------------------------------------------------
!
      IF (FrequentImpulse(ng)) THEN
!
!  Read in free-surface.
!
        CALL get_2dfld (ng, iNLM, idFsur, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % rmask,                               &
# endif
     &                  OCEAN(ng) % zetaG)
        IF (exit_flag.ne.NoError) RETURN
!
!  Read 2D momentum.
!
        CALL get_2dfld (ng, iNLM, idUbar, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % umask,                               &
# endif
     &                  OCEAN(ng) % ubarG)
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iNLM, idVbar, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
# ifdef MASKING
     &                  GRID(ng) % vmask,                               &
# endif
     &                  OCEAN(ng) % vbarG)
        IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Read in 3D momentum.
!
        CALL get_3dfld (ng, iNLM, idUvel, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % umask,                               &
#  endif
     &                  OCEAN(ng) % uG)
        IF (exit_flag.ne.NoError) RETURN

        CALL get_3dfld (ng, iNLM, idVvel, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#  ifdef MASKING
     &                  GRID(ng) % vmask,                               &
#  endif
     &                  OCEAN(ng) % vG)
        IF (exit_flag.ne.NoError) RETURN
!
!  Read in 3D tracers.
!
        DO i=1,NT(ng)
          CALL get_3dfld (ng, iNLM, idTvar(i), TLF(ng)%ncid,            &
     &                    1, TLF(ng), update(1),                        &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % tG(:,:,:,:,i))
          IF (exit_flag.ne.NoError) RETURN
        END DO
# endif
      END IF
#endif

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3)
#endif

      RETURN
      END SUBROUTINE get_data
