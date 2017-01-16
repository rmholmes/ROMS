#include "cppdefs.h"
      SUBROUTINE set_svfld (ng, model, ifield,                          &
     &                      Finp, Fout, update)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine time-interpolates requested single value from          !
!  time snapshots of input values                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     ifield     Field ID.                                             !
!     Finp       Latest two-snapshopts of value to interpolate.        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Interpolated value.                                   !
!     update     Switch indicating successful interpolation.           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: update

      integer, intent(in) :: ng, model, ifield

      real(r8), intent(in) :: Finp(2)

      real(r8), intent(out) :: Fout
!
!  Local variable declarations.
!
      logical :: Lonerec

      integer :: Tindex, it1, it2

      real(r8) :: fac, fac1, fac2
!
!----------------------------------------------------------------------
!  Set up requested field from data snapshots.
!----------------------------------------------------------------------
!
!  Get requested field information from global storage.
!
      Lonerec=Linfo(3,ifield,ng)
      Tindex=Iinfo(8,ifield,ng)
      update=.TRUE.
!
!  Set linear-interpolation factors.
!
      it1=3-Tindex
      it2=Tindex
#if defined CRAY || defined SGI
      fac1=ANINT(Tintrp(it2,ifield,ng)-time(ng))
      fac2=ANINT(time(ng)-Tintrp(it1,ifield,ng))
#else
      fac1=ANINT(Tintrp(it2,ifield,ng)-time(ng),r8)
      fac2=ANINT(time(ng)-Tintrp(it1,ifield,ng),r8)
#endif
!
!  Load time-invariant data. Time interpolation is not necessary.
!
      IF (Lonerec) THEN
         Fout=Finp(Tindex)
!
!  Time-interpolate.
!
      ELSE IF (((fac1*fac2).ge.0.0_r8).and.(fac1+fac2).gt.0.0_r8) THEN
        fac=1.0_r8/(fac1+fac2)
        fac1=fac*fac1
        fac2=fac*fac2
        Fout=fac1*Finp(it1)+fac2*Finp(it2)
!
!  Activate synchronization flag if a new time record needs to be
!  read in at the next time step.
!
        IF ((time(ng)+dt(ng)).gt.Tintrp(it2,ifield,ng)) THEN
          synchro_flag(ng)=.TRUE.
        END IF
!
!  Unable to interpolate field.  Activate error flag to quit.
!
      ELSE
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,ifield)), tdays(ng),           &
     &                      Finfo(1,ifield,ng), Finfo(2,ifield,ng),     &
     &                      Finfo(3,ifield,ng), Finfo(4,ifield,ng),     &
     &                      Tintrp(it1,ifield,ng)*sec2day,              &
     &                      Tintrp(it2,ifield,ng)*sec2day,              &
     &                      fac1*sec2day, fac2*sec2day
        END IF
  10    FORMAT (/,' SET_NGFLD  - current model time',                   &
     &          ' exceeds ending value for variable: ',a,               &
     &          /,14x,'TDAYS     = ',f15.4,                             &
     &          /,14x,'Data Tmin = ',f15.4,2x,'Data Tmax = ',f15.4,     &
     &          /,14x,'Data Tstr = ',f15.4,2x,'Data Tend = ',f15.4,     &
     &          /,14x,'TINTRP1   = ',f15.4,2x,'TINTRP2   = ',f15.4,     &
     &          /,14x,'FAC1      = ',f15.4,2x,'FAC2      = ',f15.4)
        exit_flag=2
        update=.FALSE.
      END IF
      RETURN
      END SUBROUTINE set_svfld
