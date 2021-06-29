!
!                                            ***********************
!                                            *    pmaxmin.f90      *
!                                            *   R. J. PURSER      *
!                                            *    January 2017     *
!                                            * jim.purser@noaa.gov *
!                                            ***********************
! Dependencies:
! Uses Modules pkind
!
!=============================================================================
! Find the extrema in (or at either end of) the second interval
! associated
! with a series of 3 or 4 data assumed uniformly spaced in some
! parameter
! (usually time). The assumption is that the data sample a smooth
! function
! of the parameter so that a suitable polynomial fitting of these data
! enables the extrema to be found as the turning points of the
! polynomial
! when they occur in the target interval, or just the interval end
! values
! if these are the interval extrema of the polynomial. For the case of 3
! data,
! say I, J, K, a quadratic fits all three and the second interval [J,K]
! may
! contain only a single extremum in its proper interior. In the case of
! 4
! data, say I, J, K, L, the linearly combined quadratics through I, J,
! K,
! and through J, K, L, imply a cubic valid as an approximation in the
! central
! interval, [J,K], so both a maximum and minumum can occur in the proper 
! interior of this target interval in principle (but note that this
! cubic only
! then collocates with J and K in general, not the outer pair, I and L).
!
! The interface name of both the 3-point and the 4-point versions, and
! for 
! either single, or double precision versions, is shared; it is
! "int2xtrema".
! In its argument list are put the given 3 or 4 consecutive values,
! followed
! by the returned max and min estimates for the second interval.
!
! In the case where ONLY the max value is required, a specialized
! version,
! "int2max" of int2xtrema is supplied; likewise, when ONLY the min value
! is
! required, subroutine "int2min" is provided.
!
! These elementary routines of the max and min kinds are combined in 
! ways that allow the interpolated max of min estimates to be obtained
! for an arbitrary number, n, of data through the use of calls to
! intnmax
! or intnmin, including versions that accommodate the cases where only
! parts
! of the array of data are "active" (in other words, the rest are taken
! to be
! "missing data"). In this case the max or min computation are done
! separately
! within the sub-sequences of consecutive data, and assessed to get the
! overall max or min from these partial evaluations. 
!=============================================================================
module pmaxmin
!=============================================================================
use pkind, only: sp,dp ! single and double precision real kinds
implicit none
private
public:: intnmax,intnmin,int2xtrema,int2max,int2min

interface intnmax
   module procedure intnmax_d
   module procedure intnmax_s
   module procedure intncmax_d 
   module procedure intncmax_s
end interface
interface intnmin
   module procedure intnmin_d
   module procedure intnmin_s
   module procedure intncmin_d
   module procedure intncmin_s
end interface
interface int2max
   module procedure int2max3_d
   module procedure int2max3_s
   module procedure int2max4_d
   module procedure int2max4_s
end interface
interface int2min
   module procedure int2min3_d
   module procedure int2min3_s
   module procedure int2min4_d
   module procedure int2min4_s
end interface
interface int2xtrema
   module procedure int2xtrema3_d
   module procedure int2xtrema3_s
   module procedure int2xtrema4_d   
   module procedure int2xtrema4_s
end interface
interface maxcubic 
   module procedure maxcubic_d
   module procedure maxcubic_s
end interface
interface mincubic 
   module procedure mincubic_d
   module procedure mincubic_s
end interface
interface maxmin
   module procedure maxmin_d
   module procedure maxmin_s
end interface 
interface cubic
   module procedure cubic_d
   module procedure cubic_s      
end interface

contains

!=============================================================================
subroutine intnmax_d(n,p,pmax)![intnmax]
!=============================================================================
! If n=1, return pmax as the single value of array p, and if n=2, return
! the
! maximum of the pair of values in p. For the cases where n>2, the
! following
! applies:
! By fitting n>=3 consecutive values, uniformly spaced, and given by
! array, p,
! using successive weighted quadratics (3-point formulae) in each
! interval,
! estimate and return the maximum value, pmax. For the end intervals the 
! single off-centered quadratic is constantly weighted, but for the
! interior
! intervals, the weighting of the two contending quadratics is linear
! across
! the interval, which ensures continuity of interpolated derivatives.
use pkind,   only: spi,dp
implicit none
integer(spi),         intent(in ):: n
real(dp),dimension(n),intent(in ):: p
real(dp),             intent(out):: pmax
!-----------------------------------------------------------------------------
real(dp)    :: qmax
integer(spi):: j
!=============================================================================
if (n .lt. 2) then 
        pmax=p(1)
elseif (n .eq. 2) then
        pmax=max(pmax,p(2))
endif
call int2max(p(3)  ,p(2)  ,p(1),pmax)
call int2max(p(n-2),p(n-1),p(n),qmax); pmax=max(pmax,qmax)
do j=1,n-3
   call int2max(p(j),p(j+1),p(j+2),p(j+3),qmax)
   pmax=max(pmax,qmax)
enddo
end subroutine intnmax_d
!=============================================================================
subroutine intnmax_s(n,p,pmax)![intnmax]
!=============================================================================
! Same as intnmax_d, except now in single precision.
use pkind,   only: spi,sp
implicit none
integer(spi),         intent(in ):: n
real(sp),dimension(n),intent(in ):: p
real(sp),             intent(out):: pmax
!-----------------------------------------------------------------------------
real(sp)    :: qmax
integer(spi):: j
!=============================================================================
if (n .lt. 2) then
        pmax=p(1) 
elseif (n .eq. 2) then
        pmax=max(pmax,p(2))
endif
call int2max(p(3)  ,p(2)  ,p(1),pmax)
call int2max(p(n-2),p(n-1),p(n),qmax); pmax=max(pmax,qmax)
do j=1,n-3
   call int2max(p(j),p(j+1),p(j+2),p(j+3),qmax)
   pmax=max(pmax,qmax)
enddo
end subroutine intnmax_s
!=============================================================================
subroutine intncmax_d(n,active,p,pmax,ff)![intnmax]
!==============================================================================
! For array of n real data, p, of which the active ones are marked
! TRUE in logical array "active", estimate the maximum of the
! values interpolated only within the sequences of consecutive
! active elements, and return it as pmax. However, if none of the
! data are active, report this failure with a returned failure flag,
! ff, set to TRUE (otherwise it becomes FALSE).
use pkind, only: spi,dp
implicit none
integer(spi),         intent(in ):: n
logical,dimension(n), intent(in ):: active
real(dp),dimension(n),intent(in ):: p
real(dp),             intent(out):: pmax
logical,              intent(out):: ff
!-----------------------------------------------------------------------------
real(dp),dimension(n):: q
real(dp)             :: qmax
integer(spi)         :: im,i,j,k,m
!=============================================================================
ff=.true.
j=0
do k=1,n ! <- k counts the distinct consecutive active sub-sequence
! search remaining elements for next active one:
   do i=j+1,n
   if (active(i)) then
        exit  
   endif
   enddo
   if (i .gt. n) then !<- failing this search, return
        im=i-1
   endif
! copy consecutive active elements from this point on into q:
   do j=i,n 
   if (.not. active(j)) then
        exit 
   endif 
        q(j-im)=p(j) 
   enddo
   m=j-i ! <- the count of this sequence of consecutive active elmts
   call intnmax(m,q,qmax)! <- find their max
   if (ff) then 
        ff=.false. 
        pmax=qmax!          <- first pmax is qmax
   else
        pmax=max(pmax,qmax)!<- another pmax
   endif
enddo
end subroutine intncmax_d
!=============================================================================
subroutine intncmax_s(n,active,p,pmax,ff)![intnmax]
!=============================================================================
! Like intncmax_d, except now in single precision.
use pkind, only: spi,sp
implicit none
integer(spi),         intent(in ):: n
logical,dimension(n), intent(in ):: active
real(sp),dimension(n),intent(in ):: p
real(sp),             intent(out):: pmax
logical,              intent(out):: ff
!------------------------------------------------------------------------------
real(sp),dimension(n):: q
real(sp)             :: qmax
integer(spi)         :: im,i,j,k,m
!==============================================================================
ff=.true.
j=0
do k=1,n ! <- k counts the distinct consecutive active sub-sequence
! search remaining elements for next active one:
   do i=j+1,n 
        if(active(i)) then
                exit 
        endif
   enddo
   if (i .gt. n) then!<- failing this search, return
        im=i-1
   endif
! copy consecutive active elements from this point on into q:
   do j=i,n; if(.not.active(j))exit; q(j-im)=p(j); enddo
   m=j-i ! <- the count of this sequence of consecutive active elmts
   call intnmax(m,q,qmax)! <- find their max
   if(ff)then; ff=.false.; pmax=qmax!          <- first pmax is qmax
   else;                   pmax=max(pmax,qmax)!<- another pmax
   endif
enddo
end subroutine intncmax_s

!=============================================================================
subroutine intnmin_d(n,p,pmin)![intnmin]
!=============================================================================
! If n=1, return pmin as the single value of array p, and if n=2, return
! the
! minimum of the pair of values in p. For the cases where n>2, the
! following
! applies:
! By fitting n>=3 consecutive values, uniformly spaced, and given by
! array, p,
! using successive weighted quadratics (3-point formulae) in each
! interval,
! estimate and return the minimum value, pmin. For the end intervals the
! single off-centered quadratic is constantly weighted, but for the
! interior
! intervals, the weighting of the two contending quadratics is linear
! across
! the interval, which ensures continuity of interpolated derivatives.
use pkind,   only: spi,dp
implicit none
integer(spi),         intent(in ):: n
real(dp),dimension(n),intent(in ):: p
real(dp),             intent(out):: pmin
!-----------------------------------------------------------------------------
real(dp)    :: qmin
integer(spi):: j
!=============================================================================
if(n<=2)then; pmin=p(1); if(n==2)pmin=min(pmin,p(2)); return; endif
call int2min(p(3)  ,p(2)  ,p(1),pmin)
call int2min(p(n-2),p(n-1),p(n),qmin); pmin=min(pmin,qmin)
do j=1,n-3
   call int2min(p(j),p(j+1),p(j+2),p(j+3),qmin); pmin=min(pmin,qmin)
enddo
end subroutine intnmin_d
!=============================================================================
subroutine intnmin_s(n,p,pmin)![intnmin]
!=============================================================================
! Same as intnmin_d, except now in single precision.
use pkind,   only: spi,sp
implicit none
integer(spi),         intent(in ):: n
real(sp),dimension(n),intent(in ):: p
real(sp),             intent(out):: pmin
!-----------------------------------------------------------------------------
real(sp)    :: qmin
integer(spi):: j
!=============================================================================
if(n<=2)then; pmin=p(1); if(n==2)pmin=min(pmin,p(2)); return; endif
call int2min(p(3)  ,p(2)  ,p(1),pmin)
call int2min(p(n-2),p(n-1),p(n),qmin); pmin=min(pmin,qmin)
do j=1,n-3
   call int2min(p(j),p(j+1),p(j+2),p(j+3),qmin); pmin=min(pmin,qmin)
enddo
end subroutine intnmin_s
!=============================================================================
subroutine intncmin_d(n,active,p,pmin,ff)![intnmin]
!=============================================================================
! For array of n real data, p, of which the active ones are marked
! TRUE in logical array "active", estimate the minimum of the
! values interpolated only within the sequences of consecutive
! active elements, and return it as pmin. However, if none of the
! data are active, report this failure with a returned failure flag,
! ff, set to TRUE (otherwise it becomes FALSE).
use pkind, only: spi,dp
implicit none
integer(spi),         intent(in ):: n
logical,dimension(n), intent(in ):: active
real(dp),dimension(n),intent(in ):: p
real(dp),             intent(out):: pmin
logical,              intent(out):: ff
!-----------------------------------------------------------------------------
real(dp),dimension(n):: q
real(dp)             :: qmin
integer(spi)         :: im,i,j,k,m
!=============================================================================
ff=.true.
j=0
do k=1,n ! <- k counts the distinct consecutive active sub-sequence
! search remaining elements for next active one:
   do i=j+1,n; if(active(i))exit; enddo
   if(i>n)return !<- failing this search, return
   im=i-1
! copy consecutive active elements from this point on into q:
   do j=i,n; if(.not.active(j))exit; q(j-im)=p(j); enddo
   m=j-i ! <- the count of this sequence of consecutive active elmts
   call intnmin(m,q,qmin)! <- find their min
   if(ff)then; ff=.false.; pmin=qmin!          <- first pmin is qmin
   else;                   pmin=min(pmin,qmin)!<- another pmin
   endif
enddo
end subroutine intncmin_d
!=============================================================================
subroutine intncmin_s(n,active,p,pmin,ff)![intnmin]
!=============================================================================
! Like intncmin_d, except now in single precision.
use pkind, only: spi,sp
implicit none
integer(spi),         intent(in ):: n
logical,dimension(n), intent(in ):: active
real(sp),dimension(n),intent(in ):: p
real(sp),             intent(out):: pmin
logical,              intent(out):: ff
!-----------------------------------------------------------------------------
real(sp),dimension(n):: q
real(sp)             :: qmin
integer(spi)         :: im,i,j,k,m
!=============================================================================
ff=.true.
j=0
do k=1,n ! <- k counts the distinct consecutive active sub-sequence
! search remaining elements for next active one:
   do i=j+1,n; if(active(i))exit; enddo
   if(i>n)return !<- failing this search, return
   im=i-1
! copy consecutive active elements from this point on into q:
   do j=i,n; if(.not.active(j))exit; q(j-im)=p(j); enddo
   m=j-i ! <- the count of this sequence of consecutive active elmts
   call intnmin(m,q,qmin)! <- find their min
   if(ff)then; ff=.false.; pmin=qmin!          <- first pmin is qmin
   else;                   pmin=min(pmin,qmin)!<- another pmin
   endif
enddo
end subroutine intncmin_s

!=============================================================================
subroutine int2max4_d(vi,vj,vk,vl, fmax)![int2max]
!=============================================================================
! Find the maximum value of the function interpolated by weighted
! quadratics
! within interval-2 of the 3 equal intervals whose consecutive nodes
! have
! values vi, vj, vj, vl. The maximum value in this central interval is
! fmax.
use pkind, only: dp
implicit none
real(dp),intent(in ):: vi,vj,vk,vl
real(dp),intent(out):: fmax
!-----------------------------------------------------------------------------
real(dp),parameter:: u1=1._dp,u3=3._dp,u9=9._dp,u11=11._dp,o16=u1/16._dp
real(dp)          :: a,b,c,d,vipl,viml,vjpk,vjmk
!============================================================================
vipl=vi+vl; viml=vi-vl ! <- Outer pair sums and differences
vjpk=vj+vk; vjmk=vj-vk ! <- Inner pair sums and differences
a=(-viml+ u3*vjmk)*o16 ! = (-vi + 3*vj - 3*vk +vl)/16
b=( vipl-    vjpk)*o16 ! = ( vi -   vj -   vk +vl)/16
c=( viml-u11*vjmk)*o16 ! = ( vi -11*vj +11*vk -vl)/16
d=(-vipl+ u9*vjpk)*o16 ! = (-vi + 9*vj + 9*vk -vl)/16
call maxcubic(a,b,c,d,vj,vk, fmax)
end subroutine int2max4_d
!=============================================================================
subroutine int2max4_s(vi,vj,vk,vl, fmax)![int2max]
!=============================================================================
! Like int2max4_d, except using single precision reals
use pkind, only: sp
implicit none
real(sp),intent(in ):: vi,vj,vk,vl
real(sp),intent(out):: fmax
!-----------------------------------------------------------------------------
real(sp),parameter:: u1=1._sp,u3=3._sp,u9=9._sp,u11=11._sp,o16=u1/16._sp
real(sp)          :: a,b,c,d,vipl,viml,vjpk,vjmk
!============================================================================
vipl=vi+vl; viml=vi-vl ! <- Outer pair sums and differences
vjpk=vj+vk; vjmk=vj-vk ! <- Inner pair sums and differences
a=(-viml+ u3*vjmk)*o16 ! = (-vi + 3*vj - 3*vk +vl)/16
b=( vipl-    vjpk)*o16 ! = ( vi -   vj -   vk +vl)/16
c=( viml-u11*vjmk)*o16 ! = ( vi -11*vj +11*vk -vl)/16
d=(-vipl+ u9*vjpk)*o16 ! = (-vi + 9*vj + 9*vk -vl)/16
call maxcubic(a,b,c,d,vj,vk, fmax)
end subroutine int2max4_s
!=============================================================================
subroutine int2max3_d(vi,vj,vk, fmax)![int2max]
!=============================================================================
! Find the maximum value of the function interpolated by a quadratic
! within the interval-1 of the 2 equal intervals whose consecutive nodes
! have values vi, vj, vk. The maximum value in this second interval is
! fmax.
use pkind, only: dp
implicit none
real(dp),intent(in ):: vi,vj,vk
real(dp),intent(out):: fmax
!-----------------------------------------------------------------------------
real(dp),parameter:: a=0._dp,u1=1._dp,u3=3._dp,o2=u1/2._dp,o8=u1/8._dp
real(dp)          :: b,c,d,vimj,vjmk
!=============================================================================
vimj=vi-vj
vjmk=vj-vk
b=   (vimj   -vjmk)*o8       ! = ( vi -2*vj +  vk)/8
c=           -vjmk *o2       ! = (    -4*vj +4*vk)/8
d=vj-(vimj+u3*vjmk)*o8       ! = (-vi +6*vj +3*vk)/8
call maxcubic(a,b,c,d,vj,vk, fmax)
end subroutine int2max3_d
!=============================================================================
subroutine int2max3_s(vi,vj,vk, fmax)![int2max]
!=============================================================================
! Like int2max3_d, except using single precision reals.
use pkind, only: sp
implicit none
real(sp),intent(in ):: vi,vj,vk
real(sp),intent(out):: fmax
!-----------------------------------------------------------------------------
real(sp),parameter:: a=0._sp,u1=1._sp,u3=3._sp,o2=u1/2._sp,o8=u1/8._sp
real(sp)          :: b,c,d,vimj,vjmk
!=============================================================================
vimj=vi-vj
vjmk=vj-vk
b=   (vimj   -vjmk)*o8       ! = ( vi -2*vj +  vk)/8
c=           -vjmk *o2       ! = (    -4*vj +4*vk)/8
d=vj-(vimj+u3*vjmk)*o8       ! = (-vi +6*vj +3*vk)/8
call maxcubic(a,b,c,d,vj,vk, fmax)
end subroutine int2max3_s

!=============================================================================
subroutine int2min4_d(vi,vj,vk,vl, fmin)![int2min]
!=============================================================================
! Find the minimum value of the function interpolated by weighted
! quadratics
! within interval-2 of the 3 equal intervals whose consecutive nodes
! have
! values vi, vj, vj, vl. The minimum value in this central interval is
! fmin.
use pkind, only: dp
implicit none
real(dp),intent(in ):: vi,vj,vk,vl
real(dp),intent(out):: fmin
!-----------------------------------------------------------------------------
real(dp),parameter:: u1=1._dp,u3=3._dp,u9=9._dp,u11=11._dp,o16=u1/16._dp
real(dp)          :: a,b,c,d,vipl,viml,vjpk,vjmk
!============================================================================
vipl=vi+vl; viml=vi-vl ! <- Outer pair sums and differences
vjpk=vj+vk; vjmk=vj-vk ! <- Inner pair sums and differences
a=(-viml+ u3*vjmk)*o16 ! = (-vi + 3*vj - 3*vk +vl)/16
b=( vipl-    vjpk)*o16 ! = ( vi -   vj -   vk +vl)/16
c=( viml-u11*vjmk)*o16 ! = ( vi -11*vj +11*vk -vl)/16
d=(-vipl+ u9*vjpk)*o16 ! = (-vi + 9*vj + 9*vk -vl)/16
call mincubic(a,b,c,d,vj,vk, fmin)
end subroutine int2min4_d
!=============================================================================
subroutine int2min4_s(vi,vj,vk,vl, fmin)![int2min]
!=============================================================================
! Like int2min4_d, except using single precision reals
use pkind, only: sp
implicit none
real(sp),intent(in ):: vi,vj,vk,vl
real(sp),intent(out):: fmin
!-----------------------------------------------------------------------------
real(sp),parameter:: u1=1._sp,u3=3._sp,u9=9._sp,u11=11._sp,o16=u1/16._sp
real(sp)          :: a,b,c,d,vipl,viml,vjpk,vjmk
!============================================================================
vipl=vi+vl; viml=vi-vl ! <- Outer pair sums and differences
vjpk=vj+vk; vjmk=vj-vk ! <- Inner pair sums and differences
a=(-viml+ u3*vjmk)*o16 ! = (-vi + 3*vj - 3*vk +vl)/16
b=( vipl-    vjpk)*o16 ! = ( vi -   vj -   vk +vl)/16
c=( viml-u11*vjmk)*o16 ! = ( vi -11*vj +11*vk -vl)/16
d=(-vipl+ u9*vjpk)*o16 ! = (-vi + 9*vj + 9*vk -vl)/16
call mincubic(a,b,c,d,vj,vk, fmin)
end subroutine int2min4_s
!=============================================================================
subroutine int2min3_d(vi,vj,vk, fmin)![int2min]
!=============================================================================
! Find the minimum value of the function interpolated by a quadratic
! within the interval-1 of the 2 equal intervals whose consecutive nodes
! have values vi, vj, vk. The minimum value in this second interval is
! fmin.
use pkind, only: dp
implicit none
real(dp),intent(in ):: vi,vj,vk
real(dp),intent(out):: fmin
!-----------------------------------------------------------------------------
real(dp),parameter:: a=0._dp,u1=1._dp,u3=3._dp,o2=u1/2._dp,o8=u1/8._dp
real(dp)          :: b,c,d,vimj,vjmk
!=============================================================================
vimj=vi-vj
vjmk=vj-vk
b=   (vimj   -vjmk)*o8       ! = ( vi -2*vj +  vk)/8
c=           -vjmk *o2       ! = (    -4*vj +4*vk)/8
d=vj-(vimj+u3*vjmk)*o8       ! = (-vi +6*vj +3*vk)/8
call mincubic(a,b,c,d,vj,vk, fmin)
end subroutine int2min3_d
!=============================================================================
subroutine int2min3_s(vi,vj,vk, fmin)![int2min]
!=============================================================================
! Like int2min3_d, except using single precision reals.
use pkind, only: sp
implicit none
real(sp),intent(in ):: vi,vj,vk
real(sp),intent(out):: fmin
!-----------------------------------------------------------------------------
real(sp),parameter:: a=0._sp,u1=1._sp,u3=3._sp,o2=u1/2._sp,o8=u1/8._sp
real(sp)          :: b,c,d,vimj,vjmk
!=============================================================================
vimj=vi-vj
vjmk=vj-vk
b=   (vimj   -vjmk)*o8       ! = ( vi -2*vj +  vk)/8
c=           -vjmk *o2       ! = (    -4*vj +4*vk)/8
d=vj-(vimj+u3*vjmk)*o8       ! = (-vi +6*vj +3*vk)/8
call mincubic(a,b,c,d,vj,vk, fmin)
end subroutine int2min3_s

!=============================================================================
subroutine int2xtrema4_d(vi,vj,vk,vl, fmax,fmin)![int2xtrema]
!=============================================================================
! From a uniformly-spaced series of 4 values, vi, vj, vk, vl, estimate
! the
! max and min (fmax, fmin) in interval-two, i.e., in [j,k], by linearly
! weighting, over this interval, the quadratic fit through I,J,K and the
! quadratic fit through J,K,L. This assumed cubic, valid in the central
! interval, does not collocate the end values, I and L, but it does
! ensure
! continuity of derivatives (and hence, the locations of extrema) when
! the
! same procedure (or the deficient 3-point version, int2xtrema3) is used
! for the neighboring intervals as well.
! The cubic polynomial is defined, f(x) = a*x**3 + b*x**2 + c*x + d,
! where the coordinate x is scaled and positioned to put (I,J,K,L) at
! (-3,-1,+1,+3), in order to facilitate efficient computation. We invoke
! subr. maxmin to examine whether extrema interior to the target
! interval
! (J,K) exist to override default estimates that are given by the
! interval's
! bounding values vj or vj themselves.
use pkind, only: dp
implicit none
real(dp),intent(in ):: vi,vj,vk,vl
real(dp),intent(out):: fmax,fmin
!-----------------------------------------------------------------------------
real(dp),parameter:: u1=1._dp,u3=3._dp,u9=9._dp,u11=11._dp,o16=u1/16._dp
real(dp)          :: a,b,c,d,vipl,viml,vjpk,vjmk
!============================================================================
vipl=vi+vl; viml=vi-vl ! <- Outer pair sums and differences
vjpk=vj+vk; vjmk=vj-vk ! <- Inner pair sums and differences
a=(-viml+ u3*vjmk)*o16 ! = (-vi + 3*vj - 3*vk +vl)/16
b=( vipl-    vjpk)*o16 ! = ( vi -   vj -   vk +vl)/16
c=( viml-u11*vjmk)*o16 ! = ( vi -11*vj +11*vk -vl)/16
d=(-vipl+ u9*vjpk)*o16 ! = (-vi + 9*vj + 9*vk -vl)/16
call maxmin(a,b,c,d,vj,vk, fmax,fmin)
end subroutine int2xtrema4_d
!=============================================================================
subroutine int2xtrema4_s(vi,vj,vk,vl, fmax,fmin)![int2xtrema]
!=============================================================================
! Like int2xtrema4_d, except using single precision reals
use pkind, only: sp
implicit none
real(sp),intent(in ):: vi,vj,vk,vl
real(sp),intent(out):: fmax,fmin
!-----------------------------------------------------------------------------
real(sp),parameter:: u1=1._sp,u3=3._sp,u9=9._sp,u11=11._sp,o16=u1/16._sp
real(sp)          :: a,b,c,d,vipl,viml,vjpk,vjmk
!============================================================================
vipl=vi+vl; viml=vi-vl ! <- Outer pair sums and differences
vjpk=vj+vk; vjmk=vj-vk ! <- Inner pair sums and differences
a=(-viml+ u3*vjmk)*o16 ! = (-vi + 3*vj - 3*vk +vl)/16
b=( vipl-    vjpk)*o16 ! = ( vi -   vj -   vk +vl)/16
c=( viml-u11*vjmk)*o16 ! = ( vi -11*vj +11*vk -vl)/16
d=(-vipl+ u9*vjpk)*o16 ! = (-vi + 9*vj + 9*vk -vl)/16
call maxmin(a,b,c,d,vj,vk, fmax,fmin)
end subroutine int2xtrema4_s
!=============================================================================
subroutine int2xtrema3_d(vi,vj,vk, fmax,fmin)![int2xtrema]
!=============================================================================
! From a uniformly-spaced series of 3 values, vi, vj, vk, estimate the
! max and min (fmax, fmin) in interval-two, i.e., in [j,k], by fitting a
! quadratic through the 3 values. This is a degenerate cubic of the form
! used in int2xtrema4, but with a=0, so we can still invoke maxmin for
! the search of possible extrema interior to the target interval, (J,K).
use pkind, only: dp
implicit none
real(dp),intent(in ):: vi,vj,vk
real(dp),intent(out):: fmax,fmin
!-----------------------------------------------------------------------------
real(dp),parameter:: a=0._dp,u1=1._dp,u3=3._dp,o2=u1/2._dp,o8=u1/8._dp
real(dp)          :: b,c,d,vimj,vjmk
!=============================================================================
vimj=vi-vj
vjmk=vj-vk
b=   (vimj   -vjmk)*o8       ! = ( vi -2*vj +  vk)/8
c=           -vjmk *o2       ! = (    -4*vj +4*vk)/8 
d=vj-(vimj+u3*vjmk)*o8       ! = (-vi +6*vj +3*vk)/8
call maxmin(a,b,c,d,vj,vk, fmax,fmin)
end subroutine int2xtrema3_d
!=============================================================================
subroutine int2xtrema3_s(vi,vj,vk, fmax,fmin)![int2xtrema]
!=============================================================================
! Like intxtrema3_d, except using single precision reals.
use pkind, only: sp
implicit none
real(sp),intent(in ):: vi,vj,vk
real(sp),intent(out):: fmax,fmin
!-----------------------------------------------------------------------------
real(sp),parameter:: a=0._sp,u1=1._sp,u3=3._sp,o2=u1/2._sp,o8=u1/8._sp
real(sp)          :: b,c,d,vimj,vjmk
!=============================================================================
vimj=vi-vj
vjmk=vj-vk
b=   (vimj   -vjmk)*o8       ! = ( vi -2*vj +  vk)/8
c=           -vjmk *o2       ! = (    -4*vj +4*vk)/8 
d=vj-(vimj+u3*vjmk)*o8       ! = (-vi +6*vj +3*vk)/8
call maxmin(a,b,c,d,vj,vk, fmax,fmin)
end subroutine int2xtrema3_s

!=============================================================================
subroutine maxcubic_d(a,b,c,d,vj,vk, fmax)![maxcubic]
!=============================================================================
! Find the maximum of the cubic with coefficients a, b, c, d, within the
! interval [-1,+1] at whose end points the values are vj and vk. The
! returned
! maximum is fmax.
use pkind, only: dp
implicit none
real(dp),intent(in ):: a,b,c,d   ! Cubic coefficients
real(dp),intent(in ):: vj,vk     ! End values of 2nd interval
real(dp),intent(out):: fmax      ! Estimates of max and min
!-----------------------------------------------------------------------------
real(dp),parameter:: u0=0._dp,u1=1._dp,u2=2._dp,u3=3._dp
real(dp)          :: a3i,b2,bq,cq,x1,x2,r,fmin
!=============================================================================
fmin=u0 ! Dummy value given so that subroutine cubic will not fail
fmax=max(vj,vk) ! End values are the starting default estimates
if(a==u0)then   ! The special degenerate case of a quadratic polynomial:
   b2=b*u2
   if(abs(c)<abs(b2))call cubic(a,b,c,d,-c/b2,fmax,fmin)
else
! Express the quadratic equation for the extrema in the convenient form:
! x**2 - 2*bq*x + cq = 0:
   a3i=u1/(a*u3)
   bq=-b*a3i
   cq= c*a3i
! The discriminant of this quadratic is r; if r>0 we can take its sqrt
! and find two real and distinct roots (otherwise, we don't proceed):
   r=bq*bq-cq
   if(r>u0)then
      r=sqrt(r)
      if(bq/=u0)then ! Roots are distinct
         if(bq<u0)r=-r ! Choose sign of r to be same as that of bq
         x2=bq+r   ! <- The farther root (from x=0) of the quadratic
         x1=cq/x2  ! <- The closer root (from x=0) of the quadratic
         if(abs(x1)<u1)then ! Extremum x1 occurs inside the target
            call cubic(a,b,c,d,x1, fmax,fmin)
!              .. so check whether even this farther one, x2,  does
!              also:
            if(abs(x2)<u1)call cubic(a,b,c,d,x2, fmax,fmin)
         endif
      endif
   endif
endif
end subroutine maxcubic_d
!=============================================================================
subroutine maxcubic_s(a,b,c,d,vj,vk, fmax)![maxcubic]
!=============================================================================
! Like maxcubic_d, except using single precision reals
use pkind, only: sp
implicit none
real(sp),intent(in ):: a,b,c,d   ! Cubic coefficients
real(sp),intent(in ):: vj,vk     ! End values of 2nd interval
real(sp),intent(out):: fmax      ! Estimates of max and min
!-----------------------------------------------------------------------------
real(sp),parameter:: u0=0._sp,u1=1._sp,u2=2._sp,u3=3._sp
real(sp)          :: a3i,b2,bq,cq,x1,x2,r,fmin
!=============================================================================
fmin=u0 ! Dummy value given so that subroutine cubic will not fail
fmax=max(vj,vk) ! End values are the starting default estimates
if(a==u0)then   ! The special degenerate case of a quadratic polynomial:
   b2=b*u2
   if(abs(c)<abs(b2))call cubic(a,b,c,d,-c/b2,fmax,fmin)
else
! Express the quadratic equation for the extrema in the convenient form:
! x**2 - 2*bq*x + cq = 0:
   a3i=u1/(a*u3)
   bq=-b*a3i
   cq= c*a3i
! The discriminant of this quadratic is r; if r>0 we can take its sqrt
! and find two real and distinct roots (otherwise, we don't proceed):
   r=bq*bq-cq
   if(r>u0)then
      r=sqrt(r)
      if(bq/=u0)then ! Roots are distinct
         if(bq<u0)r=-r ! Choose sign of r to be same as that of bq
         x2=bq+r   ! <- The farther root (from x=0) of the quadratic
         x1=cq/x2  ! <- The closer root (from x=0) of the quadratic
         if(abs(x1)<u1)then ! Extremum x1 occurs inside the target
            call cubic(a,b,c,d,x1, fmax,fmin)
!              .. so check whether even this farther one, x2,  does
!              also:
            if(abs(x2)<u1)call cubic(a,b,c,d,x2, fmax,fmin)
         endif
      endif
   endif
endif
end subroutine maxcubic_s

!=============================================================================
subroutine mincubic_d(a,b,c,d,vj,vk, fmin)![mincubic]
!=============================================================================
! Find the maximum of the cubic with coefficients a, b, c, d, within the
! interval [-1,+1] at whose end points the values are vj and vk. The
! returned
! maximum is fmax.
use pkind, only: dp
implicit none
real(dp),intent(in ):: a,b,c,d   ! Cubic coefficients
real(dp),intent(in ):: vj,vk     ! End values of 2nd interval
real(dp),intent(out):: fmin      ! Estimates of max and min
!-----------------------------------------------------------------------------
real(dp),parameter:: u0=0._dp,u1=1._dp,u2=2._dp,u3=3._dp
real(dp)          :: a3i,b2,bq,cq,x1,x2,r,fmax
!=============================================================================
fmax=u0 ! Dummy value given so that subroutine cubic will not fail
fmin=min(vj,vk) ! End values are the starting default estimates
if(a==u0)then   ! The special degenerate case of a quadratic polynomial:
   b2=b*u2
   if(abs(c)<abs(b2))call cubic(a,b,c,d,-c/b2,fmax,fmin)
else
! Express the quadratic equation for the extrema in the convenient form:
! x**2 - 2*bq*x + cq = 0:
   a3i=u1/(a*u3)
   bq=-b*a3i
   cq= c*a3i
! The discriminant of this quadratic is r; if r>0 we can take its sqrt
! and find two real and distinct roots (otherwise, we don't proceed):
   r=bq*bq-cq
   if(r>u0)then
      r=sqrt(r)
      if(bq/=u0)then ! Roots are distinct
         if(bq<u0)r=-r ! Choose sign of r to be same as that of bq
         x2=bq+r   ! <- The farther root (from x=0) of the quadratic
         x1=cq/x2  ! <- The closer root (from x=0) of the quadratic
         if(abs(x1)<u1)then ! Extremum x1 occurs inside the target
            call cubic(a,b,c,d,x1, fmax,fmin)
!              .. so check whether even this farther one, x2,  does
!              also:
            if(abs(x2)<u1)call cubic(a,b,c,d,x2, fmax,fmin)
         endif
      endif
   endif
endif
end subroutine mincubic_d
!=============================================================================
subroutine mincubic_s(a,b,c,d,vj,vk, fmin)![mincubic]
!=============================================================================
! Like mincubic_d, except using single precision reals
use pkind, only: sp
implicit none
real(sp),intent(in ):: a,b,c,d   ! Cubic coefficients
real(sp),intent(in ):: vj,vk     ! End values of 2nd interval
real(sp),intent(out):: fmin      ! Estimates of max and min
!-----------------------------------------------------------------------------
real(sp),parameter:: u0=0._sp,u1=1._sp,u2=2._sp,u3=3._sp
real(sp)          :: a3i,b2,bq,cq,x1,x2,r,fmax
!=============================================================================
fmax=u0 ! Dummy value given so that subroutine cubic will not fail
fmin=min(vj,vk) ! End values are the starting default estimates
if(a==u0)then   ! The special degenerate case of a quadratic polynomial:
   b2=b*u2
   if(abs(c)<abs(b2))call cubic(a,b,c,d,-c/b2,fmax,fmin)
else
! Express the quadratic equation for the extrema in the convenient form:
! x**2 - 2*bq*x + cq = 0:
   a3i=u1/(a*u3)
   bq=-b*a3i
   cq= c*a3i
! The discriminant of this quadratic is r; if r>0 we can take its sqrt
! and find two real and distinct roots (otherwise, we don't proceed):
   r=bq*bq-cq
   if(r>u0)then
      r=sqrt(r)
      if(bq/=u0)then ! Roots are distinct
         if(bq<u0)r=-r ! Choose sign of r to be same as that of bq
         x2=bq+r   ! <- The farther root (from x=0) of the quadratic
         x1=cq/x2  ! <- The closer root (from x=0) of the quadratic
         if(abs(x1)<u1)then ! Extremum x1 occurs inside the target
            call cubic(a,b,c,d,x1, fmax,fmin)
!              .. so check whether even this farther one, x2,  does
!              also:
            if(abs(x2)<u1)call cubic(a,b,c,d,x2, fmax,fmin)
         endif
      endif
   endif
endif
end subroutine mincubic_s

!=============================================================================
subroutine maxmin_d(a,b,c,d,vj,vk, fmax,fmin)![maxmin]
!=============================================================================
! Estimate the max and the min, fmax and fmin, of the cubic if they seem
! to lie
! within the given interval [-1,+1] or else take these estimates from
! the
! given values, vj or vk, at x=-1 and x=+1.
!
! The possibly degenerate cubic is assumed to be of the form,
! f(x) = a*x**3 + b*x**2 + c*x + d.
!
! If a==0, then the deficient cubic reduces to a quadratic with only one
! extremum to be examined. But in the generic case, two extrema occur as
! the roots of the quadratic obtained by differentiating the cubic. Only
! when
! the discriminant of the quadratic is positive do we need to examine
! these 
! roots (the special case where they coincide implies a stationary 
! inflexion, which does not interest us, and when the discriminant is
! negative, the roots are complex and the cubic has no real extrema at
! all).
use pkind, only: dp
implicit none
real(dp),intent(in ):: a,b,c,d   ! Cubic coefficients
real(dp),intent(in ):: vj,vk     ! End values of target interval
real(dp),intent(out):: fmax,fmin ! Estimates of max and min
!-----------------------------------------------------------------------------
real(dp),parameter:: u0=0._dp,u1=1._dp,u2=2._dp,u3=3._dp
real(dp)          :: a3i,b2,bq,cq,x1,x2,r
!=============================================================================
fmax=max(vj,vk); fmin=min(vj,vk)! End values are the starting default
if(a==u0)then ! The special degenerate case of a quadratic polynomial:
   b2=b*u2
   if(abs(c)<abs(b2))call cubic(a,b,c,d,-c/b2,fmax,fmin)
else
! Express the quadratic equation for the extrema in the convenient form:
! x**2 - 2*bq*x + cq = 0:
   a3i=u1/(a*u3)
   bq=-b*a3i
   cq= c*a3i
! The discriminant of this quadratic is r; if r>0 we can take its sqrt
! and find two real and distinct roots (otherwise, we don't proceed):
   r=bq*bq-cq
   if(r>u0)then
      r=sqrt(r)
      if(bq/=u0)then ! Roots are distinct
         if(bq<u0)r=-r ! Choose sign of r to be same as that of bq
         x2=bq+r   ! <- The farther root (from x=0) of the quadratic
         x1=cq/x2  ! <- The closer root (from x=0) of the quadratic
         if(abs(x1)<u1)then ! Extremum x1 occurs inside the target
            call cubic(a,b,c,d,x1, fmax,fmin)
!              .. so check whether even this farther one, x2,  does
!              also:
            if(abs(x2)<u1)call cubic(a,b,c,d,x2, fmax,fmin)
         endif
      endif
   endif
endif
end subroutine maxmin_d
!=============================================================================
subroutine maxmin_s(a,b,c,d,vj,vk, fmax,fmin)![maxmin]
!=============================================================================
! Like maxmin_d, except using single precision reals
use pkind, only: sp
implicit none
real(sp),intent(in ):: a,b,c,d   ! Cubic coefficients
real(sp),intent(in ):: vj,vk     ! End values of 2nd interval
real(sp),intent(out):: fmax,fmin ! Estimates of max and min
!-----------------------------------------------------------------------------
real(sp),parameter:: u0=0._sp,u1=1._sp,u2=2._sp,u3=3._sp
real(sp)          :: a3i,b2,bq,cq,x1,x2,r
!=============================================================================
fmax=max(vj,vk); fmin=min(vj,vk)! End values are the starting default
if(a==u0)then ! The special degenerate case of a quadratic polynomial:
   b2=b*u2
   if(abs(c)<abs(b2))call cubic(a,b,c,d,-c/b2,fmax,fmin)
else
! Express the quadratic equation for the extrema in the convenient form:
! x**2 - 2*bq*x + cq = 0:
   a3i=u1/(a*u3)
   bq=-b*a3i
   cq= c*a3i
! The discriminant of this quadratic is r; if r>0 we can take its sqrt
! and find two real and distinct roots (otherwise, we don't proceed):
   r=bq*bq-cq
   if(r>u0)then
      r=sqrt(r)
      if(bq/=u0)then ! Roots are distinct
         if(bq<u0)r=-r ! Choose sign of r to be same as that of bq
         x2=bq+r   ! <- The farther root (from x=0) of the quadratic
         x1=cq/x2  ! <- The closer root (from x=0) of the quadratic
         if(abs(x1)<u1)then ! Extremum x1 occurs inside the target
            call cubic(a,b,c,d,x1, fmax,fmin)
!              .. so check whether even this farther one, x2,  does
!              also:
            if(abs(x2)<u1)call cubic(a,b,c,d,x2, fmax,fmin)
         endif
      endif
   endif
endif
end subroutine maxmin_s

!=============================================================================
subroutine cubic_d(a,b,c,d,x, fmax,fmin)![cubic]
!=============================================================================
! Evaluate the cubic, ff = a*x**3 +b*x**2 +c*x +d, at given x and update
! either given max or min estimates, fmax or fmin, to contender value,
! ff, 
! if appropriate.
use pkind, only: dp
implicit none
real(dp),intent(in   ):: a,b,c,d,x
real(dp),intent(inout):: fmin,fmax
real(dp):: ff
!=============================================================================
ff=x*(x*(x*a+b)+c)+d; fmax=max(fmax,ff); fmin=min(fmin,ff)
end subroutine cubic_d
!=============================================================================
subroutine cubic_s(a,b,c,d,x, fmax,fmin)![cubic]
!=============================================================================
! Like cubic_d, except for single precision reals.
use pkind, only: sp
implicit none
real(sp),intent(in   ):: a,b,c,d,x
real(sp),intent(inout):: fmin,fmax
real(sp):: ff
!=============================================================================
ff=x*(x*(x*a+b)+c)+d; fmax=max(fmax,ff); fmin=min(fmin,ff)
end subroutine cubic_s

end module pmaxmin

            


