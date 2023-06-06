

MODULE nrtype

INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)

INTEGER, PARAMETER :: SP = KIND(1.0)
INTEGER, PARAMETER :: DP = KIND(1.0D0)

INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))

INTEGER, PARAMETER :: LGT = KIND(.true.)

! Frequently used mathematical constants (with precision to spare):

! REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp

! Derived data types for sparse matrices, single and double precision (see use in Chapter B2):

TYPE sprs2_sp
INTEGER(I4B) :: n,len
REAL(SP), DIMENSION(:), POINTER :: val
INTEGER(I4B), DIMENSION(:), POINTER :: irow
INTEGER(I4B), DIMENSION(:), POINTER :: jcol
END TYPE sprs2_sp
TYPE sprs2_dp
INTEGER(I4B) :: n,len
REAL(DP), DIMENSION(:), POINTER :: val
INTEGER(I4B), DIMENSION(:), POINTER :: irow
INTEGER(I4B), DIMENSION(:), POINTER :: jcol
END TYPE sprs2_dp
END MODULE nrtype




MODULE nrutil

! TABLE OF CONTENTS OF THE NRUTIL MODULE:
! routines that move data:
! array copy, swap, reallocate
! routines returning a location as an integer value
! i.rstloc, imaxloc, iminloc
! routines for argument checking and error handling:
! assert, assert eq, nrerror
! routines relating to polynomials and recurrences
! arth, geop, cumsum, cumprod, poly, polyterm,
! zroots unity
! routines for outer operations on vectors
! outerand, outersum, outerdi., outerprod, outerdiv
! routines for scatter-with-combine
! scatter add, scatter max
! routines for skewop erations on matrices
! diagadd, diagmult, get diag, put diag,
! unit matrix, lower triangle, upper triangle
! miscellaneous routines
! vabs

USE nrtype

! Parameters for crossover from serial to parallel algorithms (these are used only within this nrutil module):

IMPLICIT NONE
INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8

! Each NPAR2 must be ? the
! corresponding NPAR. 

INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
INTEGER(I4B), PARAMETER :: NPAR_POLY=8
INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8

! Next, generic interfaces for routines with overloaded versions. Naming conventions for appended
! codes in the names of overloaded routines are as follows: r=real, d=double precision,
! i=integer, c=complex, z=double-precision complex, h=character, l=logical. Any of
! r,d,i,c,z,h,l may be followed by v=vector or m=matrix (v,m su.xes are used only when
! needed to resolve ambiguities).
! Routines that move data:

INTERFACE array_copy
MODULE PROCEDURE array_copy_r, array_copy_d, array_copy_i
END INTERFACE
INTERFACE swap
MODULE PROCEDURE swap_i,swap_r,swap_rr,swap_rv,swap_rrv,swap_c, &
swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
masked_swap_rs,masked_swap_rv,masked_swap_rm
END INTERFACE
INTERFACE reallocate
MODULE PROCEDURE reallocate_rv,reallocate_rm,&
reallocate_iv,reallocate_im,reallocate_hv
END INTERFACE

! Routines returning a location as an integer value (i.rstloc, iminloc are not currently overloaded
! and so do not have a generic interface here):

INTERFACE imaxloc
MODULE PROCEDURE imaxloc_r,imaxloc_rr,imaxloc_i
END INTERFACE

! Routines for argument checking and error handling (nrerror is not currently overloaded):

INTERFACE assert
MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
END INTERFACE
INTERFACE assert_eq
MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
END INTERFACE

! Routines relating to polynomials and recurrences (cumprod, zroots unity are not currently
! overloaded):

INTERFACE arth
MODULE PROCEDURE arth_r, arth_d, arth_i
END INTERFACE
INTERFACE geop
MODULE PROCEDURE geop_r, geop_d, geop_i, geop_c, geop_dv
END INTERFACE
INTERFACE cumsum
MODULE PROCEDURE cumsum_r,cumsum_i
END INTERFACE
INTERFACE poly
MODULE PROCEDURE poly_rr,poly_rrv,poly_dd,poly_ddv,&
poly_rc,poly_cc,poly_msk_rrv,poly_msk_ddv
END INTERFACE
INTERFACE poly_term
MODULE PROCEDURE poly_term_rr,poly_term_cc
END INTERFACE

! Routines for outer operations on vectors (outerand, outersum, outerdiv are not currently
! overloaded):

INTERFACE outerprod

MODULE PROCEDURE outerprod_r,outerprod_d
END INTERFACE
INTERFACE outerdiff
MODULE PROCEDURE outerdiff_r,outerdiff_d,outerdiff_i
END INTERFACE

! Routines for scatter-with-combine, scatter add, scatter max:

INTERFACE scatter_add
MODULE PROCEDURE scatter_add_r,scatter_add_d
END INTERFACE
INTERFACE scatter_max
MODULE PROCEDURE scatter_max_r,scatter_max_d
END INTERFACE

! Routines for skewop erations on matrices (unit matrix, lower triangle, upper triangle not
! currently overloaded):

INTERFACE diagadd
MODULE PROCEDURE diagadd_rv,diagadd_r
END INTERFACE
INTERFACE diagmult
MODULE PROCEDURE diagmult_rv,diagmult_r
END INTERFACE
INTERFACE get_diag
MODULE PROCEDURE get_diag_rv, get_diag_dv
END INTERFACE
INTERFACE put_diag
MODULE PROCEDURE put_diag_rv, put_diag_r
END INTERFACE

! Other routines (vabs is not currently overloaded):

CONTAINS

SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)

REAL(SP), DIMENSION(:), INTENT(IN) :: src
REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
n_copied=min(size(src),size(dest))
n_not_copied=size(src)-n_copied
dest(1:n_copied)=src(1:n_copied)
END SUBROUTINE array_copy_r
SUBROUTINE array_copy_d(src,dest,n_copied,n_not_copied)
REAL(DP), DIMENSION(:), INTENT(IN) :: src
REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
n_copied=min(size(src),size(dest))
n_not_copied=size(src)-n_copied
dest(1:n_copied)=src(1:n_copied)
END SUBROUTINE array_copy_d
SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
INTEGER(I4B), DIMENSION(:), INTENT(IN) :: src
INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: dest
INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
n_copied=min(size(src),size(dest))
n_not_copied=size(src)-n_copied
dest(1:n_copied)=src(1:n_copied)
END SUBROUTINE array_copy_i
SUBROUTINE swap_i(a,b)

INTEGER(I4B), INTENT(INOUT) :: a,b
INTEGER(I4B) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_i

SUBROUTINE swap_r(a,b)
REAL(SP), INTENT(INOUT) :: a,b
REAL(SP) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_r

SUBROUTINE swap_rr(a,b)					!!! See routine on minu pandud topeltt2psuse jaoks
REAL(DP), INTENT(INOUT) :: a,b
REAL(DP) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_rr

SUBROUTINE swap_rv(a,b)
REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
REAL(SP), DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_rv

SUBROUTINE swap_rrv(a,b)				!!! See routine on minu pandud topeltt2psuse jaoks
REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
REAL(DP), DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_rrv

SUBROUTINE swap_c(a,b)
COMPLEX(SPC), INTENT(INOUT) :: a,b
COMPLEX(SPC) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_c
SUBROUTINE swap_cv(a,b)
COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_cv
SUBROUTINE swap_cm(a,b)
COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_cm
SUBROUTINE swap_z(a,b)
COMPLEX(DPC), INTENT(INOUT) :: a,b
COMPLEX(DPC) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_z
SUBROUTINE swap_zv(a,b)
COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_zv
SUBROUTINE swap_zm(a,b)
COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_zm
SUBROUTINE masked_swap_rs(a,b,mask)
REAL(SP), INTENT(INOUT) :: a,b
LOGICAL(LGT), INTENT(IN) :: mask
REAL(SP) :: swp

if (mask) then
swp=a
a=b
b=swp
end if
END SUBROUTINE masked_swap_rs
SUBROUTINE masked_swap_rv(a,b,mask)
REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
REAL(SP), DIMENSION(size(a)) :: swp
where (mask)
swp=a
a=b
b=swp
end where
END SUBROUTINE masked_swap_rv
SUBROUTINE masked_swap_rm(a,b,mask)
REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
where (mask)
swp=a
a=b
b=swp
end where
END SUBROUTINE masked_swap_rm
FUNCTION reallocate_rv(p,n)

REAL(SP), DIMENSION(:), POINTER :: p, reallocate_rv
INTEGER(I4B), INTENT(IN) :: n
INTEGER(I4B) :: nold,ierr
allocate(reallocate_rv(n),stat=ierr)
if (ierr /= 0) call &
nrerror('reallocate_rv: problem in attempt to allocate memory')
if (.not. associated(p)) RETURN
nold=size(p)
reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
deallocate(p)
END FUNCTION reallocate_rv
FUNCTION reallocate_iv(p,n)
INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
INTEGER(I4B), INTENT(IN) :: n
INTEGER(I4B) :: nold,ierr
allocate(reallocate_iv(n),stat=ierr)
if (ierr /= 0) call &
nrerror('reallocate_iv: problem in attempt to allocate memory')
if (.not. associated(p)) RETURN
nold=size(p)
reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
deallocate(p)
END FUNCTION reallocate_iv
FUNCTION reallocate_hv(p,n)
CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
INTEGER(I4B), INTENT(IN) :: n
INTEGER(I4B) :: nold,ierr
allocate(reallocate_hv(n),stat=ierr)
if (ierr /= 0) call &
nrerror('reallocate_hv: problem in attempt to allocate memory')
if (.not. associated(p)) RETURN
nold=size(p)
reallocate_hv(1:min(nold,n))=p(1:min(nold,n))

deallocate(p)
END FUNCTION reallocate_hv
FUNCTION reallocate_rm(p,n,m)
REAL(SP), DIMENSION(:,:), POINTER :: p, reallocate_rm
INTEGER(I4B), INTENT(IN) :: n,m
INTEGER(I4B) :: nold,mold,ierr
allocate(reallocate_rm(n,m),stat=ierr)
if (ierr /= 0) call &
nrerror('reallocate_rm: problem in attempt to allocate memory')
if (.not. associated(p)) RETURN
nold=size(p,1)
mold=size(p,2)
reallocate_rm(1:min(nold,n),1:min(mold,m))=&
p(1:min(nold,n),1:min(mold,m))
deallocate(p)
END FUNCTION reallocate_rm
FUNCTION reallocate_im(p,n,m)
INTEGER(I4B), DIMENSION(:,:), POINTER :: p, reallocate_im
INTEGER(I4B), INTENT(IN) :: n,m
INTEGER(I4B) :: nold,mold,ierr
allocate(reallocate_im(n,m),stat=ierr)
if (ierr /= 0) call &
nrerror('reallocate_im: problem in attempt to allocate memory')
if (.not. associated(p)) RETURN
nold=size(p,1)
mold=size(p,2)
reallocate_im(1:min(nold,n),1:min(mold,m))=&
p(1:min(nold,n),1:min(mold,m))
deallocate(p)
END FUNCTION reallocate_im

FUNCTION ifirstloc(mask)

LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
INTEGER(I4B) :: ifirstloc
INTEGER(I4B), DIMENSION(1) :: loc
loc=maxloc(merge(1,0,mask))
ifirstloc=loc(1)
if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
END FUNCTION ifirstloc
FUNCTION imaxloc_r(arr)

REAL(SP), DIMENSION(:), INTENT(IN) :: arr
INTEGER(I4B) :: imaxloc_r
INTEGER(I4B), DIMENSION(1) :: imax
imax=maxloc(arr(:))
imaxloc_r=imax(1)
END FUNCTION imaxloc_r

FUNCTION imaxloc_rr(arr)

REAL(DP), DIMENSION(:), INTENT(IN) :: arr				!!! See subroutine minu pandud, topeltt2psuse jaoks
INTEGER(I4B) :: imaxloc_rr
INTEGER(I4B), DIMENSION(1) :: imax
imax=maxloc(arr(:))
imaxloc_rr=imax(1)
END FUNCTION imaxloc_rr

FUNCTION imaxloc_i(iarr)
INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
INTEGER(I4B), DIMENSION(1) :: imax
INTEGER(I4B) :: imaxloc_i
imax=maxloc(iarr(:))
imaxloc_i=imax(1)
END FUNCTION imaxloc_i
FUNCTION iminloc(arr)

REAL(SP), DIMENSION(:), INTENT(IN) :: arr
INTEGER(I4B), DIMENSION(1) :: imin
INTEGER(I4B) :: iminloc
imin=minloc(arr(:))

iminloc=imin(1)
END FUNCTION iminloc

SUBROUTINE assert1(n1,string)

CHARACTER(LEN=*), INTENT(IN) :: string
LOGICAL, INTENT(IN) :: n1
if (.not. n1) then
write (*,*) 'nrerror: an assertion failed with this tag:', &
string
STOP 'program terminated by assert1'
end if
END SUBROUTINE assert1
SUBROUTINE assert2(n1,n2,string)
CHARACTER(LEN=*), INTENT(IN) :: string
LOGICAL, INTENT(IN) :: n1,n2
if (.not. (n1 .and. n2)) then
write (*,*) 'nrerror: an assertion failed with this tag:', &
string
STOP 'program terminated by assert2'
end if
END SUBROUTINE assert2
SUBROUTINE assert3(n1,n2,n3,string)
CHARACTER(LEN=*), INTENT(IN) :: string
LOGICAL, INTENT(IN) :: n1,n2,n3
if (.not. (n1 .and. n2 .and. n3)) then
write (*,*) 'nrerror: an assertion failed with this tag:', &
string
STOP 'program terminated by assert3'
end if
END SUBROUTINE assert3
SUBROUTINE assert4(n1,n2,n3,n4,string)
CHARACTER(LEN=*), INTENT(IN) :: string
LOGICAL, INTENT(IN) :: n1,n2,n3,n4
if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
write (*,*) 'nrerror: an assertion failed with this tag:', &
string
STOP 'program terminated by assert4'
end if
END SUBROUTINE assert4
SUBROUTINE assert_v(n,string)
CHARACTER(LEN=*), INTENT(IN) :: string
LOGICAL, DIMENSION(:), INTENT(IN) :: n
if (.not. all(n)) then
write (*,*) 'nrerror: an assertion failed with this tag:', &
string
STOP 'program terminated by assert_v'
end if
END SUBROUTINE assert_v
FUNCTION assert_eq2(n1,n2,string)

CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: n1,n2
INTEGER :: assert_eq2
if (n1 == n2) then
assert_eq2=n1
else
write (*,*) 'nrerror: an assert_eq failed with this tag:', &
string
STOP 'program terminated by assert_eq2'
end if

END FUNCTION assert_eq2
FUNCTION assert_eq3(n1,n2,n3,string)
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: n1,n2,n3
INTEGER :: assert_eq3
if (n1 == n2 .and. n2 == n3) then
assert_eq3=n1
else
write (*,*) 'nrerror: an assert_eq failed with this tag:', &
string
STOP 'program terminated by assert_eq3'
end if
END FUNCTION assert_eq3
FUNCTION assert_eq4(n1,n2,n3,n4,string)
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: n1,n2,n3,n4
INTEGER :: assert_eq4
if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
assert_eq4=n1
else
write (*,*) 'nrerror: an assert_eq failed with this tag:', &
string
STOP 'program terminated by assert_eq4'
end if
END FUNCTION assert_eq4
FUNCTION assert_eqn(nn,string)
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, DIMENSION(:), INTENT(IN) :: nn
INTEGER :: assert_eqn
if (all(nn(2:) == nn(1))) then
assert_eqn=nn(1)
else
write (*,*) 'nrerror: an assert_eq failed with this tag:', &
string
STOP 'program terminated by assert_eqn'
end if
END FUNCTION assert_eqn

SUBROUTINE nrerror(string)
CHARACTER(LEN=*), INTENT(IN) :: string
write (*,*) 'nrerror: ',string
STOP 'program terminated by nrerror'
END SUBROUTINE nrerror

FUNCTION arth_r(first,increment,n)

REAL(SP), INTENT(IN) :: first,increment
INTEGER(I4B), INTENT(IN) :: n
REAL(SP), DIMENSION(n) :: arth_r
INTEGER(I4B) :: k,k2
REAL(SP) :: temp
if (n > 0) arth_r(1)=first
if (n <= NPAR_ARTH) then
do k=2,n
arth_r(k)=arth_r(k-1)+increment
end do
else
do k=2,NPAR2_ARTH
arth_r(k)=arth_r(k-1)+increment
end do
temp=increment*NPAR2_ARTH
k=NPAR2_ARTH

do
if (k >= n) exit
k2=k+k
arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
temp=temp+temp
k=k2
end do
end if
END FUNCTION arth_r
FUNCTION arth_d(first,increment,n)
REAL(DP), INTENT(IN) :: first,increment
INTEGER(I4B), INTENT(IN) :: n
REAL(DP), DIMENSION(n) :: arth_d
INTEGER(I4B) :: k,k2
REAL(DP) :: temp
if (n > 0) arth_d(1)=first
if (n <= NPAR_ARTH) then
do k=2,n
arth_d(k)=arth_d(k-1)+increment
end do
else
do k=2,NPAR2_ARTH
arth_d(k)=arth_d(k-1)+increment
end do
temp=increment*NPAR2_ARTH
k=NPAR2_ARTH
do
if (k >= n) exit
k2=k+k
arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
temp=temp+temp
k=k2
end do
end if
END FUNCTION arth_d
FUNCTION arth_i(first,increment,n)
INTEGER(I4B), INTENT(IN) :: first,increment,n
INTEGER(I4B), DIMENSION(n) :: arth_i
INTEGER(I4B) :: k,k2,temp
if (n > 0) arth_i(1)=first
if (n <= NPAR_ARTH) then
do k=2,n
arth_i(k)=arth_i(k-1)+increment
end do
else
do k=2,NPAR2_ARTH
arth_i(k)=arth_i(k-1)+increment
end do
temp=increment*NPAR2_ARTH
k=NPAR2_ARTH
do
if (k >= n) exit
k2=k+k
arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
temp=temp+temp
k=k2
end do
end if
END FUNCTION arth_i
FUNCTION geop_r(first,factor,n)

REAL(SP), INTENT(IN) :: first,factor

INTEGER(I4B), INTENT(IN) :: n
REAL(SP), DIMENSION(n) :: geop_r
INTEGER(I4B) :: k,k2
REAL(SP) :: temp
if (n > 0) geop_r(1)=first
if (n <= NPAR_GEOP) then
do k=2,n
geop_r(k)=geop_r(k-1)*factor
end do
else
do k=2,NPAR2_GEOP
geop_r(k)=geop_r(k-1)*factor
end do
temp=factor**NPAR2_GEOP
k=NPAR2_GEOP
do
if (k >= n) exit
k2=k+k
geop_r(k+1:min(k2,n))=temp*geop_r(1:min(k,n-k))
temp=temp*temp
k=k2
end do
end if
END FUNCTION geop_r
FUNCTION geop_d(first,factor,n)
REAL(DP), INTENT(IN) :: first,factor
INTEGER(I4B), INTENT(IN) :: n
REAL(DP), DIMENSION(n) :: geop_d
INTEGER(I4B) :: k,k2
REAL(DP) :: temp
if (n > 0) geop_d(1)=first
if (n <= NPAR_GEOP) then
do k=2,n
geop_d(k)=geop_d(k-1)*factor
end do
else
do k=2,NPAR2_GEOP
geop_d(k)=geop_d(k-1)*factor
end do
temp=factor**NPAR2_GEOP
k=NPAR2_GEOP
do
if (k >= n) exit
k2=k+k
geop_d(k+1:min(k2,n))=temp*geop_d(1:min(k,n-k))
temp=temp*temp
k=k2
end do
end if
END FUNCTION geop_d
FUNCTION geop_i(first,factor,n)
INTEGER(I4B), INTENT(IN) :: first,factor,n
INTEGER(I4B), DIMENSION(n) :: geop_i
INTEGER(I4B) :: k,k2,temp
if (n > 0) geop_i(1)=first
if (n <= NPAR_GEOP) then
do k=2,n
geop_i(k)=geop_i(k-1)*factor
end do
else
do k=2,NPAR2_GEOP
geop_i(k)=geop_i(k-1)*factor
end do

temp=factor**NPAR2_GEOP
k=NPAR2_GEOP
do
if (k >= n) exit
k2=k+k
geop_i(k+1:min(k2,n))=temp*geop_i(1:min(k,n-k))
temp=temp*temp
k=k2
end do
end if
END FUNCTION geop_i
FUNCTION geop_c(first,factor,n)
COMPLEX(SP), INTENT(IN) :: first,factor
INTEGER(I4B), INTENT(IN) :: n
COMPLEX(SP), DIMENSION(n) :: geop_c
INTEGER(I4B) :: k,k2
COMPLEX(SP) :: temp
if (n > 0) geop_c(1)=first
if (n <= NPAR_GEOP) then
do k=2,n
geop_c(k)=geop_c(k-1)*factor
end do
else
do k=2,NPAR2_GEOP
geop_c(k)=geop_c(k-1)*factor
end do
temp=factor**NPAR2_GEOP
k=NPAR2_GEOP
do
if (k >= n) exit
k2=k+k
geop_c(k+1:min(k2,n))=temp*geop_c(1:min(k,n-k))
temp=temp*temp
k=k2
end do
end if
END FUNCTION geop_c
FUNCTION geop_dv(first,factor,n)
REAL(DP), DIMENSION(:), INTENT(IN) :: first,factor
INTEGER(I4B), INTENT(IN) :: n
REAL(DP), DIMENSION(size(first),n) :: geop_dv
INTEGER(I4B) :: k,k2
REAL(DP), DIMENSION(size(first)) :: temp
if (n > 0) geop_dv(:,1)=first(:)
if (n <= NPAR_GEOP) then
do k=2,n
geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
end do
else
do k=2,NPAR2_GEOP
geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
end do
temp=factor**NPAR2_GEOP
k=NPAR2_GEOP
do
if (k >= n) exit
k2=k+k
geop_dv(:,k+1:min(k2,n))=geop_dv(:,1:min(k,n-k))*&
spread(temp,2,size(geop_dv(:,1:min(k,n-k)),2))
temp=temp*temp
k=k2
end do
end if

END FUNCTION geop_dv
RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT(ans)

REAL(SP), DIMENSION(:), INTENT(IN) :: arr
REAL(SP), OPTIONAL, INTENT(IN) :: seed
REAL(SP), DIMENSION(size(arr)) :: ans
INTEGER(I4B) :: n,j
REAL(SP) :: sd
n=size(arr)
if (n == 0_i4b) RETURN
sd=0.0_sp
if (present(seed)) sd=seed
ans(1)=arr(1)+sd
if (n < NPAR_CUMSUM) then
do j=2,n
ans(j)=ans(j-1)+arr(j)
end do
else
ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
end if
END FUNCTION cumsum_r
RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT(ans)
INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
INTEGER(I4B), OPTIONAL, INTENT(IN) :: seed
INTEGER(I4B), DIMENSION(size(arr)) :: ans
INTEGER(I4B) :: n,j,sd
n=size(arr)
if (n == 0_i4b) RETURN
sd=0_i4b
if (present(seed)) sd=seed
ans(1)=arr(1)+sd
if (n < NPAR_CUMSUM) then
do j=2,n
ans(j)=ans(j-1)+arr(j)
end do
else
ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
end if
END FUNCTION cumsum_i
RECURSIVE FUNCTION cumprod(arr,seed) RESULT(ans)

REAL(SP), DIMENSION(:), INTENT(IN) :: arr
REAL(SP), OPTIONAL, INTENT(IN) :: seed
REAL(SP), DIMENSION(size(arr)) :: ans
INTEGER(I4B) :: n,j
REAL(SP) :: sd
n=size(arr)
if (n == 0_i4b) RETURN
sd=1.0_sp
if (present(seed)) sd=seed
ans(1)=arr(1)*sd
if (n < NPAR_CUMPROD) then
do j=2,n
ans(j)=ans(j-1)*arr(j)
end do
else
ans(2:n:2)=cumprod(arr(2:n:2)*arr(1:n-1:2),sd)
ans(3:n:2)=ans(2:n-1:2)*arr(3:n:2)
end if
END FUNCTION cumprod

FUNCTION poly_rr(x,coeffs)

REAL(SP), INTENT(IN) :: x
REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
REAL(SP) :: poly_rr
REAL(SP) :: pow
REAL(SP), DIMENSION(:), ALLOCATABLE :: vec
INTEGER(I4B) :: i,n,nn
n=size(coeffs)
if (n <= 0) then
poly_rr=0.0_sp
else if (n < NPAR_POLY) then
poly_rr=coeffs(n)
do i=n-1,1,-1
poly_rr=x*poly_rr+coeffs(i)
end do
else
allocate(vec(n+1))
pow=x
vec(1:n)=coeffs
do
vec(n+1)=0.0_sp
nn=ishft(n+1,-1)
vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
if (nn == 1) exit
pow=pow*pow
n=nn
end do
poly_rr=vec(1)
deallocate(vec)
end if
END FUNCTION poly_rr
FUNCTION poly_dd(x,coeffs)
REAL(DP), INTENT(IN) :: x
REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs
REAL(DP) :: poly_dd
REAL(DP) :: pow
REAL(DP), DIMENSION(:), ALLOCATABLE :: vec
INTEGER(I4B) :: i,n,nn
n=size(coeffs)
if (n <= 0) then
poly_dd=0.0_dp
else if (n < NPAR_POLY) then
poly_dd=coeffs(n)
do i=n-1,1,-1
poly_dd=x*poly_dd+coeffs(i)
end do
else
allocate(vec(n+1))
pow=x
vec(1:n)=coeffs
do
vec(n+1)=0.0_dp
nn=ishft(n+1,-1)
vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
if (nn == 1) exit
pow=pow*pow
n=nn
end do
poly_dd=vec(1)
deallocate(vec)
end if
END FUNCTION poly_dd

FUNCTION poly_rc(x,coeffs)
COMPLEX(SPC), INTENT(IN) :: x
REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
COMPLEX(SPC) :: poly_rc
COMPLEX(SPC) :: pow
COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
INTEGER(I4B) :: i,n,nn
n=size(coeffs)
if (n <= 0) then
poly_rc=0.0_sp
else if (n < NPAR_POLY) then
poly_rc=coeffs(n)
do i=n-1,1,-1
poly_rc=x*poly_rc+coeffs(i)
end do
else
allocate(vec(n+1))
pow=x
vec(1:n)=coeffs
do
vec(n+1)=0.0_sp
nn=ishft(n+1,-1)
vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
if (nn == 1) exit
pow=pow*pow
n=nn
end do
poly_rc=vec(1)
deallocate(vec)
end if
END FUNCTION poly_rc
FUNCTION poly_cc(x,coeffs)
COMPLEX(SPC), INTENT(IN) :: x
COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: coeffs
COMPLEX(SPC) :: poly_cc
COMPLEX(SPC) :: pow
COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
INTEGER(I4B) :: i,n,nn
n=size(coeffs)
if (n <= 0) then
poly_cc=0.0_sp
else if (n < NPAR_POLY) then
poly_cc=coeffs(n)
do i=n-1,1,-1
poly_cc=x*poly_cc+coeffs(i)
end do
else
allocate(vec(n+1))
pow=x
vec(1:n)=coeffs
do
vec(n+1)=0.0_sp
nn=ishft(n+1,-1)
vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
if (nn == 1) exit
pow=pow*pow
n=nn
end do
poly_cc=vec(1)
deallocate(vec)
end if
END FUNCTION poly_cc
FUNCTION poly_rrv(x,coeffs)

REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
REAL(SP), DIMENSION(size(x)) :: poly_rrv
INTEGER(I4B) :: i,n,m
m=size(coeffs)
n=size(x)
if (m <= 0) then
poly_rrv=0.0_sp
else if (m < n .or. m < NPAR_POLY) then
poly_rrv=coeffs(m)
do i=m-1,1,-1
poly_rrv=x*poly_rrv+coeffs(i)
end do
else
do i=1,n
poly_rrv(i)=poly_rr(x(i),coeffs)
end do
end if
END FUNCTION poly_rrv
FUNCTION poly_ddv(x,coeffs)
REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
REAL(DP), DIMENSION(size(x)) :: poly_ddv
INTEGER(I4B) :: i,n,m
m=size(coeffs)
n=size(x)
if (m <= 0) then
poly_ddv=0.0_dp
else if (m < n .or. m < NPAR_POLY) then
poly_ddv=coeffs(m)
do i=m-1,1,-1
poly_ddv=x*poly_ddv+coeffs(i)
end do
else
do i=1,n
poly_ddv(i)=poly_dd(x(i),coeffs)
end do
end if
END FUNCTION poly_ddv
FUNCTION poly_msk_rrv(x,coeffs,mask)
REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
REAL(SP), DIMENSION(size(x)) :: poly_msk_rrv
poly_msk_rrv=unpack(poly_rrv(pack(x,mask),coeffs),mask,0.0_sp)
END FUNCTION poly_msk_rrv
FUNCTION poly_msk_ddv(x,coeffs,mask)
REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
REAL(DP), DIMENSION(size(x)) :: poly_msk_ddv
poly_msk_ddv=unpack(poly_ddv(pack(x,mask),coeffs),mask,0.0_dp)
END FUNCTION poly_msk_ddv
RECURSIVE FUNCTION poly_term_rr(a,b) RESULT(u)

REAL(SP), DIMENSION(:), INTENT(IN) :: a
REAL(SP), INTENT(IN) :: b
REAL(SP), DIMENSION(size(a)) :: u
INTEGER(I4B) :: n,j
n=size(a)
if (n <= 0) RETURN
u(1)=a(1)
if (n < NPAR_POLYTERM) then
do j=2,n
u(j)=a(j)+b*u(j-1)
end do

else
u(2:n:2)=poly_term_rr(a(2:n:2)+a(1:n-1:2)*b,b*b)
u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
end if
END FUNCTION poly_term_rr
RECURSIVE FUNCTION poly_term_cc(a,b) RESULT(u)
COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
COMPLEX(SPC), INTENT(IN) :: b
COMPLEX(SPC), DIMENSION(size(a)) :: u
INTEGER(I4B) :: n,j
n=size(a)
if (n <= 0) RETURN
u(1)=a(1)
if (n < NPAR_POLYTERM) then
do j=2,n
u(j)=a(j)+b*u(j-1)
end do
else
u(2:n:2)=poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b)
u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
end if
END FUNCTION poly_term_cc
FUNCTION zroots_unity(n,nn)

INTEGER(I4B), INTENT(IN) :: n,nn
COMPLEX(SPC), DIMENSION(nn) :: zroots_unity
INTEGER(I4B) :: k
REAL(SP) :: theta
zroots_unity(1)=1.0
theta=TWOPI/n
k=1
do
if (k >= nn) exit
zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),SPC)
zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
zroots_unity(2:min(k,nn-k))
k=2*k
end do
END FUNCTION zroots_unity

FUNCTION outerprod_r(a,b)
REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
REAL(SP), DIMENSION(size(a),size(b)) :: outerprod_r
outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
spread(b,dim=1,ncopies=size(a))
END FUNCTION outerprod_r
FUNCTION outerprod_d(a,b)
REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
REAL(DP), DIMENSION(size(a),size(b)) :: outerprod_d
outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
spread(b,dim=1,ncopies=size(a))
END FUNCTION outerprod_d
FUNCTION outerdiv(a,b)
REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
REAL(SP), DIMENSION(size(a),size(b)) :: outerdiv
outerdiv = spread(a,dim=2,ncopies=size(b)) / &
spread(b,dim=1,ncopies=size(a))
END FUNCTION outerdiv
FUNCTION outersum(a,b)
REAL(SP), DIMENSION(:), INTENT(IN) :: a,b

REAL(SP), DIMENSION(size(a),size(b)) :: outersum
outersum = spread(a,dim=2,ncopies=size(b)) + &
spread(b,dim=1,ncopies=size(a))
END FUNCTION outersum
FUNCTION outerdiff_r(a,b)
REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
REAL(SP), DIMENSION(size(a),size(b)) :: outerdiff_r
outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
spread(b,dim=1,ncopies=size(a))
END FUNCTION outerdiff_r
FUNCTION outerdiff_d(a,b)
REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
REAL(DP), DIMENSION(size(a),size(b)) :: outerdiff_d
outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
spread(b,dim=1,ncopies=size(a))
END FUNCTION outerdiff_d
FUNCTION outerdiff_i(a,b)
INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
INTEGER(I4B), DIMENSION(size(a),size(b)) :: outerdiff_i
outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
spread(b,dim=1,ncopies=size(a))
END FUNCTION outerdiff_i
FUNCTION outerand(a,b)
LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: a,b
LOGICAL(LGT), DIMENSION(size(a),size(b)) :: outerand
outerand = spread(a,dim=2,ncopies=size(b)) .and. &
spread(b,dim=1,ncopies=size(a))
END FUNCTION outerand

SUBROUTINE scatter_add_r(dest,source,dest_index)
REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
REAL(SP), DIMENSION(:), INTENT(IN) :: source
INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
INTEGER(I4B) :: m,n,j,i
n=assert_eq2(size(source),size(dest_index),'scatter_add_r')
m=size(dest)
do j=1,n
i=dest_index(j)
if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
end do
END SUBROUTINE scatter_add_r
SUBROUTINE scatter_add_d(dest,source,dest_index)
REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
REAL(DP), DIMENSION(:), INTENT(IN) :: source
INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
INTEGER(I4B) :: m,n,j,i
n=assert_eq2(size(source),size(dest_index),'scatter_add_d')
m=size(dest)
do j=1,n
i=dest_index(j)
if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
end do
END SUBROUTINE scatter_add_d
SUBROUTINE scatter_max_r(dest,source,dest_index)
REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
REAL(SP), DIMENSION(:), INTENT(IN) :: source
INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
INTEGER(I4B) :: m,n,j,i
n=assert_eq2(size(source),size(dest_index),'scatter_max_r')
m=size(dest)
do j=1,n
i=dest_index(j)

if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
end do
END SUBROUTINE scatter_max_r
SUBROUTINE scatter_max_d(dest,source,dest_index)
REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
REAL(DP), DIMENSION(:), INTENT(IN) :: source
INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
INTEGER(I4B) :: m,n,j,i
n=assert_eq2(size(source),size(dest_index),'scatter_max_d')
m=size(dest)
do j=1,n
i=dest_index(j)
if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
end do
END SUBROUTINE scatter_max_d

SUBROUTINE diagadd_rv(mat,diag)

REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
REAL(SP), DIMENSION(:), INTENT(IN) :: diag
INTEGER(I4B) :: j,n
n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rv')
do j=1,n
mat(j,j)=mat(j,j)+diag(j)
end do
END SUBROUTINE diagadd_rv
SUBROUTINE diagadd_r(mat,diag)
REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
REAL(SP), INTENT(IN) :: diag
INTEGER(I4B) :: j,n
n = min(size(mat,1),size(mat,2))
do j=1,n
mat(j,j)=mat(j,j)+diag
end do
END SUBROUTINE diagadd_r
SUBROUTINE diagmult_rv(mat,diag)

REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
REAL(SP), DIMENSION(:), INTENT(IN) :: diag
INTEGER(I4B) :: j,n
n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
do j=1,n
mat(j,j)=mat(j,j)*diag(j)
end do
END SUBROUTINE diagmult_rv
SUBROUTINE diagmult_r(mat,diag)
REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
REAL(SP), INTENT(IN) :: diag
INTEGER(I4B) :: j,n
n = min(size(mat,1),size(mat,2))
do j=1,n
mat(j,j)=mat(j,j)*diag
end do
END SUBROUTINE diagmult_r
FUNCTION get_diag_rv(mat)

REAL(SP), DIMENSION(:,:), INTENT(IN) :: mat
REAL(SP), DIMENSION(size(mat,1)) :: get_diag_rv
INTEGER(I4B) :: j
j=assert_eq2(size(mat,1),size(mat,2),'get_diag_rv')
do j=1,size(mat,1)
get_diag_rv(j)=mat(j,j)

end do
END FUNCTION get_diag_rv
FUNCTION get_diag_dv(mat)
REAL(DP), DIMENSION(:,:), INTENT(IN) :: mat
REAL(DP), DIMENSION(size(mat,1)) :: get_diag_dv
INTEGER(I4B) :: j
j=assert_eq2(size(mat,1),size(mat,2),'get_diag_dv')
do j=1,size(mat,1)
get_diag_dv(j)=mat(j,j)
end do
END FUNCTION get_diag_dv
SUBROUTINE put_diag_rv(diagv,mat)

REAL(SP), DIMENSION(:), INTENT(IN) :: diagv
REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
INTEGER(I4B) :: j,n
n=assert_eq2(size(diagv),min(size(mat,1),size(mat,2)),'put_diag_rv')
do j=1,n
mat(j,j)=diagv(j)
end do
END SUBROUTINE put_diag_rv
SUBROUTINE put_diag_r(scal,mat)
REAL(SP), INTENT(IN) :: scal
REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
INTEGER(I4B) :: j,n
n = min(size(mat,1),size(mat,2))
do j=1,n
mat(j,j)=scal
end do
END SUBROUTINE put_diag_r
SUBROUTINE unit_matrix(mat)

REAL(SP), DIMENSION(:,:), INTENT(OUT) :: mat
INTEGER(I4B) :: i,n
n=min(size(mat,1),size(mat,2))
mat(:,:)=0.0_sp
do i=1,n
mat(i,i)=1.0_sp
end do
END SUBROUTINE unit_matrix
FUNCTION upper_triangle(j,k,extra)

INTEGER(I4B), INTENT(IN) :: j,k
INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
LOGICAL(LGT), DIMENSION(j,k) :: upper_triangle
INTEGER(I4B) :: n
n=0
if (present(extra)) n=extra
upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
END FUNCTION upper_triangle
FUNCTION lower_triangle(j,k,extra)

INTEGER(I4B), INTENT(IN) :: j,k
INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
LOGICAL(LGT), DIMENSION(j,k) :: lower_triangle
INTEGER(I4B) :: n
n=0
if (present(extra)) n=extra
lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
END FUNCTION lower_triangle

FUNCTION vabs(v)
REAL(SP), DIMENSION(:), INTENT(IN) :: v
REAL(SP) :: vabs
vabs=sqrt(dot_product(v,v))
END FUNCTION vabs
END MODULE nrutil




MODULE bessel

Contains

FUNCTION bessk_s(n,x)
USE nrtype; USE nrutil, ONLY : assert
! USE nr, ONLY : bessk0,bessk1 !!! ise
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: n
REAL(SP), INTENT(IN) :: x
REAL(SP) :: bessk_s
!!!   Returns the modi.ed Bessel function Kn(x) for positive x and n ? 2.
INTEGER(I4B) :: j
REAL(SP) :: bk,bkm,bkp,tox
! call assert(n, x, 'bessk_s args')
call assert(n >= 2, x > 0.0, 'bessk_s args')
tox=2.0_sp/x
bkm=bessk0_s(x) !!! Upward recurrence for all x...
bk=bessk1_s(x)
do j=1,n-1 !!! ...and here it is.
bkp=bkm+j*tox*bk
bkm=bk
bk=bkp
end do
bessk_s=bk
END FUNCTION bessk_s

FUNCTION bessk0_s(x)
USE nrtype; USE nrutil, ONLY : assert,poly
! USE nr, ONLY : bessi0 !!! ise
IMPLICIT NONE
REAL(SP), INTENT(IN) :: x
REAL(SP) :: bessk0_s
!!!     Returns the modi.ed Bessel function K0(x) for positive real x.
REAL(DP) :: y !!! Accumulate polynomials in double precision.
REAL(DP), DIMENSION(7) :: p = (/-0.57721566_dp,0.42278420_dp,&
0.23069756_dp,0.3488590e-1_dp,0.262698e-2_dp,0.10750e-3_dp,&
0.74e-5_dp/)
REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,-0.7832358e-1_dp,&
0.2189568e-1_dp,-0.1062446e-1_dp,0.587872e-2_dp,&
-0.251540e-2_dp,0.53208e-3_dp/)
call assert(x > 0.0, 'bessk0_s arg')
if (x <= 2.0) then      !!! Polynomial fit.
y=x*x/4.0_sp
bessk0_s=(-log(x/2.0_sp)*bessi0_s(x))+poly(y,p)
else
y=(2.0_sp/x)
bessk0_s=(exp(-x)/sqrt(x))*poly(y,q)
end if
END FUNCTION bessk0_s

FUNCTION bessk0_v(x)
USE nrtype; USE nrutil, ONLY : assert,poly
! USE nr, ONLY : bessi0 !!! ise
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(IN) :: x
REAL(SP), DIMENSION(size(x)) :: bessk0_v
REAL(DP), DIMENSION(size(x)) :: y
LOGICAL(LGT), DIMENSION(size(x)) :: mask
REAL(DP), DIMENSION(7) :: p = (/-0.57721566_dp,0.42278420_dp,&
0.23069756_dp,0.3488590e-1_dp,0.262698e-2_dp,0.10750e-3_dp,&
0.74e-5_dp/)
REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,-0.7832358e-1_dp,&
0.2189568e-1_dp,-0.1062446e-1_dp,0.587872e-2_dp,&
-0.251540e-2_dp,0.53208e-3_dp/)
call assert(all(x > 0.0), 'bessk0_v arg')
mask = (x <= 2.0)
where (mask)
y=x*x/4.0_sp
bessk0_v=(-log(x/2.0_sp)*bessi0_v(x))+poly(y,p,mask)
elsewhere
y=(2.0_sp/x)
bessk0_v=(exp(-x)/sqrt(x))*poly(y,q,.not. mask)
end where
END FUNCTION bessk0_v



FUNCTION bessk1_s(x)
USE nrtype; USE nrutil, ONLY : assert,poly
! USE nr, ONLY : bessi1 !!! ise
IMPLICIT NONE
REAL(SP), INTENT(IN) :: x
REAL(SP) :: bessk1_s
!!!     Returns the modi.ed Bessel function K1(x) for positive real x.
REAL(DP) :: y  !!!  Accumulate polynomials in double precision.
REAL(DP), DIMENSION(7) :: p = (/1.0_dp,0.15443144_dp,&
-0.67278579_dp,-0.18156897_dp,-0.1919402e-1_dp,&
-0.110404e-2_dp,-0.4686e-4_dp/)
REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,0.23498619_dp,&
-0.3655620e-1_dp,0.1504268e-1_dp,-0.780353e-2_dp,&
0.325614e-2_dp,-0.68245e-3_dp/)
call assert(x > 0.0, 'bessk1_s arg')
if (x <= 2.0) then  !!!  Polynomial fit.
y=x*x/4.0_sp
bessk1_s=(log(x/2.0_sp)*bessi1_s(x))+(1.0_sp/x)*poly(y,p)
else
y=2.0_sp/x
bessk1_s=(exp(-x)/sqrt(x))*poly(y,q)
end if
END FUNCTION bessk1_s


FUNCTION bessk1_v(x)
USE nrtype; USE nrutil, ONLY : assert,poly
! USE nr, ONLY : bessi1     !!! ise
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(IN) :: x
REAL(SP), DIMENSION(size(x)) :: bessk1_v
REAL(DP), DIMENSION(size(x)) :: y
LOGICAL(LGT), DIMENSION(size(x)) :: mask
REAL(DP), DIMENSION(7) :: p = (/1.0_dp,0.15443144_dp,&
-0.67278579_dp,-0.18156897_dp,-0.1919402e-1_dp,&
-0.110404e-2_dp,-0.4686e-4_dp/)
REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,0.23498619_dp,&
-0.3655620e-1_dp,0.1504268e-1_dp,-0.780353e-2_dp,&
0.325614e-2_dp,-0.68245e-3_dp/)
call assert(all(x > 0.0), 'bessk1_v arg')
mask = (x <= 2.0)
where (mask)
y=x*x/4.0_sp
bessk1_v=(log(x/2.0_sp)*bessi1_v(x))+(1.0_sp/x)*poly(y,p,mask)
elsewhere
y=2.0_sp/x
bessk1_v=(exp(-x)/sqrt(x))*poly(y,q,.not. mask)
end where
END FUNCTION bessk1_v



FUNCTION bessi0_s(x)
USE nrtype; USE nrutil, ONLY : poly
IMPLICIT NONE
REAL(SP), INTENT(IN) :: x
REAL(SP) :: bessi0_s
!!! Returns the modi.ed Bessel function I0(x) for any real x.
REAL(SP) :: ax
REAL(DP), DIMENSION(7) :: p = (/1.0_dp,3.5156229_dp,&
3.0899424_dp,1.2067492_dp,0.2659732_dp,0.360768e-1_dp,&
0.45813e-2_dp/)         !!! Accumulate polynomials in double precision.
REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,0.1328592e-1_dp,&
0.225319e-2_dp,-0.157565e-2_dp,0.916281e-2_dp,&
-0.2057706e-1_dp,0.2635537e-1_dp,-0.1647633e-1_dp,&
0.392377e-2_dp/)
ax=abs(x)
if (ax < 3.75) then         !!! Polynomial fit.
bessi0_s=poly(real((x/3.75_sp)**2,dp),p)
else
bessi0_s=(exp(ax)/sqrt(ax))*poly(real(3.75_sp/ax,dp),q)
end if
END FUNCTION bessi0_s

FUNCTION bessi0_v(x)
USE nrtype; USE nrutil, ONLY : poly
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(IN) :: x
REAL(SP), DIMENSION(size(x)) :: bessi0_v
REAL(SP), DIMENSION(size(x)) :: ax
REAL(DP), DIMENSION(size(x)) :: y
LOGICAL(LGT), DIMENSION(size(x)) :: mask
REAL(DP), DIMENSION(7) :: p = (/1.0_dp,3.5156229_dp,&
3.0899424_dp,1.2067492_dp,0.2659732_dp,0.360768e-1_dp,&
0.45813e-2_dp/)
REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,0.1328592e-1_dp,&
0.225319e-2_dp,-0.157565e-2_dp,0.916281e-2_dp,&
-0.2057706e-1_dp,0.2635537e-1_dp,-0.1647633e-1_dp,&
0.392377e-2_dp/)
ax=abs(x)
mask = (ax < 3.75)
where (mask)
bessi0_v=poly(real((x/3.75_sp)**2,dp),p,mask)
elsewhere
y=3.75_sp/ax
bessi0_v=(exp(ax)/sqrt(ax))*poly(real(y,dp),q,.not. mask)
end where
END FUNCTION bessi0_v




FUNCTION bessi1_s(x)
USE nrtype; USE nrutil, ONLY : poly
IMPLICIT NONE
REAL(SP), INTENT(IN) :: x
REAL(SP) :: bessi1_s
!!! Returns the modi.ed Bessel function I1(x) for any real x.
REAL(SP) :: ax
REAL(DP), DIMENSION(7) :: p = (/0.5_dp,0.87890594_dp,&
0.51498869_dp,0.15084934_dp,0.2658733e-1_dp,&
0.301532e-2_dp,0.32411e-3_dp/)
!!! Accumulate polynomials in double precision.
REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,-0.3988024e-1_dp,&
-0.362018e-2_dp,0.163801e-2_dp,-0.1031555e-1_dp,&
0.2282967e-1_dp,-0.2895312e-1_dp,0.1787654e-1_dp,&
-0.420059e-2_dp/)
ax=abs(x)
if (ax < 3.75) then     !!! Polynomial .t.
bessi1_s=ax*poly(real((x/3.75_sp)**2,dp),p)
else
bessi1_s=(exp(ax)/sqrt(ax))*poly(real(3.75_sp/ax,dp),q)
end if
if (x < 0.0) bessi1_s=-bessi1_s
END FUNCTION bessi1_s
FUNCTION bessi1_v(x)
USE nrtype; USE nrutil, ONLY : poly
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(IN) :: x
REAL(SP), DIMENSION(size(x)) :: bessi1_v
REAL(SP), DIMENSION(size(x)) :: ax
REAL(DP), DIMENSION(size(x)) :: y
LOGICAL(LGT), DIMENSION(size(x)) :: mask
REAL(DP), DIMENSION(7) :: p = (/0.5_dp,0.87890594_dp,&
0.51498869_dp,0.15084934_dp,0.2658733e-1_dp,&
0.301532e-2_dp,0.32411e-3_dp/)
REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,-0.3988024e-1_dp,&
-0.362018e-2_dp,0.163801e-2_dp,-0.1031555e-1_dp,&
0.2282967e-1_dp,-0.2895312e-1_dp,0.1787654e-1_dp,&
-0.420059e-2_dp/)
ax=abs(x)
mask = (ax < 3.75)
where (mask)
bessi1_v=ax*poly(real((x/3.75_sp)**2,dp),p,mask)
elsewhere
y=3.75_sp/ax
bessi1_v=(exp(ax)/sqrt(ax))*poly(real(y,dp),q,.not. mask)
end where
where (x < 0.0) bessi1_v=-bessi1_v
END FUNCTION bessi1_v



end module bessel




Module Besselfrac

use Realkind
use Constants

Contains


SUBROUTINE beschb_s(x,gam1,gam2,gampl,gammi)
USE nrtype
!!! USE nr, ONLY : chebev
IMPLICIT NONE
REAL(kind=rk), INTENT(IN) :: x
REAL(kind=rk), INTENT(OUT) :: gam1,gam2,gampl,gammi
INTEGER(I4B), PARAMETER :: NUSE1=5,NUSE2=5
!!! Evaluates Gamma1 and Gamma2 by Chebyshev expansion for |x| .le. 1/2. Also returns 1/gamma(1 + x) and
!!! 1/gamma(1-x). If converting to double precision, set NUSE1 = 7, NUSE2 = 8.
REAL(SP) :: xx
REAL(SP), DIMENSION(7) :: c1=(/-1.142022680371168_sp,&
    6.5165112670737e-3_sp,3.087090173086e-4_sp,-3.4706269649e-6_sp,&
    6.9437664e-9_sp,3.67795e-11_sp,-1.356e-13_sp/)
REAL(SP), DIMENSION(8) :: c2=(/1.843740587300905_sp,&
    -7.68528408447867e-2_sp,1.2719271366546e-3_sp,&
    -4.9717367042e-6_sp, -3.31261198e-8_sp,2.423096e-10_sp,&
    -1.702e-13_sp,-1.49e-15_sp/)
xx=8.0_dp*x*x-1.0_dp        !!! Multiply x by 2 to make range be 1 to 1, and then apply
                    !!! transformation for evaluating even Chebyshev
                    !!! series.
gam1=chebev_s(-1.0_sp,1.0_sp,c1(1:NUSE1),xx)
gam2=chebev_s(-1.0_sp,1.0_sp,c2(1:NUSE2),xx)
gampl=gam2-x*gam1
gammi=gam2+x*gam1
END SUBROUTINE beschb_s



SUBROUTINE beschb_v(x,gam1,gam2,gampl,gammi)
USE nrtype
!!! USE nr, ONLY : chebev
IMPLICIT NONE
REAL(DP), DIMENSION(:), INTENT(IN) :: x
REAL(DP), DIMENSION(:), INTENT(OUT) :: gam1,gam2,gampl,gammi
INTEGER(I4B), PARAMETER :: NUSE1=5,NUSE2=5
REAL(SP), DIMENSION(size(x)) :: xx
REAL(SP), DIMENSION(7) :: c1=(/-1.142022680371168_sp,&
    6.5165112670737e-3_sp,3.087090173086e-4_sp,-3.4706269649e-6_sp,&
    6.9437664e-9_sp,3.67795e-11_sp,-1.356e-13_sp/)
REAL(SP), DIMENSION(8) :: c2=(/1.843740587300905_sp,&
    -7.68528408447867e-2_sp,1.2719271366546e-3_sp,&
    -4.9717367042e-6_sp, -3.31261198e-8_sp,2.423096e-10_sp,&
    -1.702e-13_sp,-1.49e-15_sp/)
xx=8.0_dp*x*x-1.0_dp
gam1=chebev_v(-1.0_sp,1.0_sp,c1(1:NUSE1),xx)
gam2=chebev_v(-1.0_sp,1.0_sp,c2(1:NUSE2),xx)
gampl=gam2-x*gam1
gammi=gam2+x*gam1
END SUBROUTINE beschb_v



SUBROUTINE bessik(x,xnu,ri,rkk,rip,rkp)		! rk changed to rkk (19.03.06)
USE nrtype; USE nrutil, ONLY : assert,nrerror
! USE nr, ONLY : beschb
IMPLICIT NONE
REAL(kind=rk), INTENT(IN) :: x,xnu
REAL(kind=rk), INTENT(OUT) :: ri,rkk,rip,rkp	! rk changed to rkk (19.03.06)
INTEGER(I4B), PARAMETER :: MAXIT=10000 ! Enne oli 10000
REAL(kind=rk), PARAMETER :: XMIN=2.0
REAL(kind=rk), PARAMETER :: EPS=1.0e-10_dp,FPMIN=1.0e-30_dp
!!! Returns the modi.ed Bessel functions ri = Inu, rkk = Knu and their derivatives rip = I'nu, 
!!! rkp = K'nu, for positive x and for xnu = nu .ge. 0. The relative accuracy is within one or
!!! two significant digits of EPS. FPMIN is a number close to the machines smallest floatingpoint
!!! number. All internal arithmetic is in double precision. To convert the entire routine
!!! to double precision, change the REAL declaration above and decrease EPS to 10-16. Also
!!! convert the subroutine beschb.
INTEGER(I4B) :: i,l,nl
REAL(kind=rk) :: a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,&
    gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,&
    ril,ril1,rimu,rip1,ripl,ritemp,rk1,rkmu,rkmup,rktemp,&
    s,sum,sum1,x2,xi,xi2,xmu,xmu2
call assert(x > 0.0, xnu >= 0.0, 'bessik args')
nl=int(xnu+0.5_dp)          !!! nl is the number of downward recurrences
                        !!! of the Is and upward recurrences
                        !!! of Ks. xmu lies between 1/2 and 1/2.
xmu=xnu-nl
xmu2=xmu*xmu
xi=1.0_dp/x
xi2=2.0_dp*xi
h=xnu*xi                    !!! Evaluate CF1 by modi.ed Lentzs method
                        !!! (5.2). 
if (h < FPMIN) h=FPMIN
b=xi2*xnu
d=0.0
c=h
do i=1,MAXIT
    b=b+xi2
    d=1.0_dp/(b+d)              !!! Denominators cannot be zero here, so no
                        !!! need for special precautions. 
    c=b+1.0_dp/c
    del=c*d
    h=del*h
    if (abs(del-1.0_dp) < EPS) exit
end do
if (i > MAXIT) then
  if (x.gt.1e4_rk) then							!!! Added by Indrek Vurm
    rkk = sqrt(pi/(2.0_rk*x))*exp(-1.0_rk*x)				!!! Added by Indrek Vurm
    ri = 1.0_rk/sqrt(2.0_rk*x*pi)*exp(x)				!!! Added by Indrek Vurm
    rkp = -1.0_rk*sqrt(pi/2.0_rk)*exp(-1.0_rk*x)*((1.0_rk + 2.0_rk*x)/(2.0_rk*x**(3.0_rk/2.0_rk)))	!!! Added by Indrek Vurm
    rip = 1.0_rk/sqrt(2.0_rk*pi)*exp(x)*((2.0_rk*x - 1.0_rk)/(2.0_rk*x**(3.0_rk/2.0_rk)))		!!! Added by Indrek Vurm
    return								!!! Added by Indrek Vurm
  else									!!! Added by Indrek Vurm
    print *,'x=',x							!!! Added by Indrek Vurm
    call nrerror('x tooo large in bessik; try asymptotic expansion')
  end if								!!! Added by Indrek Vurm
end if
ril=FPMIN                   !!! Initialize Inu and I'nu for downward recurrence.
ripl=h*ril
ril1=ril                    !!! Store values for later rescaling.
rip1=ripl
fact=xnu*xi
do l=nl,1,-1
    ritemp=fact*ril+ripl
    fact=fact-xi
    ripl=fact*ritemp+ril
    ril=ritemp
end do
f=ripl/ril                  !!! Now have unnormalized Imu and I'mu.
if (x < XMIN) then          !!! Use series.
x2=0.5_dp*x
pimu=PI_D*xmu
if (abs(pimu) < EPS) then
fact=1.0
else
fact=pimu/sin(pimu)
end if
d=-log(x2)
e=xmu*d
if (abs(e) < EPS) then
fact2=1.0
else
fact2=sinh(e)/e
end if
call beschb_s(xmu,gam1,gam2,gampl,gammi)  !!! Chebyshev evaluation of ?1 and ?2.
ff=fact*(gam1*cosh(e)+gam2*fact2*d)     !!! f0.
sum=ff
e=exp(e)
p=0.5_dp*e/gampl                    !!! p0.
q=0.5_dp/(e*gammi)              !!! q0.
c=1.0
d=x2*x2

sum1=p
do i=1,MAXIT
    ff=(i*ff+p+q)/(i*i-xmu2)
    c=c*d/i
    p=p/(i-xmu)
    q=q/(i+xmu)
    del=c*ff
    sum=sum+del
    del1=c*(p-i*ff)
    sum1=sum1+del1
    if (abs(del) < abs(sum)*EPS) exit
end do
if (i > MAXIT) call nrerror('bessk series failed to converge')
rkmu=sum
rk1=sum1*xi2
else                            !!! Evaluate CF2 by Steeds algorithm (5.2),
                            !!! which is OK because there can be no
                            !!! zero denominators.
b=2.0_dp*(1.0_dp+x)
d=1.0_dp/b
delh=d
h=delh
q1=0.0                      !!! Initializations for recurrence (6.7.35).
q2=1.0
a1=0.25_dp-xmu2
c=a1
q=c                             !!! First term in equation (6.7.34).
a=-a1
s=1.0_dp+q*delh
do i=2,MAXIT
    a=a-2*(i-1)
    c=-a*c/i
    qnew=(q1-b*q2)/a
    q1=q2
    q2=qnew
    q=q+c*qnew
    b=b+2.0_dp
    d=1.0_dp/(b+a*d)
    delh=(b*d-1.0_dp)*delh
    h=h+delh
    dels=q*delh
    s=s+dels
    if (abs(dels/s) < EPS) exit             !!! Need only test convergence of sum, since
                                            !!! CF2 itself converges more quickly. 
end do
if (i > MAXIT) call nrerror('bessik: failure to converge in cf2')
h=a1*h
rkmu=sqrt(PI_D/(2.0_dp*x))*exp(-x)/s        !!! Omit the factor exp(x) to scale all the
                                !!! returned functions by exp(x) for x .ge. XMIN.
rk1=rkmu*(xmu+x+0.5_dp-h)*xi
end if
rkmup=xmu*xi*rkmu-rk1
rimu=xi/(f*rkmu-rkmup)                  !!! Get Imu from Wronskian.
ri=(rimu*ril1)/ril                  !!! Scale original Inu and I'nu.
rip=(rimu*rip1)/ril
do i=1,nl                           !!! Upward recurrence of Knu.
    rktemp=(xmu+i)*xi2*rk1+rkmu
    rkmu=rk1
    rk1=rktemp
end do
rkk=rkmu	! I changed it from rk to rkk (19.03.06)
rkp=xnu*xi*rkmu-rk1

END SUBROUTINE bessik


FUNCTION chebev_s(a,b,c,x)
USE nrtype; USE nrutil, ONLY : nrerror
IMPLICIT NONE
REAL(SP), INTENT(IN) :: a,b,x
REAL(SP), DIMENSION(:), INTENT(IN) :: c
REAL(SP) :: chebev_s
!!! Chebyshev evaluation: All arguments are input. c is an array of length M of Chebyshev
!!! coefficients, the first M elements of c output from chebft (which must have been called
!!! with the same a and b). The Chebyshev polynomial M
!!! k=1 ckTk.1(y)c1/2 is evaluated
!!! at a point y = [x(b+a)/2]/[(b a)/2], and the result is returned as the function value.
INTEGER(I4B) :: j,m
REAL(SP) :: d,dd,sv,y,y2
if ((x-a)*(x-b) > 0.0) call nrerror('x not in range in chebev_s')
m=size(c)
d=0.0
dd=0.0
y=(2.0_sp*x-a-b)/(b-a)                  !!! Change of variable.
y2=2.0_sp*y
do j=m,2,-1                         !!! Clenshaws recurrence.
sv=d
d=y2*d-dd+c(j)
dd=sv
end do
chebev_s=y*d-dd+0.5_sp*c(1)                 !!! Last step is different.
END FUNCTION chebev_s


FUNCTION chebev_v(a,b,c,x)
USE nrtype; USE nrutil, ONLY : nrerror
IMPLICIT NONE
REAL(SP), INTENT(IN) :: a,b
REAL(SP), DIMENSION(:), INTENT(IN) :: c,x
REAL(SP), DIMENSION(size(x)) :: chebev_v
INTEGER(I4B) :: j,m
REAL(SP), DIMENSION(size(x)) :: d,dd,sv,y,y2
if (any((x-a)*(x-b) > 0.0)) call nrerror('x not in range in chebev_v')
m=size(c)
d=0.0
dd=0.0
y=(2.0_sp*x-a-b)/(b-a)
y2=2.0_sp*y
do j=m,2,-1
sv=d
d=y2*d-dd+c(j)
dd=sv
end do
chebev_v=y*d-dd+0.5_sp*c(1)
END FUNCTION chebev_v



End Module Besselfrac





Module Expintegral

implicit none

Contains



FUNCTION expint(n,x)
USE nrtype; USE nrutil, ONLY : arth,assert,nrerror
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: n
REAL(DP), INTENT(IN) :: x
REAL(DP) :: expint
INTEGER(I4B), PARAMETER :: MAXIT=100
REAL(DP), PARAMETER :: EPS=epsilon(x),BIG=huge(x)*EPS
!   Evaluates the exponential integral En (x).
!   Parameters: MAXIT is the maximum allowed number of iterations; EPS is the desired relative
!   error, not smaller than the machine precision; BIG is a number near the largest representable
!   ﬂoating-point number; EULER (in nrtype) is Euler’s constant γ .
INTEGER(I4B) :: i,nm1
REAL(DP) :: a,b,c,d,del,fact,h
call assert(n >= 0, x >= 0.0, (x > 0.0 .or. n > 1), &
    'expint args')
if (n == 0) then			!	Special case.
    expint=exp(-x)/x
    RETURN
end if
nm1=n-1
if (x == 0.0) then			!	Another special case.
    expint=1.0_dp/nm1
else if (x > 1.0) then			!	Lentz’s algorithm (§5.2).
    b=x+n
    c=BIG
    d=1.0_dp/b
    h=d
    do i=1,MAXIT
         a=-i*(nm1+i)
         b=b+2.0_dp
         d=1.0_dp/(a*d+b)		!	Denominators cannot be zero.
         c=b+a/c
         del=c*d
         h=h*del
         if (abs(del-1.0_dp) <= EPS) exit
    end do
    if (i > MAXIT) call nrerror('expint: continued fraction failed')
    expint=h*exp(-x)
else				!	Evaluate series.
    if (nm1 /= 0) then		!	Set ﬁrst term.
         expint=1.0_dp/nm1
    else
         expint=-log(x)-EULER
    end if
    fact=1.0
    do i=1,MAXIT
         fact=-fact*x/i
         if (i /= nm1) then
           del=-fact/(i-nm1)

         else		!	ψ(n) appears here.
             del=fact*(-log(x)-EULER+sum(1.0_dp/arth(1,1,nm1)))
         end if
         expint=expint+del
         if (abs(del) < abs(expint)*EPS) exit
    end do
    if (i > MAXIT) call nrerror('expint: series failed')
end if
END FUNCTION expint


End Module Expintegral





Module dilogarithm

implicit none
Contains



      Function dilog(x,eps)

	!!! Dilogarithm. Acceptable argument values in the range -infty < x < 1.
	!!! Rel.error of the numerically calculated derivative is about 1e-8 if e = 1d-12

	real(kind=8), intent(in) :: x, eps
	real(kind=8) :: xl2, arg, gfn, dilog

	if (x.ge.0d0.and.x.lt.1.0d0) then
          xl2 = (log(1.0d0 - x))**2.0d0
          arg = 0.5d0*x/(1d0-x)
          gfn = G(arg,Eps)
          dilog = -0.5d0*xl2 + gfn
	else if (x.lt.0d0) then
	  arg = -0.5d0*x
          gfn = G(arg,eps)
          dilog = -gfn
	else if (x.eq.1.0d0) then
	  dilog = (3.1415926535897932384626433832795028841971)**2.0d0/6.0d0
	else if (x.gt.1.0d0) then
	  print *,'Module:dilogarithm, function: dilog; x>1', x
	  stop
	end if

      End Function dilog





      FUNCTION G(X,Eps)                                           
      IMPLICIT REAL*8(A-H,O-Z)                                  
      DATA XT,XTR /0.4D0,0.625D0/                               
      DATA PI,DL/3.141592653589792384D0,0.6931471805599453094D0/
      DATA PIH,PITW /1.64493406684823D0,0.822467033424113D0/    
      IF (X.LT.XT) GO TO 40                                     
      IF (X.GT.XTR) GO TO 40     
      G=PITW                      
      IF (X.EQ.0.5D0) RETURN      
      Y=2D0*X                     
      Z=DLOG(Y)                   
      G=G+DL*Z                    
      IF (X.GT.0.5D0) G=G+0.5D0*Z*Z
      IF (X.GT.0.5D0) Y=1D0/Y      
      Y=1D0-Y                      
      B=1D0                        
      U=0D0                        
      DK=1D0                       
      V=0D0                        
      A=1D0                        
   30 A=0.5D0*A                    
      V=V+A/DK                     
      DK=DK+1D0                    
      U=U+B*V/DK                   
      B=B*Y                       
      IF (B.GT.Eps) GO TO 30        
      U=U*Y*Y                     
      IF (X.GT.0.5D0) G=G-U       
      IF (X.LT.0.5D0) G=G+U       
      RETURN                      
   40 IF (X.LT.XT) Y=2D0*X        
      IF (X.GT.XTR) Y=0.5D0/X     
      G=1D0                       
      A=-Y                        
      DN=1D0                      
      B=1D0                       
   50 DN=DN+1D0                   
      B=A*B                       
      G=G+B/(DN*DN)               
      IF (DABS(B).GT.Eps) GO TO 50  
      G=Y*G                       
      IF (X.LT.XT) RETURN         
      Z=DLOG(Y)                   
            G=PIH+0.5D0*Z*Z-G     
      RETURN                      
      END function g                         

End Module dilogarithm







 
