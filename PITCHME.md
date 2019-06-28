## fortran-utils
## コードリーディング
@tkoyama010

---


### fortran-utilsとは？

[certik](https://github.com/certik)氏により2016年まで開発されていたfortranの便利なutils集です。

https://github.com/certik/fortran-utils

一部テストが通らない部分があったのでforkして変更をしました。

https://github.com/tkoyama010/fortran-utils

---


### 今回は線形代数部分を詳しく！

* Types (``dp``)
* Constants (``pi``, ``e_``, ``i_``)
* Sorting
* Saving/loading 2D arrays (``savetxt``, ``loadtxt``)
* Meshes (exponential, uniform)
* Cubic splines
* Saving/loading PPM images
* [Lapack interface (and a few simple f90 wrappers like ``eigh``, ``inv``)](https://github.com/tkoyama010/fortran-utils/blob/master/src/linalg.f90)
* HDF5 interface

---

```fortran
program test_det
use types, only: dp
use utils, only: assert
use linalg, only: det, eye
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(3,3), Adet
complex(dp) :: B(3,3), Bdet

A = real(reshape([ 1, 2, 3, 4, 5, 6, 7, 8, -9 ], shape=[3,3]))
Adet = det(A)
call assert(abs(Adet - 54.0_dp) < eps)

B = 0*i_
B = cmplx(A)
B(1,1) = i_
! B has det(B) = 147 - 93i (SymPy)
Bdet = det(B)
call assert(abs(Bdet - (147.0_dp, -93.0_dp)) < eps)

end program
```
@[1-6](まずは行列式を計算する関数detについて解説します。)
@[8-10](実数と複素数の3×3の行列と行列式のための変数を定義します。)
@[12-14](実数型の計算と確認はこのように行います。asser文を使っていますが、これは論理型を引数にした関数として定義されています。)
@[16-21](実数の行列の[1,1]の部分を虚数iにした行列の行列式はこのように計算します。)


---

```fortran
  ! determinants of real/complex square matrices:
  interface det
     module procedure ddet
     module procedure zdet
  end interface det
```
@[1-5](module procedureを使用して実数の引数を持つ場合と複素数の引数を持つ場合の処理を分けています。)

---

```fortran
  function ddet(A) result(x)

    n = size(A(1,:))
    call assert_shape(A, [n, n], "det", "A")
    allocate(At(n,n), ipiv(n))
    At = A
    call dgetrf(n, n, At, n, ipiv, info)

    x = 1.0_dp
    do i = 1,n
       if(ipiv(i) /= i) then  ! additional sign change
          x = -x*At(i,i)
       else
          x = x*At(i,i)
       endif
    end do
  end function ddet
```
@[1](実数と複素数でソースの流れは同じですので、実数の方のみを説明します。)
@[3-7](行列式計算にはLapackのLU分解のルーチンdgetrfを使用します。)
@[9-16](LU分解された行列Atの対角項を使用して行列式を計算します。)


---

```fortran
program test_diag
use types, only: dp
use utils, only: assert
use linalg, only: diag, eye
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(3,3), B(3,3)
complex(dp) :: C(3,3), D(3,3)

A = diag([1.0_dp, 1.0_dp, 1.0_dp])
call assert(maxval(abs(A-eye(3))) < eps)

A = diag([1.0_dp, 2.0_dp, 3.0_dp])
B = reshape([1,0,0,0,2,0,0,0,3], shape=[3,3])
call assert(maxval(abs(A - B)) < eps)

C = diag([1*i_, 2*i_, 3*i_])
D = i_*reshape([1,0,0,0,2,0,0,0,3], shape=[3,3])
call assert(maxval(abs(C - D)) < eps)
end program
```
@[1-6](次にdiagについて説明します。)
@[8-23](実数の行列Aをdiagで定義して別途構築した行列と比較をしています。)

---

```fortran
  ! construction of square matrices from the diagonal elements:
  interface diag
     module procedure ddiag
     module procedure zdiag
  end interface diag

  function ddiag(x) result(A)
    ! construct real matrix from diagonal elements
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: A(:,:)
    integer :: i, n

    n = size(x)
    allocate(A(n,n))
    A(:,:) = 0.0_dp
    forall(i=1:n) A(i,i) = x(i)
  end function ddiag
```
@[1-5](module procedureを使用するのは同様です。)
@[7-11](同様に実数の場合について説明します。)
@[13-17](forall文を使い対角項に値を代入してきます。)

---

```fortran
program test_eig
use types, only: dp
use utils, only: assert
use linalg, only: eig, eye
use constants, only : i_
implicit none

! test eigenvalue comutation for general matrices:

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(2, 2), B(5, 5)
complex(dp) :: AC(2, 2), lam(2), c(2, 2), r(2), n, lamb(5), cb(5, 5)
integer :: i

! test a matrix with complex eigenvalues/eigenvectors:
A = reshape([1, -1, 2, 1], shape=[2, 2])  ! lambda_i = 1 \pm \sqrt(2) i
call eig(A, lam, c)
! TODO: add test for correctness of eigenvalues
do i = 1, 2
    ! Test proper norm of eigenvectors:
    n = dot_product(c(:,i), c(:,i))
    call assert(abs(n - 1) < eps)
    ! Test that c(:, i) is an eigenvector with eigenvalue lam(i):
    r = matmul(A-lam(i)*cmplx(eye(2)), c(:, i))
    call assert(sqrt(abs(dot_product(r, r))) < eps)
end do

! test a multiple of the unit matrix:
B = 3*eye(5)
call eig(B, lamb, cb)
call assert(maxval(abs(lamb - 3.0_dp)) < eps)  ! all eigenvalues are 3
call assert(maxval(abs(cb - cmplx(eye(5)))) < eps)  ! eigenvectors are cartesian unit basis vectors

! test complex matrices:
AC = reshape([1.0_dp+0*i_, 2*i_, 3*i_, -4.0_dp+0*i_], shape=[2,2])
call eig(AC, lam, c)
call assert(all(abs(lam - cmplx([-1.0_dp, -2.0_dp])) < eps))
do i = 1, 2
    ! Test proper norm of eigenvectors:
    n = dot_product(c(:,i), c(:,i))
    call assert(abs(n - 1) < eps)
    ! Test that c(:, i) is an eigenvector with eigenvalue lam(i):
    r = matmul(AC-lam(i)*cmplx(eye(2)), c(:, i))
    call assert(sqrt(abs(dot_product(r, r))) < eps)
end do

end program
```

---

```fortran
  ! eigenvalue/-vector problem for general matrices:
  interface eig
     module procedure deig
     module procedure zeig
  end interface eig

  ! TODO: add optional switch for left or right eigenvectors in deig() and zeig()?
  subroutine deig(A, lam, c)
    real(dp), intent(in) :: A(:, :)  ! matrix for eigenvalue compuation
    complex(dp), intent(out) :: lam(:)  ! eigenvalues: A c = lam c
    complex(dp), intent(out) :: c(:, :)  ! eigenvectors: A c = lam c; c(i,j) = ith component of jth vec.
    ! LAPACK variables for DGEEV:
    real(dp), allocatable ::  At(:,:), vl(:,: ), vr(:,:), wi(:), work(:), wr(:)
    integer :: info, lda, ldvl, ldvr, lwork, n, i

    lda = size(A(:,1))
    n = size(A(1,:))
    call assert_shape(A, [n, n], "solve", "A")
    call assert_shape(c, [n, n], "solve", "c")
    ldvl = n
    ldvr = n
    lwork = 8*n  ! TODO: can this size be optimized? query first?
    allocate(At(lda,n), wr(n), wi(n), vl(ldvl,n), vr(ldvr,n), work(lwork))
    At = A

    call dgeev('N', 'V', n, At, lda, wr, wi, vl, ldvl, vr, ldvr, &
         work, lwork, info)
    if(info /= 0) then
       print *, "dgeev returned info = ", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "the QR algorithm failed to compute all the"
          print *, "eigenvalues, and no eigenvectors have been computed;"
          print *, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
          print *, "have converged."
       end if
       call stop_error('eig: dgeev error')
    end if

    lam = wr + i_*wi
    ! as DGEEV has a rather complicated way of returning the eigenvectors,
    ! it is necessary to build the complex array of eigenvectors from
    ! two real arrays:
    do i = 1,n
       if(wi(i) > 0.0) then  ! first of two conjugate eigenvalues
          c(:, i) = vr(:, i) + i_*vr(:, i+1)
       elseif(wi(i) < 0.0_dp) then  ! second of two conjugate eigenvalues
          c(:, i) = vr(:, i-1) - i_*vr(:, i)
       else
          c(:, i) = vr(:, i)
       end if
    end do
  end subroutine deig

  subroutine zeig(A, lam, c)
    complex(dp), intent(in) :: A(:, :)  ! matrix to solve eigenproblem for
    complex(dp), intent(out) :: lam(:)  ! eigenvalues: A c = lam c
    complex(dp), intent(out) :: c(:,:)  ! eigenvectors: A c = lam c; c(i,j) = ith component of jth vec.
    ! LAPACK variables:
    integer :: info, lda, ldvl, ldvr, lwork, n, lrwork
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: vl(:,:), vr(:,:), work(:)

    lda = size(A(:,1))
    n = size(A(1,:))
    call assert_shape(A, [n, n], "solve", "A")
    call assert_shape(c, [n, n], "solve", "c")
    ldvl = n
    ldvr = n
    lwork = 8*n  ! TODO: can this size be optimized? query first?
    lrwork = 2*n
    allocate(vl(ldvl,n), vr(ldvr,n), work(lwork), rwork(lrwork))
    c = A
    call zgeev('N', 'V', n, c, lda, lam, vl, ldvl, vr, ldvr, work, &
         lwork, rwork, info)
    if(info /= 0) then
       print *, "zgeev returned info = ", info
       if(info < 0) then
          print *, "the ",-info, "-th argument had an illegal value."
       else
          print *, "the QR algorithm failed to compute all the"
          print *, "eigenvalues, and no eigenvectors have been computed;"
          print *, "elements and ", info+1, ":", n, " of W contain eigenvalues which have"
          print *, "converged."
       end if
       call stop_error('eig: zgeev error')
    end if
    c = vr
  end subroutine zeig
```

---

```fortran
program test_eigh
use types, only: dp
use utils, only: assert
use linalg, only: eigh, eye
use constants, only: i_
implicit none

! test eigenvalue computation for symmetric/hermitian matrices:

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(2, 2), B(2, 2), lam(2), c(2, 2), r(2), n
complex(dp) :: AC(2, 2), BC(2, 2), cc(2, 2)
complex(dp) :: rc(2), nc
integer :: i

A = reshape([1, 0, 0, -1], [2, 2])
B = reshape([1, 0, 0, 1], [2, 2])

! test generalized eigenvalue problem with real symmetric matrices:
call eigh(A, B, lam, c)
! Test eigenvalues:
call assert(all(abs(lam - [-1, 1]) < eps))
do i = 1, 2
    ! Test that c(:, i) is an eigenvector with eigenvalue lam(i):
    r = matmul(A-lam(i)*B, c(:, i))
    call assert(sqrt(dot_product(r, r)) < eps)
    ! Test that eigenvectors are properly normalized:
    n = dot_product(c(:, i), matmul(B, c(:, i)))
    call assert(abs(n - 1) < eps)
end do

! test eigenvalue problem with real symmetric matrices:
call eigh(A, lam, c)
! Test eigenvalues:
call assert(all(abs(lam - [-1, 1]) < eps))
do i = 1, 2
    ! Test that c(:, i) is an eigenvector with eigenvalue lam(i):
    r = matmul(A, c(:, i)) - lam(i) * c(:, i)
    call assert(sqrt(dot_product(r, r)) < eps)
    ! Test that eigenvectors are properly normalized:
    n = dot_product(c(:, i), c(:, i))
    call assert(abs(n - 1) < eps)
end do

! another test for the generalized real problem:
A = reshape([2, -4, -4, 2], [2, 2])
B = reshape([2, 1, 1, 2], [2, 2])
call eigh(A, B, lam, c)
! Test eigenvalues:
call assert(all(abs(lam - [-2._dp/3, 6._dp]) < eps))
do i = 1, 2
    ! Test that c(:, i) are eigenvectors:
    r = matmul(A-lam(i)*B, c(:, i))
    call assert(sqrt(dot_product(r, r)) < eps)
    ! Test that eigenvectors are properly normalized:
    n = dot_product(c(:, i), matmul(B, c(:, i)))
    call assert(abs(n - 1) < eps)
end do

! test with complex matrices:
AC = reshape([1.0_dp+0.0_dp*i_, i_, -i_, 3.0_dp+0.0_dp*i_], [2, 2])
call eigh(AC, lam, cc)
do i = 1, 2
    ! Test that c(:, i) are eigenvectors:
    rc = matmul(AC, cc(:, i)) - lam(i) * cc(:, i)
    call assert(sqrt(abs(dot_product(rc, rc))) < eps)
    ! Test that eigenvectors are properly normalized:
    nc = dot_product(cc(:, i), cc(:, i))
    call assert(abs(nc - 1) < eps)
end do

AC = reshape([1.0_dp +0*i_, 3*i_, -3*i_, 2.0_dp-0*i_], shape=[2,2])
BC = reshape([1.0_dp + 0*i_, 1*i_, -1*i_, 2.0_dp+ 0*i_], shape=[2,2])
call eigh(AC, BC, lam, cc)
! Test eigenvalues:
! TODO: add test for correctness eigenvalues
do i = 1, 2
    ! Test that cc(:, i) are eigenvectors:
    rc = matmul(AC-lam(i)*BC, cc(:, i))
    call assert(sqrt(abs(dot_product(rc, rc))) < eps)
end do
! Test that eigenvectors are properly normalized:
BC = matmul(matmul(conjg(transpose(cc)), BC), cc)  ! should be eye(2); Z^H B Z = I
call assert(maxval(abs(BC - cmplx(eye(2)))) < eps)

end program
```

---

```fortran
  ! eigenvalue/-vector problem for real symmetric/complex hermitian matrices:
  interface eigh
     module procedure deigh_generalized
     module procedure deigh_simple
     module procedure zeigh_generalized
     module procedure zeigh_simple
  end interface eigh


  subroutine deigh_generalized(Am, Bm, lam, c)
    ! solves generalized eigen value problem for all eigenvalues and eigenvectors
    ! Am must by symmetric, Bm symmetric positive definite.
    ! Only the lower triangular part of Am and Bm is used.
    real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
    real(dp), intent(in) :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
    real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam Bm c
    real(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
    integer :: n
    ! lapack variables
    integer :: lwork, liwork, info
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: Bmt(:,:), work(:)

    ! solve
    n = size(Am,1)
    call assert_shape(Am, [n, n], "eigh", "Am")
    call assert_shape(Bm, [n, n], "eigh", "B")
    call assert_shape(c, [n, n], "eigh", "c")
    lwork = 1 + 6*n + 2*n**2
    liwork = 3 + 5*n
    allocate(Bmt(n,n), work(lwork), iwork(liwork))
    c = Am; Bmt = Bm  ! Bmt temporaries overwritten by dsygvd
    call dsygvd(1,'V','L',n,c,n,Bmt,n,lam,work,lwork,iwork,liwork,info)
    if (info /= 0) then
       print *, "dsygvd returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else if (info <= n) then
          print *, "the algorithm failed to compute an eigenvalue while working"
          print *, "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
          print *, "through", mod(info, n+1)
       else
          print *, "The leading minor of order ", info-n, &
               "of B is not positive definite. The factorization of B could ", &
               "not be completed and no eigenvalues or eigenvectors were computed."
       end if
       call stop_error('eigh: dsygvd error')
    end if
  end subroutine deigh_generalized

  subroutine deigh_simple(Am, lam, c)
    ! solves eigen value problem for all eigenvalues and eigenvectors
    ! Am must by symmetric
    ! Only the lower triangular part of Am is used.
    real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam c
    real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam c
    real(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam c; c(i,j) = ith component of jth vec.
    integer :: n
    ! lapack variables
    integer :: lwork, liwork, info
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: work(:)

    ! solve
    n = size(Am,1)
    call assert_shape(Am, [n, n], "eigh", "Am")
    call assert_shape(c, [n, n], "eigh", "c")
    lwork = 1 + 6*n + 2*n**2
    liwork = 3 + 5*n
    allocate(work(lwork), iwork(liwork))
    c = Am
    call dsyevd('V','L',n,c,n,lam,work,lwork,iwork,liwork,info)
    if (info /= 0) then
       print *, "dsyevd returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "the algorithm failed to compute an eigenvalue while working"
          print *, "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
          print *, "through", mod(info, n+1)
       end if
       call stop_error('eigh: dsyevd error')
    end if
  end subroutine deigh_simple

  subroutine zeigh_generalized(Am, Bm, lam, c)
    ! solves generalized eigen value problem for all eigenvalues and eigenvectors
    ! Am must by hermitian, Bm hermitian positive definite.
    ! Only the lower triangular part of Am and Bm is used.
    complex(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
    complex(dp), intent(in) :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
    real(dp), intent(out) :: lam(:)      ! eigenvalues: Am c = lam Bm c
    complex(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
    ! lapack variables
    integer :: info, liwork, lrwork, lwork, n
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: Bmt(:,:), work(:)

    ! solve
    n = size(Am,1)
    call assert_shape(Am, [n, n], "eigh", "Am")
    call assert_shape(Bm, [n, n], "eigh", "Bm")
    call assert_shape(c, [n, n], "eigh", "c")
    lwork = 2*n + n**2
    lrwork = 1 + 5*N + 2*n**2
    liwork = 3 + 5*n
    allocate(Bmt(n,n), work(lwork), rwork(lrwork), iwork(liwork))
    c = Am; Bmt = Bm  ! Bmt temporary overwritten by zhegvd
    call zhegvd(1,'V','L',n,c,n,Bmt,n,lam,work,lwork,rwork,lrwork,iwork,liwork,info)
    if (info /= 0) then
       print *, "zhegvd returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else if (info <= n) then
          print *, "the algorithm failed to compute an eigenvalue while working"
          print *, "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
          print *, "through", mod(info, n+1)
       else
          print *, "The leading minor of order ", info-n, &
               "of B is not positive definite. The factorization of B could ", &
               "not be completed and no eigenvalues or eigenvectors were computed."
       end if
       call stop_error('eigh: zhegvd error')
    end if
  end subroutine zeigh_generalized

  subroutine zeigh_simple(Am, lam, c)
    ! solves eigen value problem for all eigenvalues and eigenvectors
    ! Am must by symmetric
    ! Only the lower triangular part of Am is used.
    complex(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam c
    real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam c
    complex(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam c; c(i,j) = ith component of jth vec.
    ! LAPACK variables:
    integer :: info, lda, liwork, lrwork, lwork, n
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: work(:)

    ! use LAPACK's zheevd routine
    n = size(Am, 1)
    call assert_shape(Am, [n, n], "eigh", "Am")
    call assert_shape(c, [n, n], "eigh", "c")
    lda = max(1, n)
    lwork = 2*n + n**2
    lrwork = 1 + 5*n + 2*n**2
    liwork = 3 + 5*n
    allocate(work(lwork), rwork(lrwork), iwork(liwork))
    c = Am
    call zheevd("V", "L", n, c, lda, lam, work, lwork, rwork, lrwork, &
         iwork, liwork, info)
    if (info /= 0) then
       print *, "zheevd returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "the algorithm failed to compute an eigenvalue while working"
          print *, "through the submatrix lying in rows and columns through"
          print *, info/(n+1), " through ", mod(info, n+1)
       end if
       call stop_error('eigh: zheevd error')
    end if
  end subroutine zeigh_simple

```

---

```fortran
program test_eig_bigger
! You can adapt this test to test speed for bigger matrices.
use types, only: dp
use utils, only: assert
use linalg, only: eigh
implicit none

real(dp), parameter :: eps = 1e-9_dp
integer, parameter :: n = 10
real(dp) :: A(n, n), lam(n), c(n, n), r(n), norm
real(dp) :: t1, t2
integer :: co1, co2, comax, rate
integer :: i
A = 0
forall(i = 1:n) A(i, i) = 2
forall(i = 1:n-1) A(i, i+1) = -1
forall(i = 1:n-1) A(i+1, i) = -1
!To print the array (if n=10):
!print "(10(f6.2))", A
call cpu_time(t1)
call system_clock(co1, rate, comax)
call eigh(A, lam, c)
call cpu_time(t2)
call system_clock(co2, rate, comax)
print *, "Time (cpu):", t2-t1
print *, "Time (all):", (co2-co1) / real(rate, dp)
print *, "First 10 eigenvalues:"
print *, lam(:10)
do i = 1, n
    ! Test that c(:, i) are eigenvectors:
    r = matmul(A, c(:, i)) - lam(i) * c(:, i)
    call assert(sqrt(dot_product(r, r)) < eps)
    ! Test that eigenvectors are properly normalized:
    norm = dot_product(c(:, i), c(:, i))
    call assert(abs(norm - 1) < eps)
end do

end program
```

---

```fortran
program test_eigvals
use types, only: dp
use utils, only: assert
use linalg, only: eigvals, eye
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: B(5, 5)
complex(dp) :: AC(2, 2), lam(2), lamb(5)

! test a general matrix
AC = reshape([1.0_dp+0*i_, 2*i_, 3*i_, -4.0_dp+0*i_], shape=[2,2])
lam = eigvals(AC)
call assert(all(abs(lam - cmplx([-1.0_dp, -2.0_dp])) < eps))

! test a multiple of the unit matrix:
B = 3*eye(5)
lamb = eigvals(B)
call assert(maxval(abs(lamb - 3.0_dp)) < eps)  ! all eigenvalues are 3

! TODO: add more tests

end program
```

---

```fortran
  ! eigenvalues for general matrices:
  interface eigvals
     module procedure deigvals
     module procedure zeigvals
  end interface eigvals

  function deigvals(A) result(lam)
    real(dp), intent(in) :: A(:, :)  ! matrix for eigenvalue compuation
    complex(dp), allocatable :: lam(:)  ! eigenvalues: A c = lam c
    ! LAPACK variables for DGEEV:
    real(dp), allocatable ::  At(:,:), vl(:,: ), vr(:,:), wi(:), work(:), wr(:)
    integer :: info, lda, ldvl, ldvr, lwork, n

    lda = size(A(:,1))
    n = size(A(1,:))
    call assert_shape(A, [n, n], "solve", "A")
    ldvl = n
    ldvr = n
    lwork = 8*n  ! TODO: can this size be optimized? query first?
    allocate(At(lda,n), wr(n), wi(n), vl(ldvl,n), vr(ldvr,n), work(lwork), lam(n))
    At = A

    call dgeev('N', 'N', n, At, lda, wr, wi, vl, ldvl, vr, ldvr, &
         work, lwork, info)
    if(info /= 0) then
       print *, "dgeev returned info = ", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "the QR algorithm failed to compute all the"
          print *, "eigenvalues, and no eigenvectors have been computed;"
          print *, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
          print *, "have converged."
       end if
       call stop_error('eigvals: dgeev error')
    end if

    lam = wr + i_*wi
  end function deigvals

  function zeigvals(A) result(lam)
    complex(dp), intent(in) :: A(:, :)  ! matrix to solve eigenproblem for
    complex(dp), allocatable :: lam(:)  ! eigenvalues: A c = lam c
    ! LAPACK variables:
    integer :: info, lda, ldvl, ldvr, lwork, n, lrwork
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: At(:,:), vl(:,:), vr(:,:), work(:)

    lda = size(A(:,1))
    n = size(A(1,:))
    call assert_shape(A, [n, n], "solve", "A")
    ldvl = n
    ldvr = n
    lwork = 8*n  ! TODO: can this size be optimized? query first?
    lrwork = 2*n
    allocate(At(lda,n), vl(ldvl,n), vr(ldvr,n), work(lwork), rwork(lrwork), lam(n))
    At = A
    call zgeev('N', 'N', n, At, lda, lam, vl, ldvl, vr, ldvr, work, &
         lwork, rwork, info)
    if(info /= 0) then
       print *, "zgeev returned info = ", info
       if(info < 0) then
          print *, "the ",-info, "-th argument had an illegal value."
       else
          print *, "the QR algorithm failed to compute all the"
          print *, "eigenvalues, and no eigenvectors have been computed;"
          print *, "elements and ", info+1, ":", n, " of W contain eigenvalues which have"
          print *, "converged."
       end if
       call stop_error('eig: zgeev error')
    end if
  end function zeigvals

```

---

```fortran
program test_inv
use types, only: dp
use constants, only: i_
use utils, only: assert
use linalg, only: inv, eye
implicit none

real(dp), parameter :: eps = 1e-9_dp
complex(dp) :: A(2, 2), B(2, 2), C(2, 2)
real(dp) :: D(2,2), E(2,2), F(2,2)

! test for complex double routines:
A = reshape([0+0*i_, i_, 2+0*i_, 0+0*i_], [2, 2], order=[2, 1])
C = inv(A)
B = reshape([0+0*i_, 1+0*i_, -2*i_, 0+0*i_], [2, 2], order=[2, 1]) / 2
call assert(maxval(abs(C-B)) < eps)
call assert(maxval(abs(matmul(C, A)-eye(2))) < eps)
call assert(maxval(abs(matmul(A, C)-eye(2))) < eps)

A = reshape([1+0*i_, i_, 2+0*i_, 3*i_], [2, 2], order=[2, 1])
C = inv(A)
B = reshape([3+0*i_, -1+0*i_, 2*i_, -i_], [2, 2], order=[2, 1])
call assert(maxval(abs(C-B)) < eps)
call assert(maxval(abs(matmul(C, A)-eye(2))) < eps)
call assert(maxval(abs(matmul(A, C)-eye(2))) < eps)

A = reshape([0+0*i_, i_, -i_, 0+0*i_], [2, 2], order=[2, 1])
C = inv(A)
call assert(maxval(abs(matmul(C, A)-eye(2))) < eps)
call assert(maxval(abs(matmul(A, C)-eye(2))) < eps)

! tests for double precision routine:
D = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], shape=[2,2])
E = reshape([-2.0_dp, 1.0_dp, 1.5_dp, -0.5_dp], shape=[2,2])
F = inv(D)
call assert(maxval(abs(E - F)) < eps)
call assert(maxval(abs(matmul(F, D) - cmplx(eye(2)))) < eps)
call assert(maxval(abs(matmul(E, D) - cmplx(eye(2)))) < eps)

D = reshape([0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], shape=[2,2])  ! is its own inverse
F = inv(D)
call assert(maxval(abs(D - F)) < eps)
call assert(maxval(abs(matmul(F, D) - cmplx(eye(2)))) < eps)

end program
```

---

```fortran
  ! matrix inversion for real/complex matrices:
  interface inv
     module procedure dinv
     module procedure zinv
  end interface inv


  function dinv(Am) result(Bm)
    real(dp), intent(in) :: Am(:,:)  ! matrix to be inverted
    real(dp) :: Bm(size(Am, 1), size(Am, 2))   ! Bm = inv(Am)
    real(dp), allocatable :: Amt(:,:), work(:)  ! temporary work arrays

    ! LAPACK variables:
    integer ::  info, lda, n, lwork, nb
    integer, allocatable :: ipiv(:)

    ! use LAPACK's dgetrf and dgetri
    n = size(Am(1, :))
    call assert_shape(Am, [n, n], "inv", "Am")
    lda = n
    nb = ilaenv(1, 'DGETRI', "UN", n, -1, -1, -1)  ! TODO: check UN param
    lwork = n*nb
    if (nb < 1) nb = max(1, n)
    allocate(Amt(n,n), work(lwork), ipiv(n))
    Amt = Am
    call dgetrf(n, n, Amt, lda, ipiv, info)
    if(info /= 0) then
       print *, "dgetrf returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, and division by zero will occur if it is used"
          print *, "to solve a system of equations."
       end if
       call stop_error('inv: dgetrf error')
    end if

    call dgetri(n, Amt, n, ipiv, work, lwork, info)
    if (info /= 0) then
       print *, "dgetri returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; the matrix is"
          print *, "singular and its inverse could not be computed."
       end if
       call stop_error('inv: dgetri error')
    end if
    Bm = Amt

  end function dinv

  function zinv(Am) result(Bm)
    ! Inverts the general complex matrix Am
    complex(dp), intent(in) :: Am(:,:)   ! Matrix to be inverted
    complex(dp) :: Bm(size(Am, 1), size(Am, 2))   ! Bm = inv(Am)
    integer :: n, nb
    ! lapack variables
    integer :: lwork, info
    complex(dp), allocatable:: Amt(:,:), work(:)
    integer, allocatable:: ipiv(:)

    n = size(Am, 1)
    call assert_shape(Am, [n, n], "inv", "Am")
    nb = ilaenv(1, 'ZGETRI', "UN", n, -1, -1, -1)  ! TODO: check UN param
    if (nb < 1) nb = max(1, n)
    lwork = n*nb
    allocate(Amt(n,n), ipiv(n), work(lwork))
    Amt = Am
    call zgetrf(n, n, Amt, n, ipiv, info)
    if (info /= 0) then
       print *, "zgetrf returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, and division by zero will occur if it is used"
          print *, "to solve a system of equations."
       end if
       call stop_error('inv: zgetrf error')
    end if
    call zgetri(n, Amt, n, ipiv, work, lwork, info)
    if (info /= 0) then
       print *, "zgetri returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; the matrix is"
          print *, "singular and its inverse could not be computed."
       end if
       call stop_error('inv: zgetri error')
    end if
    Bm = Amt
  end function zinv

```

---

```fortran
program test_lstsq
use types, only: dp
use utils, only: assert
use linalg, only: lstsq, eye
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(10,10), b(10), x(10)
real(dp) :: AA(5,3), bb(5), xx(3)
complex(dp) :: AC(5,3), bc(5), xc(3)
complex(dp) :: C(3,3), d(3), y(3)

! solve square systems of equations, recycle the solve() tests:

! test the eye*x = 0 with solution x = 0
A = eye(10)
b = 0.0_dp
x = lstsq(A,b)
call assert(maxval(x) < eps)

! test eye*x = 1, with solution x(:) = 1
b = 1.0_dp
x = lstsq(A,b)
call assert(maxval(abs(x-1.0_dp)) < eps)

! test i*eye*y = -1 with solution y(:) = i
d = cmplx(-1.0_dp + 0*i_)
C = i_*cmplx(eye(3))
y = lstsq(C, d)
call assert(maxval(abs(y - i_)) < eps)

! perform least squares computations:

! test overdetermined systems:

! use model y = ax^2 + bx + c; choose 5 points (x=0,\pm 1,\pm 2) exactly on the parabola with a=b=c=1
! expected solution: [a, b, c] = [1, 1, 1]
AA = reshape([0, 0, 1, 1, -1, 1, 1, 1, 1, 4, -2, 1, 4, 2, 1], order=[2,1], shape=[5,3])
bb = [1, 1, 3, 3, 7]  ! RHS of eq. for x=0,-1,+1,-2,+2
xx = lstsq(AA, bb)
call assert(maxval(abs(xx - [1, 1, 1])) < eps)

! use same model above, but now use complex data: a=b=1, c=i
AC = reshape([0, 0, 1, 1, -1, 1, 1, 1, 1, 4, -2, 1, 4, 2, 1], order=[2,1], shape=[5,3])
bc = [i_, i_, 2+i_, 2+i_, 6+i_]  ! RHS of eq. for x=0,-1,+1,-2,+2
xc = lstsq(AC, bc)
call assert(maxval(abs(xc - [1+0*i_, 1+0*i_, i_])) < eps)

end program
```

---

```fortran
  ! least square solutions the real/complex systems of equations of possibly non-square shape:
  interface lstsq
     module procedure dlstsq
     module procedure zlstsq
  end interface lstsq

  function dlstsq(A, b) result(x)
    ! compute least square solution to A x = b for real A, b
    real(dp), intent(in) :: A(:,:), b(:)
    real(dp), allocatable :: x(:)
    ! LAPACK variables:
    integer :: info, ldb, lwork, m, n, rank
    real(dp) :: rcond
    real(dp), allocatable :: work(:), At(:,:), Bt(:,:)
    integer, allocatable :: jpvt(:)

    m = size(A(:,1)) ! = lda
    n = size(A(1,:))
    ldb = size(b)
    allocate(x(n), At(m,n), Bt(ldb,1), jpvt(n), work(1))
    call dgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
         -1, info)  ! query optimal workspace size
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))  ! allocate with ideal size
    rcond = 0.0_dp
    jpvt(:) = 0
    Bt(:,1) = b(:)  ! only one right-hand side
    At(:,:) = A(:,:)
    call dgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
         lwork, info)
    if(info /= 0) then
       print *, "dgelsy returned info = ", info
       print *, "the ", -info, "-th argument had an illegal value"
       call stop_error('lstsq: dgelsy error')
    endif
    x(:) = Bt(1:n,1)
  end function dlstsq

  function zlstsq(A, b) result(x)
    ! compute least square solution to A x = b for complex A, b
    complex(dp), intent(in) :: A(:,:), b(:)
    complex(dp), allocatable :: x(:)
    ! LAPACK variables:
    integer :: info, ldb, lwork, m, n, rank
    real(dp) :: rcond
    complex(dp), allocatable :: At(:,:), Bt(:,:), work(:)
    real(dp), allocatable :: rwork(:)
    integer, allocatable :: jpvt(:)

    m = size(A(:,1)) ! = lda
    n = size(A(1,:))
    ldb = size(b)
    allocate(x(n), At(m,n), Bt(ldb,1), jpvt(n), work(1), rwork(2*n))
    call zgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
         -1, rwork, info)  ! query optimal workspace size
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))  ! allocate with ideal size
    rcond = 0.0_dp
    jpvt(:) = 0
    Bt(:,1) = b(:)  ! only one right-hand side
    At(:,:) = A(:,:)
    call zgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
         lwork, rwork, info)
    if(info /= 0) then
       print *, "zgelsy returned info = ", info
       print *, "the ", -info, "-th argument had an illegal value"
       call stop_error('lstsq: zgelsy error')
    endif
    x(:) = Bt(1:n,1)
  end function zlstsq
```

---

```fortran
program test_solve
use types, only: dp
use utils, only: assert
use linalg, only: solve, eye
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(10,10), b(10), x(10)
complex(dp) :: C(3,3), d(3), y(3)

! test the eye*x = 0 with solution x = 0
A = eye(10)
b = 0.0_dp
x = solve(A,b)
call assert(maxval(x) < eps)

! test eye*x = 1, with solution x(:) = 1
b = 1.0_dp
x = solve(A,b)
call assert(maxval(abs(x-1.0_dp)) < eps)

! test i*eye*y = -1 with solution y(:) = i
d = cmplx(-1.0_dp + 0*i_)
C = i_*cmplx(eye(3))
y = solve(C, d)
call assert(maxval(abs(y - i_)) < eps)

! TODO: implement complex tests
end program
```

---

```fortran
  ! solution to linear systems of equation with real/complex coefficients:
  interface solve
     module procedure dsolve
     module procedure zsolve
  end interface solve


  function dsolve(A, b) result(x)
    ! solves a system of equations A x = b with one right hand side
    real(dp), intent(in) :: A(:,:)  ! coefficient matrix A
    real(dp), intent(in) :: b(:)  ! right-hand-side A x = b
    real(dp), allocatable :: x(:)
    ! LAPACK variables:
    real(dp), allocatable :: At(:,:), bt(:,:)
    integer :: n, info, lda
    integer, allocatable :: ipiv(:)

    n = size(A(1,:))
    lda = size(A(:, 1))  ! TODO: remove lda (which is = n!)
    call assert_shape(A, [n, n], "solve", "A")
    allocate(At(lda,n), bt(n,1), ipiv(n), x(n))
    At = A
    bt(:,1) = b(:)
    call dgesv(n, 1, At, lda, ipiv, bt, n, info)
    if(info /= 0) then
       print *, "dgesv returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, so the solution could not be computed."
       end if
       call stop_error('inv: dgesv error')
    endif
    x = bt(:,1)
  end function dsolve

  function zsolve(A, b) result(x)
    ! solves a system of equations A x = b with one right hand side
    complex(dp), intent(in) :: A(:,:)  ! coefficient matrix A
    complex(dp), intent(in) :: b(:)  ! right-hand-side A x = b
    complex(dp), allocatable :: x(:)
    ! LAPACK variables:
    complex(dp), allocatable :: At(:,:), bt(:,:)
    integer :: n, info, lda
    integer, allocatable :: ipiv(:)

    n = size(A(1,:))
    lda = size(A(:, 1))  ! TODO: remove lda here, too
    call assert_shape(A, [n, n], "solve", "A")
    allocate(At(lda,n), bt(n,1), ipiv(n), x(n))
    At = A
    bt(:,1) = b(:)
    call zgesv(n, 1, At, lda, ipiv, bt, n, info)
    if(info /= 0) then
       print *, "zgesv returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, so the solution could not be computed."
       end if
       call stop_error('inv: zgesv error')
    endif
    x = bt(:,1)
  end function zsolve
```

---

```fortran
program test_svd
use types, only: dp
use utils, only: assert
use linalg, only: svd, eye, diag
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: B(3, 2), s(2), sc(2), U(3, 3), Vtransp(2, 2), sigma(3, 2)
complex(dp) :: AC(3, 2), Uc(3, 3), Vtranspc(2, 2)

! test if the matrix' reconstruction from the SVD is faithful

AC = reshape([1+i_, i_, 0*i_, -i_, 2*i_, 2+0*i_], shape=[3,2], order=[2,1])
call svd(AC, sc, Uc, Vtranspc)
sigma = 0.0_dp
sigma(:2,:2) = diag(sc)
call assert(maxval(abs(AC - matmul(Uc, matmul(sigma, Vtranspc)))) < eps)

B = reshape([1, 2, 3, 4, 5, 6], shape=[3,2], order=[2,1])
call svd(B, s, U, Vtransp)
sigma(:2,:2) = diag(s)
call assert(maxval(abs(B - matmul(U, matmul(sigma, Vtransp)))) < eps)

end program
```

---

```fortran
  ! singular value decomposition of real/complex matrices:
  interface svd
     module procedure dsvd
     module procedure zsvd
  end interface svd

  subroutine dsvd(A, s, U, Vtransp)
    ! compute the singular value decomposition A = U sigma Vtransp of a
    ! real m x n matrix A
    ! U is m x m
    ! Vtransp is n x n
    ! s has size min(m, n) --> sigma matrix is (n x m) with sigma_ii = s_i
    real(dp), intent(in) :: A(:,:)
    real(dp), intent(out) :: s(:), U(:,:), Vtransp(:,:)
    ! LAPACK related:
    integer :: info, lwork, m, n, ldu
    real(dp), allocatable :: work(:), At(:,:)

    ! TODO: check shapes here and in other routines?

    m = size(A(:,1))  ! = lda
    n = size(A(1,:))
    ldu = m
    allocate(At(m,n))
    At(:,:) = A(:,:)  ! use a temporary as dgesvd destroys its input

    call assert_shape(U, [m, m], "svd", "U")
    call assert_shape(Vtransp, [n, n], "svd", "Vtransp")

    ! query optimal lwork and allocate workspace:
    allocate(work(1))
    call dgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, -1, info)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))

    call dgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, lwork, info)
    if(info /= 0) then
       print *, "dgesvd returned info = ", info
       if(info < 0) then
          print *, "the ", -info, "-th argument had an illegal value"
       else
          print *, "DBDSQR did not converge, there are ", info
          print *, "superdiagonals of an intermediate bidiagonal form B"
          print *, "did not converge to zero. See the description of WORK"
          print *, "in DGESVD's man page for details."
       endif
       call stop_error('svd: dgesvd error')
    endif
  end subroutine dsvd

  subroutine zsvd(A, s, U, Vtransp)
    ! compute the singular value decomposition A = U sigma V^H of a
    ! complex m x m matrix A
    ! U is m x min(m, n)
    ! Vtransp is n x n
    ! sigma is m x n with with sigma_ii = s_i
    ! note that this routine returns V^H, not V!
    complex(dp), intent(in) :: A(:,:)
    real(dp), intent(out) :: s(:)
    complex(dp), intent(out) :: U(:,:), Vtransp(:,:)
    ! LAPACK related:
    integer :: info, lwork, m, n, ldu, lrwork
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: work(:), At(:,:)

    ! TODO: check shapes here and in other routines?

    m = size(A(:,1))  ! = lda
    n = size(A(1,:))
    ldu = m
    lrwork = 5*min(m,n)
    allocate(rwork(lrwork), At(m,n))
    At(:,:) = A(:,:)  ! use a temporary as zgesvd destroys its input

    call assert_shape(U, [m, m], "svd", "U")
    call assert_shape(Vtransp, [n, n], "svd", "Vtransp")

    ! query optimal lwork and allocate workspace:
    allocate(work(1))
    call zgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, -1,&
         rwork, info)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))

    call zgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, &
         lwork, rwork, info)
    if(info /= 0) then
       print *, "zgesvd returned info = ", info
       if(info < 0) then
          print *, "the ", -info, "-th argument had an illegal value"
       else
          print *, "ZBDSQR did not converge, there are ", info
          print *, "superdiagonals of an intermediate bidiagonal form B"
          print *, "did not converge to zero. See the description of WORK"
          print *, "in DGESVD's man page for details."
       endif
       call stop_error('svd: zgesvd error')
    endif
  end subroutine zsvd
```

---

```fortran
program test_svdvals
use types, only: dp
use utils, only: assert
use linalg, only: svdvals, eye
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: B(3, 2), s(2), sc(2)
complex(dp) :: AC(3, 2)

! test a general complex matrix
AC = reshape([1+i_, i_, 0*i_, -i_, 2*i_, 2+0*i_], shape=[3,2], order=[2,1])
sc = svdvals(AC)
! svdvals from SciPy: 3.0269254467476365, 1.6845540477620835
call assert(maxval(abs(sc - [3.0269254467476365_dp, 1.6845540477620835_dp])) < eps)

B = reshape([1, 2, 3, 4, 5, 6], shape=[3,2], order=[2,1])
s = svdvals(B)
! svdvals from SciPy: 9.5255180915651092, 0.51430058065864448
call assert(maxval(abs(s - [9.5255180915651092_dp, 0.51430058065864448_dp])) < eps)

! note that SciPy uses the same LAPACK routines. thus, comparing with 
! SciPy's results does not test correctness of the results, instead, the
! wrapper is tested.
end program
```

---

```fortran
  ! singular values of real/complex matrices:
  interface svdvals
     module procedure dsvdvals
     module procedure zsvdvals
  end interface svdvals
  function dsvdvals(A) result(s)
    ! compute singular values s_i of a real m x n matrix A
    real(dp), intent(in) :: A(:,:)
    real(dp), allocatable :: s(:)
    ! LAPACK related:
    integer :: info, lwork, m, n
    real(dp), allocatable :: work(:), At(:,:)
    real(dp) :: u(1,1), vt(1,1)  ! not used if only s is to be computed

    m = size(A(:,1))  ! = lda
    n = size(A(1,:))
    allocate(At(m,n), s(min(m,n)))
    At(:,:) = A(:, :)  ! A is overwritten in dgesvd

    ! query optimal lwork and allocate workspace:
    allocate(work(1))
    call dgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, -1, info)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))

    call dgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, lwork, info)
    if(info /= 0) then
       print *, "dgesvd returned info = ", info
       if(info < 0) then
          print *, "the ", -info, "-th argument had an illegal value"
       else
          print *, "DBDSQR did not converge, there are ", info
          print *, "superdiagonals of an intermediate bidiagonal form B"
          print *, "did not converge to zero. See the description of WORK"
          print *, "in DGESVD's man page for details."
       endif
       call stop_error('svdvals: dgesvd error')
    endif
  end function dsvdvals

  function zsvdvals(A) result(s)
    ! compute singular values s_i of a real m x n matrix A
    complex(dp), intent(in) :: A(:,:)
    real(dp), allocatable :: s(:)
    ! LAPACK related:
    integer :: info, lwork, m, n, lrwork
    complex(dp), allocatable :: work(:), At(:,:)
    real(dp), allocatable :: rwork(:)
    complex(dp) :: u(1,1), vt(1,1)  ! not used if only s is to be computed

    m = size(A(:,1))  ! = lda
    n = size(A(1,:))
    lrwork = 5*min(m,n)
    allocate(At(m,n), s(min(m,n)), rwork(lrwork))
    At(:,:) = A(:,:)  ! A is overwritten in zgesvd!

    ! query optimal lwork and allocate workspace:
    allocate(work(1))
    call zgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, -1, rwork, info)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))

    call zgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, lwork, rwork, info)
    if(info /= 0) then
       print *, "zgesvd returned info = ", info
       if(info < 0) then
          print *, "the ", -info, "-th argument had an illegal value"
       else
          print *, "ZBDSQR did not converge, there are ", info
          print *, "superdiagonals of an intermediate bidiagonal form B"
          print *, "did not converge to zero. See the description of RWORK"
          print *, "in ZGESVD's man page for details."
       endif
       call stop_error('svdvals: zgesvd error')
    endif
  end function zsvdvals
```

---

```fortran
program test_trace
use types, only: dp
use utils, only: assert
use linalg, only: trace
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(3,3)
complex(dp) :: B(3,3)

A = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9], shape=[3,3], order=[2,1])
call assert(abs(trace(A) - 15.0_dp) < eps)

B = i_*reshape([1, 2, 3, 4, 5, 6, 7, 8, 9], shape=[3,3], order=[2,1])
call assert(abs(trace(B) - 15*i_) < eps)

end program
```

---

```fortran
  ! trace of real/complex matrices:
  interface trace
     module procedure dtrace
     module procedure ztrace
  end interface trace


  ! TODO: add optional axis parameter in both xtrace() functions
  function dtrace(A) result(t)
    ! return trace along the main diagonal
    real(dp), intent(in) :: A(:,:)
    real(dp) :: t
    integer :: i

    t = 0.0_dp
    do i = 1,minval(shape(A))
       t = t + A(i,i)
    end do
  end function dtrace

  function ztrace(A) result(t)
    ! return trace along the main diagonal
    complex(dp), intent(in) :: A(:,:)
    complex(dp) :: t
    integer :: i

    t = 0*i_
    do i = 1,minval(shape(A))
       t = t + A(i,i)
    end do
  end function ztrace
```

---

#### まとめ

