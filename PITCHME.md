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


### 今回は線形代数部分をkwsk⚠

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

---

```fortran
```

---

```fortran
  ! determinants of real/complex square matrices:
  interface det
     module procedure ddet
     module procedure zdet
  end interface det
```

---

```fortran
```

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

---

```fortran
```

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
```

---

#### まとめ

