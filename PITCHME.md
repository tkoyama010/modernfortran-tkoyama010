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
  ! determinants of real/complex square matrices:
  interface det
     module procedure ddet
     module procedure zdet
  end interface det
```

---

#### まとめ

