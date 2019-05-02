# FortranをPythonに近づけろ
## fortran-utilsを使ってみた

---


### はじめに

なぜその言語を使うの?


+++


### はじめに

言語の仕様ではなくライブラリが使えるからでは？


+++


### はじめに

<img src="https://3.bp.blogspot.com/-CR-N5Kb5OOw/WoO1nOR-xHI/AAAAAAAAXKw/WJ0tWg108Zw_Va8Z1JCpnIHrDZz5YwSJQCLcBGAs/s1600/logo-stack-python.png" width="600" height="300">

が使えなかったらPythonを使うだろうか？


---


### Fortranを使ってもらうには

Numpy/Scipyみたいな仕様のライブラリがあればよいのではないか。


+++


### Fortranを使ってもらうには

今回はそれに近いfortran-utilsを使ってみます。


---


### fortran-utilsとは？

[certik](https://github.com/certik)氏により2016年まで開発されていたfortranの便利なutils集です。

https://github.com/certik/fortran-utils

一部テストが通らない部分があったのでforkして変更をしました。

https://github.com/tkoyama010/fortran-utils

---


### 全体構成

* Types (``dp``)
* Constants (``pi``, ``e_``, ``i_``)
* Sorting
* Saving/loading 2D arrays (``savetxt``, ``loadtxt``)
* Meshes (exponential, uniform)
* Cubic splines
* Saving/loading PPM images
* Lapack interface (and a few simple f90 wrappers like ``eigh``, ``inv``)
* HDF5 interface


---


### Types (``dp``)


---


### Constants in Numpy/Scipy

Euler's equation $e^{i\pi}+1=0$

```
import numpy as np

print(abs(np.e**(np.pi*1j) + 1) < 10**(-5))

```


+++


### Constants in fortran-utils

Euler's equation $e^{i\pi}+1=0$

```
program test_constants
use types, only: dp
use constants, only: pi, e_, i_
use utils, only: assert
implicit none

! Euler's identity:
call assert(abs(e_**(pi*i_) + 1) < 1e-15_dp)

end program
```


---


### Sorting


```
program test_sort
use types, only: dp
use utils, only: assert
use sorting, only: sort, sortpairs
implicit none
integer :: a(5)
real(dp) :: b(5), c(5), vec(2, 5)
a = [4, 3, 2, 1, 5]
call sort(a)
call assert(all(a == [1, 2, 3, 4, 5]))
end program
```


---


### Sorting


```
program test_argsort
use types, only: dp
use utils, only: assert
use sorting, only: argsort
implicit none
integer :: a(5)
a = [4, 3, 2, 1, 5]
call assert(all(argsort(a) == [4, 3, 2, 1, 5]))
end program
```


---


### Saving/loading 2D arrays




---


### Saving/loading 2D arrays


```
program test_loadtxt
use types, only: dp
use utils, only: loadtxt, savetxt, assert
implicit none

real(dp) :: d(3, 2), e(2, 3)
real(dp), allocatable :: d2(:, :)
d = reshape([1, 2, 3, 4, 5, 6], [3, 2])
call savetxt("tmp.dat", d)
call loadtxt("tmp.dat", d2)
call assert(all(shape(d2) == [3, 2]))
call assert(all(d == d2))

e = reshape([1, 2, 3, 4, 5, 6], [2, 3])
call savetxt("tmp.dat", e)
call loadtxt("tmp.dat", d2)
call assert(all(shape(d2) == [2, 3]))
call assert(all(e == d2))

end program
```


---


### Lapack interface (det)


---


### Lapack interface (det)


```
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


### Meshes (exponential, uniform)


---


### Meshes (exponential, uniform)


---


### Cubic splines


---


### Cubic splines


---


### Saving/loading PPM images


---


### Saving/loading PPM images


---


### HDF5 interface


---


### HDF5 interface


---


### まとめ

言語の表現を近づければみんなが使ってくれる？

<img src="https://upload.wikimedia.org/wikipedia/commons/2/29/Pieter_Bruegel_the_Elder_-_The_Tower_of_Babel_%28Rotterdam%29_-_Google_Art_Project.jpg" width="400" height="400">
