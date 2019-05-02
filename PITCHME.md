## FortranをPythonに
## 近づけろ
### fortran-utilsを使ってみた
@tkoyama010

---


### Introduction

なぜその言語を使うの?


+++


### Introduction

言語の仕様ではなくライブラリが使えるからでは？


+++?https://3.bp.blogspot.com/-CR-N5Kb5OOw/WoO1nOR-xHI/AAAAAAAAXKw/WJ0tWg108Zw_Va8Z1JCpnIHrDZz5YwSJQCLcBGAs/s1600/logo-stack-python.png


### Introduction


ライブラリが使えなかったらPythonを使うだろうか？
(画像はhttps://3.bp.blogspot.com より)


---


### Numpy/Scipyみたいな仕様のライブラリがあればよいのではないか。


+++


### 今回はそれに近いfortran-utilsを使ってみます。


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
import numpy as np
a = [4, 3, 2, 1, 5]
a = np.sort(a)
print(np.all(a == [1, 2, 3, 4, 5]))
```


+++


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


```
import numpy as np

d = np.array([[1, 3, 5], [2, 4, 6]])
np.savetxt("tmp.dat", d)
d2 = np.loadtxt("tmp.dat")
print(np.all(d == d2))
print(d.shape == d2.shape)
```


+++


### Saving/loading 2D arrays


```
program test_loadtxt
use types, only: dp
use utils, only: loadtxt, savetxt, assert
implicit none

real(dp) :: d(3, 2)
real(dp), allocatable :: d2(:, :)
d = reshape([1, 2, 3, 4, 5, 6], [3, 2])
call savetxt("tmp.dat", d)
call loadtxt("tmp.dat", d2)
call assert(all(shape(d2) == [3, 2]))
call assert(all(d == d2))

end program
```


---


### Lapack interface (det)

```
import numpy as np
import numpy.linalg as LA

A = [[ 1, 4, 7], [ 2, 5, 8], [ 3, 6, -9]]
detA = LA.det(A)
print(abs(detA - 54.0) < 1.0e-09)
```

+++


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

A = real(reshape([ 1, 2, 3, 4, 5, 6, 7, 8, -9 ], shape=[3,3]))
Adet = det(A)
call assert(abs(Adet - 54.0_dp) < eps)

end program
```


---


### eigen value problem

```
import numpy
import numpy.linalg as LA

B = 3.0*np.eye(5)
lamb, cb = LA.eig(B)
print(np.max((lamb-3) < 1.0e-09))
print(np.max((cb-np.eye(5)) < 1.0e-09))
```


+++


### eigen value problem

```
program test_eig
use types, only: dp
use utils, only: assert
use linalg, only: eig, eye
use constants, only : i_
implicit none

! test eigenvalue comutation for general matrices:

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: B(5, 5)
complex(dp) :: lamb(5), cb(5, 5)

! test a multiple of the unit matrix:
B = 3*eye(5)
call eig(B, lamb, cb)
call assert(maxval(abs(lamb - 3.0_dp)) < eps)  ! all eigenvalues are 3
call assert(maxval(abs(cb - cmplx(eye(5)))) < eps)  ! eigenvectors are cartesian unit basis vectors

end program
```


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


### まとめ?https://upload.wikimedia.org/wikipedia/commons/2/29/Pieter_Bruegel_the_Elder_-_The_Tower_of_Babel_%28Rotterdam%29_-_Google_Art_Project.jpg

言語の表現を近づければみんなが使ってくれる？

(画像はhttps://upload.wikimedia.org より)
