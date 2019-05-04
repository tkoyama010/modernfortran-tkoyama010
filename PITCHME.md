## FortranをPythonに
## 近づけろ
### fortran-utilsを使ってみた
@tkoyama010

---


### なぜその言語を使うの?


---


### 言語の仕様ではなくライブラリが使えるからでは？


---


### ライブラリが使えなかったらPythonを使うだろうか？


---


### Numpy/Scipyみたいな仕様のライブラリがあればよいのではないか。


---


### 今回はそれに近いfortran-utilsを使ってみた。


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


### Sorting in Numpy/Scipy


```
import numpy as np
a = [4, 3, 2, 1, 5]
a = np.sort(a)
print(np.all(a == [1, 2, 3, 4, 5]))
```


+++


### Sorting in fortran-utils


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


### ArgSorting in Numpy/Scipy


```
import numpy as np
a = [4, 3, 2, 1, 5]
print(np.all(np.argsort(a) == [3, 2, 1, 0, 4]))
```


+++


### ArgSorting in fortran-utils


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


### Saving/loading 2D arrays in Numpy/Scipy


```
import numpy as np

d = np.array([[1, 3, 5], [2, 4, 6]])
np.savetxt("tmp.dat", d)
d2 = np.loadtxt("tmp.dat")
print(np.all(d == d2))
print(d.shape == d2.shape)
```


+++


### Saving/loading 2D arrays in fortran-utils


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


### Lapack interface (det) in Numpy/Scipy

```
import numpy as np
import numpy.linalg as LA

A = [[ 1, 4, 7], [ 2, 5, 8], [ 3, 6, -9]]
detA = LA.det(A)
print(abs(detA - 54.0) < 1.0e-09)
```

+++


### Lapack interface (det) in fortran-utils


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


### diag in Numpy/Scipy


```
import numpy as np

eps = 1.0e-9

A = np.diag([1.0, 1.0, 1.0])
print(np.max(A-np.eye(3)) < eps)
```


+++


### diag in fortran-utils


```
program test_diag
use types, only: dp
use utils, only: assert
use linalg, only: diag, eye
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(3,3)

A = diag([1.0_dp, 1.0_dp, 1.0_dp])
call assert(maxval(abs(A-eye(3))) < eps)

end program
```


---


### eigen value problem in Numpy/Scipy

```
import numpy
import numpy.linalg as LA

eps = 1.0e-09
B = 3.0*np.eye(5)
lamb, cb = LA.eig(B)
print(np.max((lamb-3) < eps))
print(np.max((cb-np.eye(5)) < eps))
```


+++


### eigen value problem in fortran-utils

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


### inv in Numpy/Scipy


```
import numpy as np
import numpy.linalg as LA

eps = 1.0e-09

D = np.array([[0.0, 1.0], [1.0, 0.0]])
F = LA.inv(D)
print(np.max(abs(D.dot(F)-np.eye(2))) < eps)
```


+++


### inv in fortran-utils


```
program test_inv
use types, only: dp
use utils, only: assert
use linalg, only: inv, eye
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: D(2,2), F(2,2)

D = reshape([0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], shape=[2,2])  ! is its own inverse
F = inv(D)
call assert(maxval(abs(matmul(F, D) - cmplx(eye(2)))) < eps)

end program
```


---


### solve in Numpy/Scipy


```
import numpy as np
import numpy.linalg as LA

eps = 1.0e-09

d = [-1.0, -1.0, -1.0]
C = 1j*np.eye(3)
y = LA.solve(C, d)
print(np.max(np.abs(y-1j)) < eps)
```


+++


### solve in fortran-utils


```
program test_solve
use types, only: dp
use utils, only: assert
use linalg, only: solve, eye
use constants, only : i_
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(10,10), b(10), x(10)
complex(dp) :: C(3,3), d(3), y(3)

! test i*eye*y = -1 with solution y(:) = i
d = cmplx(-1.0_dp + 0*i_)
C = i_*cmplx(eye(3))
y = solve(C, d)
call assert(maxval(abs(y - i_)) < eps)

end program
```

---

#### まとめ

<img src="https://upload.wikimedia.org/wikipedia/commons/2/29/Pieter_Bruegel_the_Elder_-_The_Tower_of_Babel_%28Rotterdam%29_-_Google_Art_Project.jpg" width="400" height="400">

(https://upload.wikimedia.org より)
