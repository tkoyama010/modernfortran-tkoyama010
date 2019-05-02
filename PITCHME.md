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

<img src="https://www.fullstackpython.com/img/logos/numpy.jpg" width="200" height="100">

<img src="https://www.fullstackpython.com/img/logos/scipy.png" width="200" height="100">

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


### Constants (``pi``, ``e_``, ``i_``)

- Eulerの公式 $e^{i\pi}+1=0$

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
@[8]


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
call assert(all(argsort([4, 3, 2, 1, 5]) == [4, 3, 2, 1, 5]))
call assert(all(argsort([10, 9, 8, 7, 6]) == [5, 4, 3, 2, 1]))
call assert(all(argsort([1, -1]) == [2, 1]))
call assert(all(argsort([1, 2, 2, 2, 3]) == [1, 2, 3, 4, 5]))
call assert(all(argsort([2, 2, 2, 3, 1]) == [5, 2, 3, 1, 4]))
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


### Lapack interface (and a few simple f90 wrappers like ``eigh``, ``inv``)


---


### Lapack interface (and a few simple f90 wrappers like ``eigh``, ``inv``)


---


### HDF5 interface


---


### HDF5 interface


---


### まとめ

- Numpy/Scipyは数値計算系がPythonを使う主な理由です。
- Fortranでもfortran-utilsを使えばNumpy/Scipyと同じことができる？
- fortran-utilsを充実させればみんながFortranを使ってくれる！？
- というかもうFortranの文法をPythonと同一に進化させればいい（この行はネタです）。
