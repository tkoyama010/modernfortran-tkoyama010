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
@[8](Eulerの公式)
$e^{ i\pi } + 1 = 0$


---


### Sorting


---


### Sorting


---


### Saving/loading 2D arrays (``savetxt``, ``loadtxt``)


---


### Saving/loading 2D arrays (``savetxt``, ``loadtxt``)


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
