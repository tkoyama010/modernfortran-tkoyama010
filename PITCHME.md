## fortran-utilsコードリーディング
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


#### まとめ

