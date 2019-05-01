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

![Numpy](https://camo.githubusercontent.com/b3ffb735f64eea19dd6790720c0d7e8db71931aa/68747470733a2f2f63646e2e7261776769742e636f6d2f6e756d70792f6e756d70792f6d61737465722f6272616e64696e672f69636f6e732f6e756d70796c6f676f2e737667)
![Sicpy](https://www.fullstackpython.com/img/logos/scipy.png)

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
