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


#### 今回は一部を詳しくコードリーディングします

* Types (``dp``)
* Constants (``pi``, ``e_``, ``i_``)
* [Sorting](https://github.com/tkoyama010/fortran-utils/blob/master/src/sorting.f90)
* Saving/loading 2D arrays (``savetxt``, ``loadtxt``)
* Meshes (exponential, uniform)
* Cubic splines
* Saving/loading PPM images
* [Lapack interface (and a few simple f90 wrappers like ``eigh``, ``inv``)](https://github.com/tkoyama010/fortran-utils/blob/master/src/linalg.f90)
* HDF5 interface

---

```fortran
a = [4, 3, 2, 1, 5]
call sort(a)
call assert(all(a == [1, 2, 3, 4, 5]))
a = [5, 4, 3, 2, 1]
call sort(a)
call assert(all(a == [1, 2, 3, 4, 5]))

b = [5, 4, 3, 2, 1]
call sort(b)
call assert(all(b == [1, 2, 3, 4, 5]))
```
@[2,3,4](まずはsort関数について解説します。)
@[3,4,5](テストコードではasser文を使っていますが、これは論理型を引数にした関数として定義されています。)

---

```fortran
! overload sort
interface sort
    module procedure sortNums, sortINums, sortVecs
end interface

subroutine sortNums(nums)
! sorts array of numbers, nums, from smallest to largest
real(dp), intent(inout):: nums(:)   ! array of numbers
nums = nums(argsort(nums))
end subroutine

subroutine sortINums(nums)
! sorts array of inegers, nums, from smallest to largest
integer, intent(inout):: nums(:)    ! array of numbers
nums = nums(argsort(nums))
end subroutine
```
@[1-5](module procedureを使用してREALの配列とINTEGERの配列ベクトルの配列でそれぞれ処理を行います。)
@[6-11](REALの配列を扱う関数内部では、配列をargsort関数を使用して並び替えています。)
@[12-16](INTEGERの配列を扱う関数でも処理は同様です。)

---

```fortran
subroutine sortVecs(vecs)
! sorts array of vectors, vecs, by length, from smallest to largest
real(dp), intent(inout):: vecs(:,:) ! array of vectors: vecs(i,j) = ith comp. of jth vec.
real(dp) len2(size(vecs,2))         ! array of squares of vector lengths
integer i
do i=1,size(len2)
   len2(i)=dot_product(vecs(:,i),vecs(:,i))
end do
call sortpairs(len2,vecs)
end subroutine
```
@[1-5](ベクトルをソートする関数はこちらです。)
@[6-9](ベクトル距離の2乗を計算し、sortpairs関数を使用してベクトルを並び変えます。)

---

```fortran
! overload argsort
interface argsort
    module procedure iargsort, rargsort
end interface
```
@[1-4](先程説明したargsortについても説明します。)
@[1-4](argsortもmodule procedureを使用してREALの配列とINTEGERの配列ベクトルの配列でそれぞれ処理を行います。)

---

```fortran
function iargsort(a) result(b)
a2 = a
N=size(a)
do i = 1, N
    b(i) = i
end do
do i = 1, N-1
    ! find ith smallest in 'a'
    imin = minloc(a2(i:),1) + i - 1

    ! swap to position i in 'a' and 'b', if not already there
    if (imin /= i) then
        temp = a2(i); a2(i) = a2(imin); a2(imin) = temp
        temp = b(i); b(i) = b(imin); b(imin) = temp
    end if
end do
end function
```
@[1-17](INTEGER配列のsort関数はこちらです。)
@[9](fortranのminloc関数を使って配列を並び替えています。)
@[1-17](ソートの方法については指定できません。)

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
@[1-6](次に行列式を計算する関数detについて解説します。)
@[8-10](実数と複素数の3×3の行列と行列式のための変数を定義します。)
@[12-14](実数型の計算と確認はこのように行います。)
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
! test a multiple of the unit matrix:
B = 3*eye(5)
call eig(B, lamb, cb)
call assert(maxval(abs(lamb - 3.0_dp)) < eps)  ! all eigenvalues are 3
call assert(maxval(abs(cb - cmplx(eye(5)))) < eps)  ! eigenvectors are cartesian unit basis vectors
```
@[1-5](単位行列に3を掛けた行列の固有値と固有ベクトルを計算します。)
@[1-5](固有値は3となります。)

---

```fortran
  ! eigenvalue/-vector problem for general matrices:
  interface eig
     module procedure deig
     module procedure zeig
  end interface eig

  ! TODO: add optional switch for left or right eigenvectors in deig() and zeig()?
  subroutine deig(A, lam, c)

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
```
@[1-5](module procedureを使用して実数の引数を持つ場合と複素数の引数を持つ場合の処理を分けています。)
@[7-13](まずはAとcが正方行列であるかを確認します。)
@[14-21](その後Lapackのdgeevを呼び出して固有値を計算します。)

---

```fortran
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
```
@[1-12](infoの値を確認することでDGEEVのエラー処理を行います。)
@[13](エラーがなければ先程計算したwrとwiから固有値を計算します。)

---

```fortran
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
```
@[1-13](DGEEVの固有ベクトルは実部と虚部で分けて得られます。)
@[5-11](与え方が特殊で固有値の虚部の符号に固有ベクトルの虚部の符号も合わせる必要があります。)

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
@[1-6](traceについて説明します。)
@[8-10](traceは正方行列に対して対角行列の和を計算したものです。)
@[12-16](他の場合と同様に確認を行います。この場合のtraceの値は15です。)


---

```fortran
  ! trace of real/complex matrices:
  interface trace
     module procedure dtrace
     module procedure ztrace
  end interface trace

  function dtrace(A) result(t)
    t = 0.0_dp
    do i = 1,minval(shape(A))
       t = t + A(i,i)
    end do
  end function dtrace
```
@[1-6](今回も実数と複素数は同じ形なので実数について説明します。)
@[7-12](traceについては対角部分を足し合わせるだけで実現ができます。)

---

#### まとめ
- sort関数と線形代数系の関数についてコードリーディングをしました。
- argsort関数はソートの種類を選ぶなどの改善ができるかもしれません。
- テンプレートなどを使えば実数と複素数の部分は一般化できるかもしれません。ただ、Lapackのルーチンが異なるので難しいかもしれません。
- TODOが結構あるのでそこを基に改良を加えられるかも。

