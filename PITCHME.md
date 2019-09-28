## fortran-utilsにPRを送ってみた@size[2.0em](🤠) 
[@tkoyama010](https://twitter.com/tkoyama010)

---


### fortran-utilsとは？@size[2.0em](🤠) 

[certik](https://github.com/certik)氏により2016年まで開発されていたfortranの便利なutils集です。

https://github.com/certik/fortran-utils

一部テストが通らない部分があったのでforkして変更をしました。

https://github.com/tkoyama010/fortran-utils

---


### fortran-utilsとは？@size[2.0em](🤠) 

[certik](https://github.com/certik)氏により2016年まで開発されていたfortranの便利なutils集です。

https://github.com/certik/fortran-utils

~~一部テストが通らない部分があったのでforkして変更をしました。~~

~~https://github.com/tkoyama010/fortran-utils~~

---

### 反省文@size[2.0em](🙇) 

けっこう適当にやってました@size[2.0em](💦)

@size[2.0em](👼)よい子なのでまじめに修正して[PR](https://github.com/certik/fortran-utils/pull/24)作りました

---

### masterのテスト状況@size[2.0em](🤔) 

[テスト結果(2019/9/29現在)](https://github.com/tkoyama010/fortran-utils/runs/232185941)

[テストに関して以前出されたIssue](https://github.com/certik/fortran-utils/issues/19)

WARNINGが1つでも出たらNGになるコンパイルオプションです@size[2.0em](😱) 

---

### masterのテスト状況@size[2.0em](🤔) 

[テスト結果(2019/9/29現在)](https://github.com/tkoyama010/fortran-utils/runs/232185941)

[テストに関して以前出されたIssue](https://github.com/certik/fortran-utils/issues/19)

このテストが通るまで直し続けます。ご指摘お願いします@size[2.0em](🙇) 

---

### DO文の修正@size[2.0em](🏃) 

[コミット](https://github.com/certik/fortran-utils/pull/24/commits/1995866a2b802476838dbee847a6fe4f7e60c249)

[テスト結果](https://github.com/tkoyama010/fortran-utils/runs/232200514)

---

### DO文の修正@size[2.0em](🏃) 

[コミット](https://github.com/tkoyama010/fortran-utils/commit/7f4c355f8867da451be6192d50896eff95187035)

[テスト結果](https://github.com/tkoyama010/fortran-utils/runs/232207124)

---

### 関数の未使用@size[2.0em](🏃) 

[コミット](https://github.com/tkoyama010/fortran-utils/commit/5e3ac66a3be635819f4d1be6a84bcb878da3a4d8)

[テスト結果](https://github.com/tkoyama010/fortran-utils/runs/232213983)

---

### 実数値の比較のエラー@size[2.0em](🏃) 

[コミット](https://github.com/tkoyama010/fortran-utils/commit/95aef3deae9c1a61e33c40dd4ebb7823430e74f4)

[テスト結果](https://travis-ci.com/tkoyama010/fortran-utils/builds/128705063)

❌この修正は誤りです(後ででてきます)

---

### 実数値の比較のエラー@size[2.0em](🏃) 

[コミット](https://github.com/tkoyama010/fortran-utils/commit/fe79d568bdc217a04e114f4b6f43e64fe93d2c99)

[テスト結果](https://github.com/tkoyama010/fortran-utils/runs/232289476)

❌この修正は誤りです(後ででてきます)

[テスト結果](https://travis-ci.com/tkoyama010/fortran-utils/builds/128721873)

---

### 実数値の比較のエラー@size[2.0em](🏃) 

[コミット](https://github.com/tkoyama010/fortran-utils/commit/c4c39f6a9fc64638f99f1c0f554ff0b3128bf9d9)
[テスト結果](https://github.com/tkoyama010/fortran-utils/runs/232416922)
