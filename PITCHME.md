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

[テスト結果](https://github.com/tkoyama010/fortran-utils/runs/232207124)

[コミット](https://github.com/tkoyama010/fortran-utils/commit/1995866a2b802476838dbee847a6fe4f7e60c249)

---

### 関数の未使用@size[2.0em](🏃) 

[テスト結果](https://github.com/tkoyama010/fortran-utils/runs/232207124)

[コミット](https://github.com/tkoyama010/fortran-utils/commit/5e3ac66a3be635819f4d1be6a84bcb878da3a4d8)

---

### 負になる可能性があるインデックス@size[2.0em](🏃) 

[テスト結果](https://travis-ci.com/tkoyama010/fortran-utils/builds/128721873)

[コミット](https://github.com/tkoyama010/fortran-utils/commit/c4c39f6a9fc64638f99f1c0f554ff0b3128bf9d9)

---

### 実数値が等しいかの比較@size[2.0em](🏃) 

[テスト結果](https://github.com/tkoyama010/fortran-utils/runs/236940764)

[コミット](https://github.com/tkoyama010/fortran-utils/commit/bcbdd5dc10da180db00c001a64215928006183bd)

---

### 倍精度実数からの単精度複素数の作成@size[2.0em](🏃) 

[テスト結果](https://github.com/tkoyama010/fortran-utils/runs/239593677)
[コミット](https://github.com/tkoyama010/fortran-utils/commit/1ea516cb35924167d16d9ecab72a9339ce9bd0ce)

---

### レビュー@size[2.0em](👨‍💻) 

[レビュー](https://github.com/certik/fortran-utils/pull/24)

---

### まとめ

- TravisのエラーをもとにPR
- 新しい[issue](https://github.com/certik/fortran-utils/issues/25)を作成

---

### Tips
https://twitter.com/OndrejCertik/status/1150507549822558208
