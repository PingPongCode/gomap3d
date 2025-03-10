# gomap3d Cgo C++支持

## 简介

main.cpp为windows平台下调用示例

mainLinux.cpp为linux平台下调用示例

coordinateSystem可以填写wgs84以及cgcs2000

## windows下使用

+ 编译动态库 go build -o lib/main.dll -buildmode=c-shared main.go
+ 交叉编译 g++ main.cpp lib/main.dll -o main.exe
+ 运行 main.exe

## Linux下使用

+ 编译动态库 go build -o lib/libmain.so -buildmode=c-shared main.go
+ 交叉编译 g++ mainLinux.cpp lib/libmain.so -o mainLinux.out
+ 运行 ./mainLinux.out
