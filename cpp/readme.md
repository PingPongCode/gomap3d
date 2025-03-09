# gomap3d Cgo C++支持

## windows下使用

+ 编译动态库 go build -o main.dll -buildmode=c-shared main.go
+ 交叉编译 g++ main.cpp main.dll
