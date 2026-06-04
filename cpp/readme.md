# gomap3d C/C++ 支持

本目录提供两种在 C/C++ 中调用 gomap3d 速度/位置转换函数的方式：

- **CGo 动态链接库** — 编译为 DLL (Windows) 或 SO (Linux)
- **纯 C 头文件库** — 零依赖，直接 `#include`

## 方式一：CGo 动态链接库

### Windows

```cmd
rem 编译 DLL（在 cpp/ 目录下执行）
set CGO_ENABLED=1
go build -buildmode=c-shared -o gomap3d.dll .

rem 编译 C++ 测试程序
g++ -std=c++11 -I. -o tests/test_gomap3d.exe tests/main.cpp gomap3d.dll

rem 运行
tests\test_gomap3d.exe
```

### Linux

```bash
# 编译 SO
CGO_ENABLED=1 go build -buildmode=c-shared -o libgomap3d.so .

# 编译 C++ 测试程序
g++ -std=c++11 -I. -o tests/test_gomap3d tests/main.cpp ./libgomap3d.so

# 运行
LD_LIBRARY_PATH=. ./tests/test_gomap3d
```

生成的 `gomap3d.h` 包含所有导出函数的声明，C/C++ 中 `#include "gomap3d.h"` 即可使用。

## 方式二：纯 C 头文件（推荐）

`velocity_c.h` 是纯 C 语言实现的**头文件库**，所有函数均为 `static inline`，无需链接任何外部库。

将以下代码保存为 `test.c`：

```c
#include "velocity_c.h"

int main() {
    double vx, vy, vz;
    ecefVel2eciVel_c(1500, 2500, 800,     // ECEF速度
                     5000000, 3000000, 2000000,  // ECEF位置
                     1780472536.0,         // Unix时间戳
                     &vx, &vy, &vz);
    return 0;
}
```

编译：`gcc -std=c99 -I.. -o test_c test.c -lm`（Windows CMD 中文乱码可加 `-finput-charset=utf-8 -fexec-charset=gbk`）

> 注意：使用前需将上方示例代码保存为 `test.c`

## 速度转换函数清单

| C 函数（velocity_c.h）  | CGo 导出（gomap3d.dll） | 说明                             |
| ----------------------- | ----------------------- | -------------------------------- |
| `ecefVel2eciVel_c()`  | `ECEFVel2ECIVel`      | ECEF 速度 → ECI 速度            |
| `eciVel2ecefVel_c()`  | `ECIVel2ECEFVel`      | ECI 速度 → ECEF 速度            |
| `enuVel2ecefVel_c()`  | `ENUVel2ECEFVel`      | ENU 速度 → ECEF 速度            |
| `ecefVel2enuVel_c()`  | `ECEFVel2ENUVel`      | ECEF 速度 → ENU 速度            |
| `aerDeriv2enuVel_c()` | `AERDeriv2ENUVel`     | AER 变化率 → ENU 速度           |
| `enuVel2aerDeriv_c()` | `ENUVel2AERDeriv`     | ENU 速度 → AER 变化率           |
| `aerDeriv2eciVel_c()` | `AERDeriv2ECIVel`     | AER 变化率 → ECI 速度（一站式） |
| `geodetic2ecef_c()`   | -                       | 经纬高 → ECEF（纯 C 辅助函数）  |
| `ecef2enu_c()`        | -                       | ECEF → ENU（纯 C 辅助函数）     |
| `enu2ecef_c()`        | -                       | ENU → ECEF（纯 C 辅助函数）     |
| `enu2aer_c()`         | -                       | ENU → AER（纯 C 辅助函数）      |
| `juliandate_c()`      | `juliandate`          | Unix 时间戳 → 儒略日            |
| `greenwichsrt_c()`    | `greenwichsrt`        | 儒略日 → 格林威治恒星时         |

## 测试

测试文件位于 `tests/` 目录：

```bash
cd tests

# C++ 测试（需先编译 DLL）
g++ -std=c++11 -I.. -o test_gomap3d.exe main.cpp ../gomap3d.dll
./test_gomap3d.exe

# 纯 C 测试（无需 DLL）
gcc -std=c99 -I.. -o test_velocity_c.exe test_velocity_c.c -lm
./test_velocity_c.exe
```

## 注意事项

- `mainLinux.cpp` 为 Linux 平台调用示例，位于 `linux/` 目录
- 椭球体参数 `coordinateSystem` 可填写 `wgs84` 或 `cgcs2000`
