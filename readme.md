# gomap3d

Go 语言实现的多坐标系转换库，支持天文学/航天领域常用坐标系转换与速度矢量转换。

通过 CGo 提供 C/C++ 调用支持，详见 `cpp/` 文件夹。

## 特性

- **位置坐标转换**：5 种坐标系互转
  - 站心坐标系 (AER)
  - 东北天坐标系 (ENU)
  - 地心地固坐标系 (ECEF)
  - 地心惯性坐标系 (ECI)
  - 大地坐标系 (LLA/Geodetic)

- **速度矢量转换**：6 种速度转换（含正逆方向）
  - ECEF ↔ ECI 速度（含地球自转修正）
  - ENU ↔ ECEF 速度
  - AER 变化率（RAE 导数）↔ ENU 速度
  - 一站式：AER 变化率 ↔ ECI 速度

- **多种参考椭球体**
  - WGS-84
  - CGCS2000
  - 月球
  - 火星

- **精确天文计算**
  - 儒略日计算
  - 格林威治恒星时
  - ECI/ECEF 时变转换

- **C/C++ 支持**
  - CGo 动态链接库 (DLL/SO)
  - 纯 C 头文件库（零依赖，直接 `#include`）

## 安装

```bash
go get github.com/PingPongCode/gomap3d
```

## 使用示例

### 单位

角度相关字段（经度、纬度、方位角、俯仰角）单位均为 **度 (°)**，长度单位均为 **米 (m)**。

### 基本坐标转换

```go
package main

import (
	"fmt"
	"time"
	"github.com/PingPongCode/gomap3d"
)

func main() {
	// 创建WGS84椭球体 (也可使用cgcs2000)
	ell, _ := gomap3d.NewEllipsoid("wgs84")

	// 大地坐标 (北京)
	beijing := gomap3d.Geodetic{
		Latitude:  39.9042,
		Longitude: 116.4074,
		Altitude:  43.5,
		Ell:       ell,
	}

	// 转换为ECEF
	ecef := beijing.ToECEF()
	fmt.Printf("ECEF坐标: %.2f, %.2f, %.2f\n", ecef.X, ecef.Y, ecef.Z)

	// 转换为ENU (以上海为参考点)
	shanghai := gomap3d.Geodetic{
		Latitude:  31.2304,
		Longitude: 121.4737,
		Altitude:  4.0,
		Ell:       ell,
	}
	enu := beijing.ToENU(shanghai)
	fmt.Printf("ENU坐标: 东%.2fm, 北%.2fm, 上%.2fm\n", enu.East, enu.North, enu.Up)

	// 时间相关转换 (ECI)
	t := time.Date(2023, 6, 15, 12, 0, 0, 0, time.UTC)
	eci := ecef.ToECI(t)
	fmt.Printf("ECI坐标: %.2f, %.2f, %.2f\n", eci.X, eci.Y, eci.Z)
}
```

### 速度矢量转换

```go
package main

import (
	"fmt"
	"time"
	"github.com/PingPongCode/gomap3d"
)

func main() {
	t := time.Date(2026, 6, 3, 7, 42, 16, 0, time.UTC)

	// ECEF速度 → ECI速度
	vx, vy, vz := gomap3d.ECEFVel2ECIVel(
		1500, 2500, 800,   // ECEF速度 (m/s)
		5000000, 3000000, 2000000, // ECEF位置 (m)
		t,
	)
	fmt.Printf("ECI速度: %.2f, %.2f, %.2f m/s\n", vx, vy, vz)

	// ENU速度 → ECEF速度
	vx2, vy2, vz2 := gomap3d.ENUVel2ECEFVel(
		500, -300, 200,  // ENU速度 (m/s)
		40.0, 125.0,     // 测站经纬度
	)
	fmt.Printf("ECEF速度: %.2f, %.2f, %.2f m/s\n", vx2, vy2, vz2)

	// AER变化率 (RAE导数) → ECI速度 (一站式)
	vx3, vy3, vz3 := gomap3d.AERDeriv2ECIVel(
		1200000, 180, 45,  // 斜距(m), 方位角(°), 俯仰角(°)
		3000, 0, 0.01,     // 斜距变化率, 方位角变化率, 俯仰角变化率
		40.0, 116.4, t,    // 测站经纬度, 时间
	)
	fmt.Printf("ECI速度: %.2f, %.2f, %.2f m/s\n", vx3, vy3, vz3)
}
```

## 函数清单

### 位置坐标转换 (base.go)

```go
func ENU2AER(e, n, u float64) (az, el, srange float64)
func AER2ENU(az, el, srange float64) (e, n, u float64)
func Geodetic2ECEF(lat, lon, alt float64, ell *Ellipsoid) (x, y, z float64)
func ECEF2Geodetic(x, y, z float64, ell *Ellipsoid) (lat, lon, alt float64)
func ECEF2ENU(x, y, z, lat0, lon0, h0 float64, ell *Ellipsoid) (e, n, u float64)
func ENU2ECEF(e, n, u, lat0, lon0, h0 float64, ell *Ellipsoid) (x, y, z float64)
func ECI2ECEF(x, y, z float64, t time.Time) (xEcef, yEcef, zEcef float64)
func ECEF2ECI(x, y, z float64, t time.Time) (xEci, yEci, zEci float64)
```

### 位置转换类型方法

每种坐标类型均提供 `ToXXX()` 方法，支持链式调用：

```go
aer.ToENU()
enu.ToECEF(ref Geodetic)
ecef.ToECI(t time.Time)
eci.ToGeodetic()
geodetic.ToAER(ref Geodetic)
// ... 详见各 .go 文件
```

### 速度矢量转换 (velocity.go)

```go
// ECEF速度 ↔ ECI速度
func ECEFVel2ECIVel(vx, vy, vz, x, y, z float64, t time.Time) (vxEci, vyEci, vzEci float64)
func ECIVel2ECEFVel(vx, vy, vz, x, y, z float64, t time.Time) (vxEcef, vyEcef, vzEcef float64)

// ENU速度 ↔ ECEF速度
func ENUVel2ECEFVel(eVel, nVel, uVel, latDeg, lonDeg float64) (vx, vy, vz float64)
func ECEFVel2ENUVel(vx, vy, vz, latDeg, lonDeg float64) (e, n, u float64)

// AER变化率 ↔ ENU速度
func AERDeriv2ENUVel(R, azDeg, elDeg, dR, dAzDeg, dElDeg float64) (e, n, u float64)
func ENUVel2AERDeriv(eVel, nVel, uVel, R, azDeg, elDeg float64) (dR, dAzDeg, dElDeg float64)

// 一站式：AER变化率 ↔ ECI速度
func AERDeriv2ECIVel(R, azDeg, elDeg, dR, dAzDeg, dElDeg, latDeg, lonDeg float64, t time.Time) (vx, vy, vz float64)
func ECIVel2AERDeriv(vx, vy, vz, Rx, Ry, Rz, latDeg, lonDeg float64, t time.Time) (dR, dAzDeg, dElDeg float64)
```

### 天文计算 (base.go)

```go
func juliandate(t time.Time) float64
func greenwichsrt(jd float64) float64
```

## C/C++ 支持

本库支持两种方式在 C/C++ 代码中使用：

### 方式一：CGo 动态链接库

通过 CGo 编译为 DLL (Windows) 或 SO (Linux) 供 C/C++ 程序调用。

**Windows（CMD）：**
```cmd
cd cpp
set CGO_ENABLED=1
go build -buildmode=c-shared -o gomap3d.dll .
```

**Linux（bash）：**
```bash
cd cpp
CGO_ENABLED=1 go build -buildmode=c-shared -o libgomap3d.so .
```

生成的 `gomap3d.h` 头文件可直接在 C/C++ 中 `#include` 使用。

### 方式二：纯 C 头文件 (推荐)

`cpp/velocity_c.h` 是纯 C 语言的函数实现，**零外部依赖**，直接 `#include` 即可使用。所有函数均为 `static inline`，编译后无额外链接开销。

将以下代码保存为 `test.c` 并编译：

```c
#include "velocity_c.h"

int main() {
    double vx, vy, vz;
    ecefVel2eciVel_c(1500, 2500, 800,
                     5000000, 3000000, 2000000,
                     1780472536.0,  // Unix时间戳
                     &vx, &vy, &vz);
    return 0;
}
```

编译：`gcc -std=c99 -I cpp -o test test.c -lm`（Windows CMD 中文乱码可加 `-finput-charset=utf-8 -fexec-charset=gbk`）

C 头文件支持的函数清单：

| C 函数 | 对应 Go 函数 | 说明 |
|--------|-------------|------|
| `ecefVel2eciVel_c()` | `ECEFVel2ECIVel` | ECEF速度 → ECI速度 |
| `eciVel2ecefVel_c()` | `ECIVel2ECEFVel` | ECI速度 → ECEF速度 |
| `enuVel2ecefVel_c()` | `ENUVel2ECEFVel` | ENU速度 → ECEF速度 |
| `ecefVel2enuVel_c()` | `ECEFVel2ENUVel` | ECEF速度 → ENU速度 |
| `aerDeriv2enuVel_c()` | `AERDeriv2ENUVel` | AER变化率 → ENU速度 |
| `enuVel2aerDeriv_c()` | `ENUVel2AERDeriv` | ENU速度 → AER变化率 |
| `aerDeriv2eciVel_c()` | `AERDeriv2ECIVel` | AER变化率 → ECI速度 |
| `juliandate_c()` | `juliandate` | Unix时间戳 → 儒略日 |
| `greenwichsrt_c()` | `greenwichsrt` | 儒略日 → 格林威治恒星时 |

## 测试

所有测试用例覆盖位置坐标转换与速度矢量转换的正向/逆向可逆性验证。

```bash
# 运行全部测试
go test .

# 仅运行速度转换测试
go test -run "TestECEF|TestENU|TestAER|TestFull|TestKnown" -v .

# C 头文件测试
cd cpp/tests
gcc -std=c99 -I.. -o test_velocity_c test_velocity_c.c -lm && ./test_velocity_c
```

## 贡献

欢迎提交 Issue 和 PR。提交代码前请确保：

1. 通过所有测试 `go test .`
2. 添加新功能的测试用例
3. 更新相关文档

## 许可证

MIT License
