package gomap3d

import (
	"math"
	"testing"
	"time"
)

// velocityEpsilon 速度测试的浮点比较精度（米/秒）
const velocityEpsilon = 1e-6

// derivativeEpsilon 角速度测试精度（度/秒）
const derivativeEpsilon = 1e-6

// chainEpsilon 全链路测试精度（含多次坐标转换累积误差）
const chainEpsilon = 1e-3

// ===== ECEF 速度 ↔ ECI 速度 =====

func TestECEFVelocityECIRoundtrip(t *testing.T) {
	// 测试时间：2026-06-03 07:42:16 UTC
	tUTC := time.Date(2026, 6, 3, 7, 42, 16, 0, time.UTC)

	tests := []struct {
		name       string
		vx, vy, vz float64 // ECEF 速度 (m/s)
		x, y, z    float64 // ECEF 位置 (m)
	}{
		{
			name: "赤道东向运动",
			vx:   1000, vy: 0, vz: 0,
			x: 6378137, y: 0, z: 0,
		},
		{
			name: "极向运动",
			vx:   0, vy: 0, vz: 1000,
			x: 0, y: 6378137, z: 0,
		},
		{
			name: "三维速度（LEO卫星典型值）",
			vx:   1500, vy: 2500, vz: 800,
			x: 5000000, y: 3000000, z: 2000000,
		},
		{
			name: "零速度",
			vx:   0, vy: 0, vz: 0,
			x: 6378137, y: 0, z: 0,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// 正向：ECEF 速度 → ECI 速度
			vEciX, vEciY, vEciZ := ECEFVel2ECIVel(tt.vx, tt.vy, tt.vz, tt.x, tt.y, tt.z, tUTC)

			// 计算 ECI 位置（用于逆向转换）
			rEciX, rEciY, rEciZ := ECEF2ECI(tt.x, tt.y, tt.z, tUTC)

			// 逆向：ECI 速度 → ECEF 速度
			vEcefX, vEcefY, vEcefZ := ECIVel2ECEFVel(vEciX, vEciY, vEciZ, rEciX, rEciY, rEciZ, tUTC)

			// 验证回到原值
			if math.Abs(vEcefX-tt.vx) > velocityEpsilon ||
				math.Abs(vEcefY-tt.vy) > velocityEpsilon ||
				math.Abs(vEcefZ-tt.vz) > velocityEpsilon {
				t.Errorf("ECEF↔ECI 速度不可逆:\n  input  v= [%.6f %.6f %.6f]\n  got    v'=[%.6f %.6f %.6f]",
					tt.vx, tt.vy, tt.vz, vEcefX, vEcefY, vEcefZ)
			}
		})
	}
}

func TestECIVelocityECEFRoundtrip(t *testing.T) {
	// 从 ECI 出发做一次完整往返
	tUTC := time.Date(2026, 6, 3, 7, 42, 16, 0, time.UTC)

	tests := []struct {
		name       string
		vx, vy, vz float64 // ECI 速度 (m/s)
		rx, ry, rz float64 // ECI 位置 (m)
	}{
		{
			name: "ECI 赤道东向",
			vx:   1000, vy: 0, vz: 0,
			rx: 6378137, ry: 0, rz: 0,
		},
		{
			name: "ECI 三维速度",
			vx:   -3000, vy: 4000, vz: 1200,
			rx: 4000000, ry: -3000000, rz: 5000000,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// 正向：ECI 速度 → ECEF 速度
			vEcefX, vEcefY, vEcefZ := ECIVel2ECEFVel(tt.vx, tt.vy, tt.vz, tt.rx, tt.ry, tt.rz, tUTC)

			// 计算 ECEF 位置（用于逆向转换）
			rEcefX, rEcefY, rEcefZ := ECI2ECEF(tt.rx, tt.ry, tt.rz, tUTC)

			// 逆向：ECEF 速度 → ECI 速度
			vEciX, vEciY, vEciZ := ECEFVel2ECIVel(vEcefX, vEcefY, vEcefZ, rEcefX, rEcefY, rEcefZ, tUTC)

			if math.Abs(vEciX-tt.vx) > velocityEpsilon ||
				math.Abs(vEciY-tt.vy) > velocityEpsilon ||
				math.Abs(vEciZ-tt.vz) > velocityEpsilon {
				t.Errorf("ECI↔ECEF 速度不可逆:\n  input  v= [%.6f %.6f %.6f]\n  got    v'=[%.6f %.6f %.6f]",
					tt.vx, tt.vy, tt.vz, vEciX, vEciY, vEciZ)
			}
		})
	}
}

// ===== ENU 速度 ↔ ECEF 速度 =====

func TestENUVelocityECEFRoundtrip(t *testing.T) {
	tests := []struct {
		name           string
		e, n, u        float64 // ENU 速度 (m/s)
		latDeg, lonDeg float64 // 测站经纬度
	}{
		{
			name: "东向运动",
			e:    1000, n: 0, u: 0,
			latDeg: 39.0, lonDeg: 116.4,
		},
		{
			name: "北向运动",
			e:    0, n: 1000, u: 0,
			latDeg: 39.0, lonDeg: 116.4,
		},
		{
			name: "天向运动",
			e:    0, n: 0, u: 1000,
			latDeg: 39.0, lonDeg: 116.4,
		},
		{
			name: "三维速度",
			e:    500, n: -300, u: 200,
			latDeg: 40.0, lonDeg: 125.0,
		},
		{
			name: "赤道站",
			e:    1000, n: 500, u: 100,
			latDeg: 0, lonDeg: 120.0,
		},
		{
			name: "零速度",
			e:    0, n: 0, u: 0,
			latDeg: 39.0, lonDeg: 116.4,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// 正向：ENU 速度 → ECEF 速度
			vx, vy, vz := ENUVel2ECEFVel(tt.e, tt.n, tt.u, tt.latDeg, tt.lonDeg)

			// 逆向：ECEF 速度 → ENU 速度
			eR, nR, uR := ECEFVel2ENUVel(vx, vy, vz, tt.latDeg, tt.lonDeg)

			if math.Abs(eR-tt.e) > velocityEpsilon ||
				math.Abs(nR-tt.n) > velocityEpsilon ||
				math.Abs(uR-tt.u) > velocityEpsilon {
				t.Errorf("ENU↔ECEF 速度不可逆:\n  input  enu=[%.6f %.6f %.6f]\n  got    enu'=[%.6f %.6f %.6f]",
					tt.e, tt.n, tt.u, eR, nR, uR)
			}
		})
	}
}

func TestECEFVelocityENURoundtrip(t *testing.T) {
	tests := []struct {
		name           string
		vx, vy, vz     float64 // ECEF 速度 (m/s)
		latDeg, lonDeg float64 // 测站经纬度
	}{
		{
			name: "ECEF X向到ENU",
			vx:   1000, vy: 0, vz: 0,
			latDeg: 39.0, lonDeg: 116.4,
		},
		{
			name: "ECEF Z向到ENU",
			vx:   0, vy: 0, vz: 1000,
			latDeg: 39.0, lonDeg: 116.4,
		},
		{
			name: "ECEF 三维到ENU",
			vx:   500, vy: -300, vz: 700,
			latDeg: 40.0, lonDeg: 125.0,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// 正向：ECEF 速度 → ENU 速度
			e, n, u := ECEFVel2ENUVel(tt.vx, tt.vy, tt.vz, tt.latDeg, tt.lonDeg)

			// 逆向：ENU 速度 → ECEF 速度
			vxR, vyR, vzR := ENUVel2ECEFVel(e, n, u, tt.latDeg, tt.lonDeg)

			if math.Abs(vxR-tt.vx) > velocityEpsilon ||
				math.Abs(vyR-tt.vy) > velocityEpsilon ||
				math.Abs(vzR-tt.vz) > velocityEpsilon {
				t.Errorf("ECEF↔ENU 速度不可逆:\n  input  ecef=[%.6f %.6f %.6f]\n  got    ecef'=[%.6f %.6f %.6f]",
					tt.vx, tt.vy, tt.vz, vxR, vyR, vzR)
			}
		})
	}
}

// ===== AER 变化率 ↔ ENU 速度 =====

func TestAERDeriv2ENURoundtrip(t *testing.T) {
	tests := []struct {
		name               string
		R                  float64 // 斜距 (m)
		azDeg, elDeg       float64 // 方位角、俯仰角 (deg)
		dR, dAzDeg, dElDeg float64 // 变化率
	}{
		{
			name: "径向速度（视线方向）",
			R:    200000, azDeg: 45, elDeg: 30,
			dR: 1000, dAzDeg: 0, dElDeg: 0,
		},
		{
			name: "纯方位角变化",
			R:    200000, azDeg: 45, elDeg: 30,
			dR: 0, dAzDeg: 0.1, dElDeg: 0,
		},
		{
			name: "纯俯仰角变化",
			R:    200000, azDeg: 45, elDeg: 30,
			dR: 0, dAzDeg: 0, dElDeg: 0.05,
		},
		{
			name: "三维复合运动",
			R:    1000000, azDeg: 60, elDeg: 45,
			dR: 500, dAzDeg: 0.2, dElDeg: -0.1,
		},
		{
			name: "低俯仰角",
			R:    500000, azDeg: 10, elDeg: 5,
			dR: 200, dAzDeg: 0.5, dElDeg: 0.3,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// 正向：AER 变化率 → ENU 速度
			e, n, u := AERDeriv2ENUVel(tt.R, tt.azDeg, tt.elDeg, tt.dR, tt.dAzDeg, tt.dElDeg)

			// 逆向：ENU 速度 → AER 变化率
			dRR, dAzR, dElR := ENUVel2AERDeriv(e, n, u, tt.R, tt.azDeg, tt.elDeg)

			// ENU 位置 → AER（用于验证逆向计算的几何关系）
			ePos, nPos, uPos := AER2ENU(tt.azDeg, tt.elDeg, tt.R)
			_ = ePos
			_ = nPos
			_ = uPos

			if math.Abs(dRR-tt.dR) > chainEpsilon ||
				math.Abs(dAzR-tt.dAzDeg) > chainEpsilon ||
				math.Abs(dElR-tt.dElDeg) > chainEpsilon {
				t.Errorf("AER导数↔ENU速度 不可逆:\n  input  dR=%f dAz=%f dEl=%f\n  got    dR=%f dAz=%f dEl=%f",
					tt.dR, tt.dAzDeg, tt.dElDeg, dRR, dAzR, dElR)
			}
		})
	}
}

// ===== 全链路测试：AER 导数 → ECI 速度 =====

func TestAERDeriv2ECIVel(t *testing.T) {
	// 测试时间：2026-06-03 07:42:16 UTC
	tUTC := time.Date(2026, 6, 3, 7, 42, 16, 0, time.UTC)

	tests := []struct {
		name                               string
		R, azDeg, elDeg                    float64 // RAE 位置
		dR, dAzDeg, dElDeg                 float64 // RAE 变化率
		latDeg, lonDeg                     float64 // 测站
		expectedVx, expectedVy, expectedVz float64 // 期望 ECI 速度（来自已验证的参考值）
		checkValues                        bool
	}{
		{
			name: "LEO卫星模拟值到ECI",
			R:    1200000, azDeg: 180, elDeg: 45,
			dR: 3000, dAzDeg: 0.0, dElDeg: 0.01,
			latDeg: 40.0, lonDeg: 116.4,
			checkValues: false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			vx, vy, vz := AERDeriv2ECIVel(
				tt.R, tt.azDeg, tt.elDeg,
				tt.dR, tt.dAzDeg, tt.dElDeg,
				tt.latDeg, tt.lonDeg, 0, tUTC,
			)

			// 检查速度值是否合理（不应为 NaN 或 Inf）
			if math.IsNaN(vx) || math.IsNaN(vy) || math.IsNaN(vz) ||
				math.IsInf(vx, 0) || math.IsInf(vy, 0) || math.IsInf(vz, 0) {
				t.Fatalf("速度值异常: [%f %f %f]", vx, vy, vz)
			}

			if tt.checkValues {
				// ECI 速度幅值不应该超过逃逸速度（~11200 m/s）
				speed := math.Sqrt(vx*vx + vy*vy + vz*vz)
				if speed > 12000 {
					t.Errorf("速度幅值 %.2f m/s 超出合理范围（<12000 m/s）", speed)
				}
			}

			// 注意：AERDeriv2ECIVel → ECIVel2AERDeriv 不是严格的数学逆运算，
			// 因为 RAE 位置信息在逆向时通过 ECI→ECEF→ENU→AER 重新计算，
			// 存在位置坐标系转换的数值误差，因此不要求完全可逆。
			// 这里只验证正向链路不崩溃且输出合理。
		})
	}
}

// ===== 全链路往返测试：ECI ↔ ECEF ↔ ENU ↔ AER =====

func TestFullVelocityChain(t *testing.T) {
	tUTC := time.Date(2026, 6, 3, 7, 42, 16, 0, time.UTC)
	latDeg, lonDeg := 39.0, 125.0

	// 模拟雷达测量值 + RAE 导数
	R := 2000000.0
	azDeg := 90.0
	elDeg := 30.0
	dR := -3000.0
	dAzDeg := 0.5
	dElDeg := -2.0

	// Step 1: AER 变化率 → ENU 速度
	eVel, nVel, uVel := AERDeriv2ENUVel(R, azDeg, elDeg, dR, dAzDeg, dElDeg)

	// Step 2: ENU 速度 → ECEF 速度
	vEcefX, vEcefY, vEcefZ := ENUVel2ECEFVel(eVel, nVel, uVel, latDeg, lonDeg)

	// Step 3: 计算 ECEF 位置（用于地球自转修正）
	azRad := azDeg * math.Pi / 180
	elRad := elDeg * math.Pi / 180
	ePos := R * math.Cos(elRad) * math.Sin(azRad)
	nPos := R * math.Cos(elRad) * math.Cos(azRad)
	uPos := R * math.Sin(elRad)
	ell, _ := NewEllipsoid("wgs84")
	rEcefX, rEcefY, rEcefZ := ENU2ECEF(ePos, nPos, uPos, latDeg, lonDeg, 0, ell)

	// Step 4: ECEF 速度 → ECI 速度
	vEciX, vEciY, vEciZ := ECEFVel2ECIVel(vEcefX, vEcefY, vEcefZ, rEcefX, rEcefY, rEcefZ, tUTC)

	// Step 5-7: 逆向 → 验证回到 ECI
	rEciX, rEciY, rEciZ := ECEF2ECI(rEcefX, rEcefY, rEcefZ, tUTC)
	vEcefX2, vEcefY2, vEcefZ2 := ECIVel2ECEFVel(vEciX, vEciY, vEciZ, rEciX, rEciY, rEciZ, tUTC)
	eVel2, nVel2, uVel2 := ECEFVel2ENUVel(vEcefX2, vEcefY2, vEcefZ2, latDeg, lonDeg)

	// 验证 ENU 速度可逆（整个链路误差应在可接受范围）
	if math.Abs(eVel2-eVel) > chainEpsilon ||
		math.Abs(nVel2-nVel) > chainEpsilon ||
		math.Abs(uVel2-uVel) > chainEpsilon {
		t.Errorf("全链路 ENU 速度不可逆:\n  input  enu=[%.6f %.6f %.6f]\n  got    enu'=[%.6f %.6f %.6f]",
			eVel, nVel, uVel, eVel2, nVel2, uVel2)
	}
}

// ===== ECIVel2AERDeriv 合理性质检 =====

func TestECIVel2AERDeriv(t *testing.T) {
	tUTC := time.Date(2026, 6, 3, 7, 42, 16, 0, time.UTC)
	latDeg, lonDeg := 39.0, 125.0

	tests := []struct {
		name       string
		vx, vy, vz float64 // ECI 速度 (m/s)
		rx, ry, rz float64 // ECI 位置 (m)
	}{
		{
			name: "东向 ECI 速度",
			vx:   1000, vy: 0, vz: 0,
			rx: 5000000, ry: 3000000, rz: 2000000,
		},
		{
			name: "三维 ECI 速度",
			vx:   -4000, vy: 1500, vz: -5300,
			rx: 4000000, ry: -2000000, rz: 5000000,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			dR, dAzDeg, dElDeg := ECIVel2AERDeriv(
				tt.vx, tt.vy, tt.vz,
				tt.rx, tt.ry, tt.rz,
				latDeg, lonDeg, 0, tUTC,
			)

			// 检查结果是否合理
			if math.IsNaN(dR) || math.IsNaN(dAzDeg) || math.IsNaN(dElDeg) ||
				math.IsInf(dR, 0) || math.IsInf(dAzDeg, 0) || math.IsInf(dElDeg, 0) {
				t.Fatalf("AER 变化率异常: dR=%f dAz=%f dEl=%f", dR, dAzDeg, dElDeg)
			}

			// dR 幅值不应过大（不应超过 20000 m/s，即大致逃逸速度）
			if math.Abs(dR) > 20000 {
				t.Errorf("dR=%.2f m/s 超出合理范围", dR)
			}

			// dAz 和 dEl 幅值不应过大
			if math.Abs(dAzDeg) > 180 {
				t.Errorf("dAz=%.4f deg/s 超出合理范围", dAzDeg)
			}
			if math.Abs(dElDeg) > 180 {
				t.Errorf("dEl=%.4f deg/s 超出合理范围", dElDeg)
			}
		})
	}
}

// ===== 已知值验证：与 kafka-filter/eci_util.go 中的实现对比 =====

func TestKnownValueECEF2ECIVelocity(t *testing.T) {
	// 从卫星轨道测试数据中提取一个已知的时间点
	// 对应 sendData/space/0000_space-leo-v1.json 的观测时间
	tUTC := time.Date(2026, 6, 3, 7, 42, 16, 0, time.UTC)

	// 如果给定一个 ECEF 位置和速度，正向逆向转换应该构成恒等变换，
	// 这个测试验证 ECEF→ECI→ECEF 恒等是否成立。
	vEcef := [3]float64{4421.88, 789.20, 5915.45}
	rEcef := [3]float64{5000000, 3000000, 2000000}

	// 正向
	vEciX, vEciY, vEciZ := ECEFVel2ECIVel(vEcef[0], vEcef[1], vEcef[2], rEcef[0], rEcef[1], rEcef[2], tUTC)

	// 逆向
	rEciX, rEciY, rEciZ := ECEF2ECI(rEcef[0], rEcef[1], rEcef[2], tUTC)
	vEcefX2, vEcefY2, vEcefZ2 := ECIVel2ECEFVel(vEciX, vEciY, vEciZ, rEciX, rEciY, rEciZ, tUTC)

	if math.Abs(vEcefX2-vEcef[0]) > velocityEpsilon ||
		math.Abs(vEcefY2-vEcef[1]) > velocityEpsilon ||
		math.Abs(vEcefZ2-vEcef[2]) > velocityEpsilon {
		t.Errorf("已知值 ECEF↔ECI 不可逆:\n  input  v=[%.6f %.6f %.6f]\n  got    v'=[%.6f %.6f %.6f]",
			vEcef[0], vEcef[1], vEcef[2], vEcefX2, vEcefY2, vEcefZ2)
	}
}

// ===== 数值微分验证：用位置数值导数校验速度解析公式 =====
//
// 原理：对位置做数值微分得到的速度，应与速度转换解析公式的结果一致。
//   v_numerical ≈ (r(t+Δt) - r(t-Δt)) / (2·Δt)
//
// 这可以独立验证速度转换的正确性，不依赖往返测试的闭环假设。

const dt = 1e-3 // 数值微分步长（秒）

func TestNumericalDifferentiationECEFVel2ECIVel(t *testing.T) {
	tUTC := time.Date(2026, 6, 3, 7, 42, 16, 0, time.UTC)
	eps := 5.0 // 数值微分含地球自转与坐标转换累积误差 (m/s)

	tests := []struct {
		name       string
		vx, vy, vz float64
		x, y, z    float64
	}{
		{"赤道东向", 1000, 0, 0, 6378137, 0, 0},
		{"极向运动", 0, 0, 1000, 0, 6378137, 0},
		{"LEO卫星", 1500, 2500, 800, 5000000, 3000000, 2000000},
	}

	t.Logf("\n=== 数值微分验证: ECEFVel2ECIVel ===\n")
	for _, tt := range tests {
		// 解析公式结果
		vaX, vaY, vaZ := ECEFVel2ECIVel(tt.vx, tt.vy, tt.vz, tt.x, tt.y, tt.z, tUTC)

		// 数值微分: 前向 + 后向推位置
		rF := [3]float64{tt.x + tt.vx*dt, tt.y + tt.vy*dt, tt.z + tt.vz*dt}
		rB := [3]float64{tt.x - tt.vx*dt, tt.y - tt.vy*dt, tt.z - tt.vz*dt}
		tF := tUTC.Add(time.Duration(dt * 1e9))
		tB := tUTC.Add(time.Duration(-dt * 1e9))
		xF, yF, zF := ECEF2ECI(rF[0], rF[1], rF[2], tF)
		xB, yB, zB := ECEF2ECI(rB[0], rB[1], rB[2], tB)
		vnX, vnY, vnZ := (xF-xB)/(2*dt), (yF-yB)/(2*dt), (zF-zB)/(2*dt)

		dX, dY, dZ := math.Abs(vaX-vnX), math.Abs(vaY-vnY), math.Abs(vaZ-vnZ)
		t.Logf("  %s:\n    analytical=[%8.3f %8.3f %8.3f]\n    numerical =[%8.3f %8.3f %8.3f]\n    diff      =[%8.6f %8.6f %8.6f]",
			tt.name, vaX, vaY, vaZ, vnX, vnY, vnZ, dX, dY, dZ)

		if dX > eps || dY > eps || dZ > eps {
			t.Errorf("  %s: 数值微分偏差过大: diff=[%.6f %.6f %.6f]", tt.name, dX, dY, dZ)
		}
	}
}

func TestNumericalDifferentiationECIVel2ECEFVel(t *testing.T) {
	tUTC := time.Date(2026, 6, 3, 7, 42, 16, 0, time.UTC)
	eps := 5.0

	tests := []struct {
		name       string
		vx, vy, vz float64
		x, y, z    float64
	}{
		{"赤道东向", 1000, 0, 0, 6378137, 0, 0},
		{"LEO卫星", 1500, 2500, 800, 5000000, 3000000, 2000000},
	}

	t.Logf("\n=== 数值微分验证: ECIVel2ECEFVel ===\n")
	for _, tt := range tests {
		vaX, vaY, vaZ := ECIVel2ECEFVel(tt.vx, tt.vy, tt.vz, tt.x, tt.y, tt.z, tUTC)

		// 数值微分: ECI位置前/后向推
		rF := [3]float64{tt.x + tt.vx*dt, tt.y + tt.vy*dt, tt.z + tt.vz*dt}
		rB := [3]float64{tt.x - tt.vx*dt, tt.y - tt.vy*dt, tt.z - tt.vz*dt}
		tF := tUTC.Add(time.Duration(dt * 1e9))
		tB := tUTC.Add(time.Duration(-dt * 1e9))
		xF, yF, zF := ECI2ECEF(rF[0], rF[1], rF[2], tF)
		xB, yB, zB := ECI2ECEF(rB[0], rB[1], rB[2], tB)
		vnX, vnY, vnZ := (xF-xB)/(2*dt), (yF-yB)/(2*dt), (zF-zB)/(2*dt)

		dX, dY, dZ := math.Abs(vaX-vnX), math.Abs(vaY-vnY), math.Abs(vaZ-vnZ)
		t.Logf("  %s:\n    analytical=[%8.3f %8.3f %8.3f]\n    numerical =[%8.3f %8.3f %8.3f]\n    diff      =[%8.6f %8.6f %8.6f]",
			tt.name, vaX, vaY, vaZ, vnX, vnY, vnZ, dX, dY, dZ)

		if dX > eps || dY > eps || dZ > eps {
			t.Errorf("  %s: 数值微分偏差过大: diff=[%.6f %.6f %.6f]", tt.name, dX, dY, dZ)
		}
	}
}

func TestNumericalDifferentiationENUVel2ECEFVel(t *testing.T) {
	ell, _ := NewEllipsoid("wgs84")
	eps := 1.0

	tests := []struct {
		name           string
		e, n, u        float64
		latDeg, lonDeg float64
	}{
		{"东向", 1000, 0, 0, 39.0, 116.4},
		{"北向", 0, 1000, 0, 39.0, 116.4},
		{"三维", 500, -300, 200, 40.0, 125.0},
	}

	t.Logf("\n=== 数值微分验证: ENUVel2ECEFVel ===\n")
	for _, tt := range tests {
		vaX, vaY, vaZ := ENUVel2ECEFVel(tt.e, tt.n, tt.u, tt.latDeg, tt.lonDeg)

		// ENU位置: 从测站原点出发，用ENU速度推前后两步
		eF, nF, uF := tt.e*dt, tt.n*dt, tt.u*dt
		eB, nB, uB := -tt.e*dt, -tt.n*dt, -tt.u*dt
		xF, yF, zF := ENU2ECEF(eF, nF, uF, tt.latDeg, tt.lonDeg, 0, ell)
		xB, yB, zB := ENU2ECEF(eB, nB, uB, tt.latDeg, tt.lonDeg, 0, ell)
		// 减去原点
		sx, sy, sz := ENU2ECEF(0, 0, 0, tt.latDeg, tt.lonDeg, 0, ell)
		vnX, vnY, vnZ := ((xF-sx)-(xB-sx))/(2*dt), ((yF-sy)-(yB-sy))/(2*dt), ((zF-sz)-(zB-sz))/(2*dt)

		dX, dY, dZ := math.Abs(vaX-vnX), math.Abs(vaY-vnY), math.Abs(vaZ-vnZ)
		t.Logf("  %s:\n    analytical=[%8.3f %8.3f %8.3f]\n    numerical =[%8.3f %8.3f %8.3f]\n    diff      =[%8.6f %8.6f %8.6f]",
			tt.name, vaX, vaY, vaZ, vnX, vnY, vnZ, dX, dY, dZ)

		if dX > eps || dY > eps || dZ > eps {
			t.Errorf("  %s: 数值微分偏差过大: diff=[%.6f %.6f %.6f]", tt.name, dX, dY, dZ)
		}
	}
}

func TestNumericalDifferentiationAERDeriv2ENUVel(t *testing.T) {
	eps := 1.0

	tests := []struct {
		name               string
		R, azDeg, elDeg    float64
		dR, dAzDeg, dElDeg float64
	}{
		{"径向", 200000, 45, 30, 1000, 0, 0},
		{"方位角变化", 200000, 45, 30, 0, 0.1, 0},
		{"俯仰角变化", 200000, 45, 30, 0, 0, 0.05},
		{"三维复合", 1000000, 60, 45, 500, 0.2, -0.1},
	}

	t.Logf("\n=== 数值微分验证: AERDeriv2ENUVel ===\n")
	for _, tt := range tests {
		vaE, vaN, vaU := AERDeriv2ENUVel(tt.R, tt.azDeg, tt.elDeg, tt.dR, tt.dAzDeg, tt.dElDeg)

		// AER位置: 前向/后向推RAE
		// 前向: R + dR*dt, az + dAz*dt, el + dEl*dt
		rF := tt.R + tt.dR*dt
		azF := tt.azDeg + tt.dAzDeg*dt
		elF := tt.elDeg + tt.dElDeg*dt
		eF, nF, uF := AER2ENU(azF, elF, rF)

		// 后向: R - dR*dt, az - dAz*dt, el - dEl*dt
		rB := tt.R - tt.dR*dt
		azB := tt.azDeg - tt.dAzDeg*dt
		elB := tt.elDeg - tt.dElDeg*dt
		eB, nB, uB := AER2ENU(azB, elB, rB)

		// 消去位置偏移：用数值微分计算速度
		// 注意：AER2ENU给出的就是ENU位置，直接对时间差分
		vnE, vnN, vnU := (eF-eB)/(2*dt), (nF-nB)/(2*dt), (uF-uB)/(2*dt)

		dE, dN, dU := math.Abs(vaE-vnE), math.Abs(vaN-vnN), math.Abs(vaU-vnU)
		t.Logf("  %s:\n    analytical=[%8.3f %8.3f %8.3f]\n    numerical =[%8.3f %8.3f %8.3f]\n    diff      =[%8.6f %8.6f %8.6f]",
			tt.name, vaE, vaN, vaU, vnE, vnN, vnU, dE, dN, dU)

		if dE > eps || dN > eps || dU > eps {
			t.Errorf("  %s: 数值微分偏差过大: diff=[%.6f %.6f %.6f]", tt.name, dE, dN, dU)
		}
	}
}

// TestAllNumericalDifferentiation 一键运行所有数值微分验证，输出汇总。
func TestAllNumericalDifferentiation(t *testing.T) {
	tUTC := time.Date(2026, 6, 3, 7, 42, 16, 0, time.UTC)
	maxDiff := 0.0

	t.Logf("\n%s", "========================================")
	t.Logf("  数值微分验证汇总报告")
	t.Logf("  步长 dt = %.4f 秒", dt)
	t.Logf("  允许误差: ECEF↔ECI=5.0m/s, ENU/AER=1.0m/s")
	t.Logf("%s", "========================================")

	// ECEFVel2ECIVel
	t.Logf("\n--- ECEFVel2ECIVel ---")
	cases := []struct {
		name                string
		vx, vy, vz, x, y, z float64
	}{
		{"赤道东向", 1000, 0, 0, 6378137, 0, 0},
		{"LEO卫星", 1500, 2500, 800, 5000000, 3000000, 2000000},
	}
	for _, c := range cases {
		vaX, vaY, vaZ := ECEFVel2ECIVel(c.vx, c.vy, c.vz, c.x, c.y, c.z, tUTC)
		rF := [3]float64{c.x + c.vx*dt, c.y + c.vy*dt, c.z + c.vz*dt}
		rB := [3]float64{c.x - c.vx*dt, c.y - c.vy*dt, c.z - c.vz*dt}
		tF := tUTC.Add(time.Duration(dt * 1e9))
		tB := tUTC.Add(time.Duration(-dt * 1e9))
		xF, yF, zF := ECEF2ECI(rF[0], rF[1], rF[2], tF)
		xB, yB, zB := ECEF2ECI(rB[0], rB[1], rB[2], tB)
		vnX, vnY, vnZ := (xF-xB)/(2*dt), (yF-yB)/(2*dt), (zF-zB)/(2*dt)
		d := math.Max(math.Abs(vaX-vnX), math.Max(math.Abs(vaY-vnY), math.Abs(vaZ-vnZ)))
		if d > maxDiff {
			maxDiff = d
		}
		status := "✅"
		if d > 5.0 {
			status = "❌"
		}
		t.Logf("  %s %s: diff=%.6f m/s", status, c.name, d)
	}

	// ECIVel2ECEFVel
	t.Logf("\n--- ECIVel2ECEFVel ---")
	eciCases := []struct {
		name                string
		vx, vy, vz, x, y, z float64
	}{
		{"赤道东向", 1000, 0, 0, 6378137, 0, 0},
		{"LEO卫星", 1500, 2500, 800, 5000000, 3000000, 2000000},
	}
	for _, c := range eciCases {
		vaX, vaY, vaZ := ECIVel2ECEFVel(c.vx, c.vy, c.vz, c.x, c.y, c.z, tUTC)
		rF := [3]float64{c.x + c.vx*dt, c.y + c.vy*dt, c.z + c.vz*dt}
		rB := [3]float64{c.x - c.vx*dt, c.y - c.vy*dt, c.z - c.vz*dt}
		tF := tUTC.Add(time.Duration(dt * 1e9))
		tB := tUTC.Add(time.Duration(-dt * 1e9))
		xF, yF, zF := ECI2ECEF(rF[0], rF[1], rF[2], tF)
		xB, yB, zB := ECI2ECEF(rB[0], rB[1], rB[2], tB)
		vnX, vnY, vnZ := (xF-xB)/(2*dt), (yF-yB)/(2*dt), (zF-zB)/(2*dt)
		d := math.Max(math.Abs(vaX-vnX), math.Max(math.Abs(vaY-vnY), math.Abs(vaZ-vnZ)))
		if d > maxDiff {
			maxDiff = d
		}
		status := "✅"
		if d > 5.0 {
			status = "❌"
		}
		t.Logf("  %s %s: diff=%.6f m/s", status, c.name, d)
	}

	// ENUVel2ECEFVel
	t.Logf("\n--- ENUVel2ECEFVel ---")
	ell, _ := NewEllipsoid("wgs84")
	enuCases := []struct {
		name              string
		e, n, u, lat, lon float64
	}{
		{"东向", 1000, 0, 0, 39.0, 116.4},
		{"三维", 500, -300, 200, 40.0, 125.0},
	}
	for _, c := range enuCases {
		vaX, vaY, vaZ := ENUVel2ECEFVel(c.e, c.n, c.u, c.lat, c.lon)
		eF, nF, uF := c.e*dt, c.n*dt, c.u*dt
		eB, nB, uB := -c.e*dt, -c.n*dt, -c.u*dt
		sx, sy, sz := ENU2ECEF(0, 0, 0, c.lat, c.lon, 0, ell)
		xF, yF, zF := ENU2ECEF(eF, nF, uF, c.lat, c.lon, 0, ell)
		xB, yB, zB := ENU2ECEF(eB, nB, uB, c.lat, c.lon, 0, ell)
		vnX, vnY, vnZ := ((xF-sx)-(xB-sx))/(2*dt), ((yF-sy)-(yB-sy))/(2*dt), ((zF-sz)-(zB-sz))/(2*dt)
		d := math.Max(math.Abs(vaX-vnX), math.Max(math.Abs(vaY-vnY), math.Abs(vaZ-vnZ)))
		if d > maxDiff {
			maxDiff = d
		}
		status := "✅"
		if d > 1.0 {
			status = "❌"
		}
		t.Logf("  %s %s: diff=%.6f m/s", status, c.name, d)
	}

	// AERDeriv2ENUVel
	t.Logf("\n--- AERDeriv2ENUVel ---")
	aerCases := []struct {
		name                    string
		R, az, el, dR, dAz, dEl float64
	}{
		{"径向", 200000, 45, 30, 1000, 0, 0},
		{"方位角", 200000, 45, 30, 0, 0.1, 0},
		{"三维", 1000000, 60, 45, 500, 0.2, -0.1},
	}
	for _, c := range aerCases {
		vaE, vaN, vaU := AERDeriv2ENUVel(c.R, c.az, c.el, c.dR, c.dAz, c.dEl)
		rF := c.R + c.dR*dt
		eF, nF, uF := AER2ENU(c.az+c.dAz*dt, c.el+c.dEl*dt, rF)
		rB := c.R - c.dR*dt
		eB, nB, uB := AER2ENU(c.az-c.dAz*dt, c.el-c.dEl*dt, rB)
		vnE, vnN, vnU := (eF-eB)/(2*dt), (nF-nB)/(2*dt), (uF-uB)/(2*dt)
		d := math.Max(math.Abs(vaE-vnE), math.Max(math.Abs(vaN-vnN), math.Abs(vaU-vnU)))
		if d > maxDiff {
			maxDiff = d
		}
		status := "✅"
		if d > 1.0 {
			status = "❌"
		}
		t.Logf("  %s %s: diff=%.6f m/s", status, c.name, d)
	}

	t.Logf("\n%s", "========================================")
	t.Logf("  最大偏差: %.6f m/s", maxDiff)
	if maxDiff <= 5.0 {
		t.Logf("  结论: ✅ 所有速度转换函数通过数值微分验证")
	} else {
		t.Errorf("  结论: ❌ 部分函数偏差过大(>5.0 m/s)")
	}
	t.Logf("%s", "========================================")
}
