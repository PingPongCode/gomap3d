package gomap3d

import (
	"math"
	"time"
)

// ===== 地球自转常数 =====
const (
	// We 地球自转角速度 (rad/s)
	We = 7.2921150e-5
)

// ===== ECEF 速度 ↔ ECI 速度 =====

// ECEFVel2ECIVel 将 ECEF 速度转换为 ECI 速度。
//
// 公式: v_eci = M^T · (v_ecef + ω × r_ecef)
// 默认 M = Rz(GMST)（GMST 模式），SetGMSTMode(false) 后 M = GCRF2ITRF
func ECEFVel2ECIVel(vx, vy, vz, x, y, z float64, t time.Time) (vxEci, vyEci, vzEci float64) {
	jd := juliandate(t)
	var M [3][3]float64
	if useSimpleGMST {
		gmst := greenwichsrt(jd)
		cosG := math.Cos(gmst)
		sinG := math.Sin(gmst)
		M = [3][3]float64{
			{cosG, sinG, 0},
			{-sinG, cosG, 0},
			{0, 0, 1},
		}
	} else {
		M = GCRF2ITRF(jd)
		M = transpose(M)
	}

	// ω × r_ecef
	wxr := [3]float64{
		-We * y,
		We * x,
		0,
	}

	// v_ecef + ω × r_ecef
	v := [3]float64{
		vx + wxr[0],
		vy + wxr[1],
		vz + wxr[2],
	}

	// M^T · (v_ecef + ω × r_ecef)
	result := multiplyMatrixVector(M, v)
	return result[0], result[1], result[2]
}

// ECIVel2ECEFVel 将 ECI 速度转换为 ECEF 速度。
//
// 公式: v_ecef = M · v_eci - ω × r_ecef
// 其中 r_ecef = M · r_eci
// 默认 M = Rz(GMST)（GMST 模式），SetGMSTMode(false) 后 M = GCRF2ITRF
func ECIVel2ECEFVel(vx, vy, vz, x, y, z float64, t time.Time) (vxEcef, vyEcef, vzEcef float64) {
	jd := juliandate(t)
	var M [3][3]float64
	if useSimpleGMST {
		gmst := greenwichsrt(jd)
		cosG := math.Cos(gmst)
		sinG := math.Sin(gmst)
		// Rz(+GMST)
		M = [3][3]float64{
			{cosG, -sinG, 0},
			{sinG, cosG, 0},
			{0, 0, 1},
		}
	} else {
		M = GCRF2ITRF(jd)
	}

	// r_ecef = M · r_eci
	rEci := [3]float64{x, y, z}
	rEcef := multiplyMatrixVector(M, rEci)

	// ω × r_ecef
	wxr := [3]float64{
		-We * rEcef[1],
		We * rEcef[0],
		0,
	}

	// v_ecef_rot = M · v_eci
	vEci := [3]float64{vx, vy, vz}
	vRot := multiplyMatrixVector(M, vEci)

	// v_ecef = M · v_eci - ω × r_ecef
	return vRot[0] - wxr[0], vRot[1] - wxr[1], vRot[2] - wxr[2]
}

// ===== ENU 速度 ↔ ECEF 速度 =====

// ENUVel2ECEFVel 将 ENU 速度转换为 ECEF 速度。
//
// 原理: 将速度矢量 (eVel, nVel, uVel) 视为 ENU 坐标点，
// 通过 ENU→ECEF 转换得到 ECEF 坐标，再减去原点得到纯速度矢量。
//
// 输入:
//   - eVel, nVel, uVel: ENU 速度 (m/s)
//   - latDeg, lonDeg: 测站经纬度 (度)
//
// 输出: ECEF 速度 (vx, vy, vz) (m/s)
func ENUVel2ECEFVel(eVel, nVel, uVel, latDeg, lonDeg float64) (vx, vy, vz float64) {
	ell, _ := NewEllipsoid("wgs84")
	// 站址原点对应的 ECEF 坐标
	sx, sy, sz := ENU2ECEF(0, 0, 0, latDeg, lonDeg, 0, ell)
	// 含速度偏移的 ECEF 坐标
	x, y, z := ENU2ECEF(eVel, nVel, uVel, latDeg, lonDeg, 0, ell)
	return x - sx, y - sy, z - sz
}

// ECEFVel2ENUVel 将 ECEF 速度转换为 ENU 速度。
//
// 原理: 将站址原点加上 ECEF 速度矢量视为一个 ECEF 坐标点，
// 通过 ECEF→ENU 转换得到 ENU 坐标。
//
// 输入:
//   - vx, vy, vz: ECEF 速度 (m/s)
//   - latDeg, lonDeg: 测站经纬度 (度)
//
// 输出: ENU 速度 (eV, nV, uV) (m/s)
func ECEFVel2ENUVel(vx, vy, vz, latDeg, lonDeg float64) (eV, nV, uV float64) {
	ell, _ := NewEllipsoid("wgs84")
	// 站址原点对应的 ECEF 坐标
	sx, sy, sz := ENU2ECEF(0, 0, 0, latDeg, lonDeg, 0, ell)
	// 原点 + ECEF 速度矢量 → 计算 ENU 分量
	e, n, u := ECEF2ENU(sx+vx, sy+vy, sz+vz, latDeg, lonDeg, 0, ell)
	return e, n, u
}

// ===== AER 变化率（RAE 导数）↔ ENU 速度 =====

// AERDeriv2ENUVel 将 AER 变化率（RAE 导数）转换为 ENU 速度。
//
// 输入:
//   - R: 斜距 (m)
//   - azDeg: 方位角 (deg)
//   - elDeg: 俯仰角 (deg)
//   - dR: 斜距变化率 (m/s)
//   - dAzDeg: 方位角变化率 (deg/s)
//   - dElDeg: 俯仰角变化率 (deg/s)
//
// 输出: ENU 速度 (e, n, u) (m/s)
//
// 推导:
//   e  = R·cos(El)·sin(Az)
//   n  = R·cos(El)·cos(Az)
//   u  = R·sin(El)
//   ė = Ṙ·cos(El)·sin(Az) - R·sin(El)·Ėl·sin(Az) + R·cos(El)·cos(Az)·Åz
//   ṅ = Ṙ·cos(El)·cos(Az) - R·sin(El)·Ėl·cos(Az) - R·cos(El)·sin(Az)·Åz
//   u̇ = Ṙ·sin(El) + R·cos(El)·Ėl
func AERDeriv2ENUVel(R, azDeg, elDeg, dR, dAzDeg, dElDeg float64) (e, n, u float64) {
	azRad := azDeg * math.Pi / 180
	elRad := elDeg * math.Pi / 180
	dAzRad := dAzDeg * math.Pi / 180
	dElRad := dElDeg * math.Pi / 180

	cosEl := math.Cos(elRad)
	sinEl := math.Sin(elRad)
	cosAz := math.Cos(azRad)
	sinAz := math.Sin(azRad)

	e = dR*cosEl*sinAz - R*sinEl*dElRad*sinAz + R*cosEl*cosAz*dAzRad
	n = dR*cosEl*cosAz - R*sinEl*dElRad*cosAz - R*cosEl*sinAz*dAzRad
	u = dR*sinEl + R*cosEl*dElRad
	return
}

// ENUVel2AERDeriv 将 ENU 速度转换为 AER 变化率（RAE 导数）。
//
// 输入:
//   - e, n, u: ENU 速度 (m/s)
//   - R: 斜距 (m)
//   - azDeg: 方位角 (deg)
//   - elDeg: 俯仰角 (deg)
//
// 输出:
//   - dR: 斜距变化率 (m/s)
//   - dAzDeg: 方位角变化率 (deg/s)
//   - dElDeg: 俯仰角变化率 (deg/s)
//
// 推导:
//   R = √(e² + n² + u²)
//   Ṙ  = eVel·cos(El)·sin(Az) + nVel·cos(El)·cos(Az) + uVel·sin(El)
//        = (ePos·eVel + nPos·nVel + uPos·uVel) / R
//       即速度矢量在视线方向上的投影
//
//   Ėl = (uVel - Ṙ·sin(El)) / (R·cos(El))
//   Åz = (eVel·cos(Az) - nVel·sin(Az)) / (R·cos(El))
func ENUVel2AERDeriv(eVel, nVel, uVel, R, azDeg, elDeg float64) (dR, dAzDeg, dElDeg float64) {
	azRad := azDeg * math.Pi / 180
	elRad := elDeg * math.Pi / 180
	cosEl := math.Cos(elRad)
	sinEl := math.Sin(elRad)
	cosAz := math.Cos(azRad)
	sinAz := math.Sin(azRad)

	// Ṙ = eVel·cos(El)·sin(Az) + nVel·cos(El)·cos(Az) + uVel·sin(El)
	//   速度矢量在视线单位方向上的投影
	dR = eVel*cosEl*sinAz + nVel*cosEl*cosAz + uVel*sinEl

	// Ėl = (uVel - Ṙ·sin(El)) / (R·cos(El))
	if math.Abs(cosEl) > 1e-12 {
		dElRad := (uVel - dR*sinEl) / (R * cosEl)
		dElDeg = dElRad * 180 / math.Pi
	}

	// Åz = (eVel·cos(Az) - nVel·sin(Az)) / (R·cos(El))
	if math.Abs(cosEl) > 1e-12 {
		dAzRad := (eVel*cosAz - nVel*sinAz) / (R * cosEl)
		dAzDeg = dAzRad * 180 / math.Pi
	}

	return
}

// ===== 便捷方法：一条龙 RAE 导数 → ECI 速度 =====

// AERDeriv2ECIVel 将 AER 变化率（RAE 导数）一步转换为 ECI 速度。
//
// 转换链路: AER 变化率 → ENU 速度 → ECEF 速度 → ECI 速度
//
// 输入:
//   - R: 斜距 (m)
//   - azDeg: 方位角 (deg)
//   - elDeg: 俯仰角 (deg)
//   - dR: 斜距变化率 (m/s)
//   - dAzDeg: 方位角变化率 (deg/s)
//   - dElDeg: 俯仰角变化率 (deg/s)
//   - latDeg, lonDeg: 测站经纬度 (度)
//   - alt: 测站海拔 (m)
//   - t: UTC 时间
//
// 输出: ECI 速度 (vx, vy, vz) (m/s)
func AERDeriv2ECIVel(R, azDeg, elDeg, dR, dAzDeg, dElDeg, latDeg, lonDeg, alt float64, t time.Time) (vx, vy, vz float64) {
	// 1. AER 变化率 → ENU 速度
	eVel, nVel, uVel := AERDeriv2ENUVel(R, azDeg, elDeg, dR, dAzDeg, dElDeg)

	// 2. ENU 速度 → ECEF 速度
	vxE, vyE, vzE := ENUVel2ECEFVel(eVel, nVel, uVel, latDeg, lonDeg)

	// 3. 计算 ECEF 位置（用于地球自转修正）
	ell, _ := NewEllipsoid("wgs84")
	azRad := azDeg * math.Pi / 180
	elRad := elDeg * math.Pi / 180
	ex := R * math.Cos(elRad) * math.Sin(azRad)
	nx := R * math.Cos(elRad) * math.Cos(azRad)
	ux := R * math.Sin(elRad)
	px, py, pz := ENU2ECEF(ex, nx, ux, latDeg, lonDeg, alt, ell)

	// 4. ECEF 速度 → ECI 速度
	return ECEFVel2ECIVel(vxE, vyE, vzE, px, py, pz, t)
}

// ECIVel2AERDeriv 将 ECI 速度一步反算为 AER 变化率（RAE 导数）。
//
// 转换链路: ECI 速度 → ECEF 速度 → ENU 速度 → AER 变化率
//
// 输入:
//   - vx, vy, vz: ECI 速度 (m/s)
//   - Rx, Ry, Rz: ECI 位置 (m)（用于地球自转修正和 RAE 计算）
//   - latDeg, lonDeg: 测站经纬度 (度)
//   - alt: 测站海拔 (m)
//   - t: UTC 时间
//
// 输出:
//   - dR: 斜距变化率 (m/s)
//   - dAzDeg: 方位角变化率 (deg/s)
//   - dElDeg: 俯仰角变化率 (deg/s)
func ECIVel2AERDeriv(vx, vy, vz, Rx, Ry, Rz, latDeg, lonDeg, alt float64, t time.Time) (dR, dAzDeg, dElDeg float64) {
	// 1. ECI 速度 → ECEF 速度
	vxE, vyE, vzE := ECIVel2ECEFVel(vx, vy, vz, Rx, Ry, Rz, t)

	// 2. ECEF 速度 → ENU 速度
	eVel, nVel, uVel := ECEFVel2ENUVel(vxE, vyE, vzE, latDeg, lonDeg)

	// 3. 计算目标到测站的 RAE 参数（用于确定几何关系）
	ell, _ := NewEllipsoid("wgs84")

	// ECI 位置 → ECEF 位置
	rEcefX, rEcefY, rEcefZ := ECI2ECEF(Rx, Ry, Rz, t)

	// ECEF 位置 → ENU 位置（使用实际测站海拔）
	en, nn, un := ECEF2ENU(rEcefX, rEcefY, rEcefZ, latDeg, lonDeg, alt, ell)

	// ENU 位置 → RAE
	azDeg, elDeg, srange := ENU2AER(en, nn, un)

	// 4. ENU 速度 → AER 变化率
	return ENUVel2AERDeriv(eVel, nVel, uVel, srange, azDeg, elDeg)
}
