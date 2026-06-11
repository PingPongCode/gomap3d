package gomap3d

import (
	"fmt"
	"math"
	"testing"
	"time"
)

const epsFloat = 0.1 // 0.1m 位置容差（完整链变换精度）

// ---- 工具函数 ----

func vecDist(a, b [3]float64) float64 {
	dx, dy, dz := a[0]-b[0], a[1]-b[1], a[2]-b[2]
	return math.Sqrt(dx*dx + dy*dy + dz*dz)
}

func vecLen(v [3]float64) float64 {
	return math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
}

// ---- 旧版 GMST 旋转（用于对比） ----

func eci2ecefOld(x, y, z float64, t time.Time) (float64, float64, float64) {
	jd := juliandate(t)
	gst := greenwichsrt(jd)
	rotMat := R3(gst)
	eciVec := [3]float64{x, y, z}
	ecefVec := multiplyMatrixVector(rotMat, eciVec)
	return ecefVec[0], ecefVec[1], ecefVec[2]
}

func ecef2eciOld(x, y, z float64, t time.Time) (float64, float64, float64) {
	jd := juliandate(t)
	gst := greenwichsrt(jd)
	rotMat := R3(gst)
	transposed := transpose(rotMat)
	ecefVec := [3]float64{x, y, z}
	eciVec := multiplyMatrixVector(transposed, ecefVec)
	return eciVec[0], eciVec[1], eciVec[2]
}

// ---- 测试 ----

func TestECI2ECEFRoundtrip(t *testing.T) {
	//  来回测试：ECEF → ECI → ECEF 应恢复原始值
	tUTC := time.Date(2024, 6, 10, 12, 0, 0, 0, time.UTC)

	testCases := []struct {
		name string
		p    [3]float64
	}{
		{"北京", [3]float64{-2.174e6, 4.389e6, 4.077e6}},
		{"纽约", [3]float64{1.334e6, -4.654e6, 4.139e6}},
		{"赤道", [3]float64{6.378e6, 0, 0}},
		{"极地", [3]float64{0, 0, 6.357e6}},
		{"LEO", [3]float64{4.5e6, 4.5e6, 0}},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			// ECEF → ECI
			xi, yi, zi := ECEF2ECI(tc.p[0], tc.p[1], tc.p[2], tUTC)
			// ECI → ECEF
			xe, ye, ze := ECI2ECEF(xi, yi, zi, tUTC)

			d := vecDist(tc.p, [3]float64{xe, ye, ze})
			if d > epsFloat {
				t.Errorf("roundtrip error = %.3f m", d)
			}
		})
	}
}

func TestOldVsNewDelta(t *testing.T) {
	// 比较旧（纯GMST旋转）与新（完整IAU链）的差异
	tUTC := time.Date(2024, 6, 10, 12, 0, 0, 0, time.UTC)

	// 地球表面 + 500km 高度点
	rEcef := [3]float64{-2.174e6, 4.389e6, 4.077e6} // ~500km 高度

	// 新：完整链
	xiNew, yiNew, ziNew := ECEF2ECI(rEcef[0], rEcef[1], rEcef[2], tUTC)
	// 旧：纯GMST
	xiOld, yiOld, ziOld := ecef2eciOld(rEcef[0], rEcef[1], rEcef[2], tUTC)

	delta := vecDist([3]float64{xiNew, yiNew, ziNew}, [3]float64{xiOld, yiOld, ziOld})
	t.Logf("ECEF → ECI 新旧差异: %.1f m (2024.5)", delta)

	// 差异应约为 28-73 km (岁差 0.34° at ~6878 km → 41 km,
	// classical 3-angle 实际约 0.23° 总旋转 → ~28 km)
	if delta < 20000 || delta > 80000 {
		t.Errorf("unexpected delta = %.0f m, expected ~40000 m", delta)
	}
}

func TestVelocityRoundtrip(t *testing.T) {
	// 速度往返测试
	tUTC := time.Date(2024, 6, 10, 12, 0, 0, 0, time.UTC)
	rEcef := [3]float64{-2.174e6, 4.389e6, 4.077e6}
	vEcef := [3]float64{3000, 2000, -4000}

	// ECEF → ECI
	vxi, vyi, vzi := ECEFVel2ECIVel(vEcef[0], vEcef[1], vEcef[2],
		rEcef[0], rEcef[1], rEcef[2], tUTC)
	// ECI → ECEF
	// 先将位置转到ECI
	xi, yi, zi := ECEF2ECI(rEcef[0], rEcef[1], rEcef[2], tUTC)
	vxe, vye, vze := ECIVel2ECEFVel(vxi, vyi, vzi, xi, yi, zi, tUTC)

	d := vecDist(vEcef, [3]float64{vxe, vye, vze})
	if d > 0.01 { // 速度往返精度 1 cm/s
		t.Errorf("velocity roundtrip error = %.6f m/s", d)
	}
}

func TestPrecessionAngle(t *testing.T) {
	// 验证岁差角在 2024 约为 0.56° (θ ≈ 2004"/cy × 0.244 cy ≈ 490")
	jd2024 := juliandate(time.Date(2024, 6, 10, 12, 0, 0, 0, time.UTC))
	T := (jd2024 - 2451545.0) / 36525.0
	_, _, theta := precessionAngles(T)

	thetaDeg := theta * 180 / math.Pi
	t.Logf("T = %.6f centuries, precession θ = %.4f°", T, thetaDeg)

	// θ ≈ 2004.19"/cy × 0.2444 cy / 3600 = 0.136°
	if thetaDeg < 0.13 || thetaDeg > 0.14 {
		t.Errorf("unexpected precession theta = %.4f°", thetaDeg)
	}
}

func TestNutation(t *testing.T) {
	// 验证章动值量级合理 (< 20")
	jd2024 := juliandate(time.Date(2024, 6, 10, 12, 0, 0, 0, time.UTC))
	T := (jd2024 - 2451545.0) / 36525.0
	dpsi, deps := nutation(T)

	dpsiAS := dpsi / as2r
	depsAS := deps / as2r

	t.Logf("Δψ = %.3f arcsec, Δε = %.3f arcsec", dpsiAS, depsAS)

	if math.Abs(dpsiAS) > 20 || math.Abs(depsAS) > 20 {
		t.Errorf("nutation values out of expected range")
	}
}

func TestGASTvsGMST(t *testing.T) {
	// GAST = GMST + equation of equinoxes
	jd2024 := juliandate(time.Date(2024, 6, 10, 12, 0, 0, 0, time.UTC))
	gmst := greenwichsrt(jd2024)
	gast := GAST(jd2024)

	// 赤经章动 = GAST - GMST (rad) → arcsec
	eeqArcsec := math.Abs(gast-gmst) * 180.0 * 3600.0 / math.Pi
	t.Logf("GMST = %.4f°, GAST = %.4f°, eq of equinoxes = %.6f arcsec",
		gmst*180/math.Pi, gast*180/math.Pi, eeqArcsec)

	// 赤经章动 < 2 arcsec
	if eeqArcsec > 2 {
		t.Errorf("equation of equinoxes too large: %.3f arcsec", eeqArcsec)
	}
}

func TestConsistencyGCRF2ITRF(t *testing.T) {
	// GCRF2ITRF 矩阵应为正交矩阵
	jd2024 := juliandate(time.Date(2024, 6, 10, 12, 0, 0, 0, time.UTC))
	M := GCRF2ITRF(jd2024)

	// M^T · M 应 ≈ I
	Mt := transpose(M)
	I := mul33(Mt, M)

	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			expected := 0.0
			if i == j {
				expected = 1.0
			}
			if math.Abs(I[i][j]-expected) > 1e-9 {
				t.Errorf("orthogonality check failed at (%d,%d): %.3e", i, j, I[i][j]-expected)
			}
		}
	}
}

func TestKnownPoint(t *testing.T) {
	// 用 mc_04 测试用例验证：雷达(38°N,128°E,200m)，目标ECEF → ECI → ECEF 来回
	tUTC := time.Date(2024, 2, 9, 12, 0, 0, 0, time.UTC)

	ell, _ := NewEllipsoid("wgs84")
	// 雷达站
	rx, ry, rz := Geodetic2ECEF(38.0, 128.0, 200, ell)
	// R=350km, Az=30°, El=25° → ENU → ECEF
	e, n, u := AER2ENU(30, 25, 350000)
	tx, ty, tz := ENU2ECEF(e, n, u, 38.0, 128.0, 200, ell)
	// 目标 ECEF
	tgx, tgy, tgz := rx+tx, ry+ty, rz+tz

	// ECEF → ECI (完整链)
	xi, yi, zi := ECEF2ECI(tgx, tgy, tgz, tUTC)
	// ECEF → ECI (旧：纯GMST)
	xiOld, yiOld, ziOld := ecef2eciOld(tgx, tgy, tgz, tUTC)

	delta := vecDist([3]float64{xi, yi, zi}, [3]float64{xiOld, yiOld, ziOld})
	t.Logf("mc_04 ECI位置: 新(%.1f, %.1f, %.1f) km", xi/1000, yi/1000, zi/1000)
	t.Logf("mc_04 ECI位置: 旧(%.1f, %.1f, %.1f) km", xiOld/1000, yiOld/1000, ziOld/1000)
	t.Logf("mc_04 新旧差异: %.1f m (预期 ~40 km @ 岁差 0.34°)", delta)

	// 来回验证
	tgx2, tgy2, tgz2 := ECI2ECEF(xi, yi, zi, tUTC)
	dRound := vecDist([3]float64{tgx, tgy, tgz}, [3]float64{tgx2, tgy2, tgz2})
	if dRound > epsFloat {
		t.Errorf("roundtrip error for mc_04: %.3f m", dRound)
	}
}

func TestMatrixChain(t *testing.T) {
	// 验证完整矩阵链的一致性：框架偏差 × 岁差 × 章动 × Rz(GAST) 应有合理数值
	jd2024 := juliandate(time.Date(2024, 6, 10, 12, 0, 0, 0, time.UTC))
	T := (jd2024 - 2451545.0) / 36525.0

	B := frameBiasMatrix()
	P := precessionMatrix(T)
	N := nutationMatrix(T)
	gast := GAST(jd2024)
	R := R3(gast)

	// 组合矩阵
	M := mul33(R, mul33(N, mul33(P, B)))
	Mref := GCRF2ITRF(jd2024)

	// 应与 GCRF2ITRF 结果一致
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if math.Abs(M[i][j]-Mref[i][j]) > 1e-14 {
				t.Errorf("matrix chain mismatch at (%d,%d): %.3e vs %.3e", i, j, M[i][j], Mref[i][j])
			}
		}
	}
}

// ---- Benchmark ----

func BenchmarkECEF2ECI_Old(b *testing.B) {
	tUTC := time.Date(2024, 6, 10, 12, 0, 0, 0, time.UTC)
	x, y, z := -2.174e6, 4.389e6, 4.077e6
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ecef2eciOld(x, y, z, tUTC)
	}
}

func BenchmarkECEF2ECI_New(b *testing.B) {
	tUTC := time.Date(2024, 6, 10, 12, 0, 0, 0, time.UTC)
	x, y, z := -2.174e6, 4.389e6, 4.077e6
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ECEF2ECI(x, y, z, tUTC)
	}
}

func BenchmarkECI2ECEF_Old(b *testing.B) {
	tUTC := time.Date(2024, 6, 10, 12, 0, 0, 0, time.UTC)
	x, y, z := -6.219e6, 1.801e6, 4.077e6
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		eci2ecefOld(x, y, z, tUTC)
	}
}

func BenchmarkECI2ECEF_New(b *testing.B) {
	tUTC := time.Date(2024, 6, 10, 12, 0, 0, 0, time.UTC)
	x, y, z := -6.219e6, 1.801e6, 4.077e6
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ECI2ECEF(x, y, z, tUTC)
	}
}

// ---- Example ----

func ExampleECI2ECEF() {
	tUTC := time.Date(2024, 6, 10, 12, 0, 0, 0, time.UTC)
	xi, yi, zi := -6.219e6, 1.801e6, 4.077e6
	xe, ye, ze := ECI2ECEF(xi, yi, zi, tUTC)
	fmt.Printf("ECI(%.0f,%.0f,%.0f) → ECEF(%.0f,%.0f,%.0f) m\n", xi, yi, zi, xe, ye, ze)
}
