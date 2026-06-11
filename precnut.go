package gomap3d

import (
	"math"
)

// ============================================================
// IAU 2006 Precession + IAU 2000B Nutation + Frame Bias
// 实现 Aerospace Toolbox 等效的 GCRF ↔ ITRF 完整变换链
//
// 变换链: ECI(GCRF) → [frame bias] → [precession] → [nutation] → [GAST rotation] → ECEF(ITRF)
// 参考: Vallado "Fundamentals of Astrodynamics", IERS Conventions 2010
// ============================================================

const (
	// arcsec 到 rad 统一转换
	as2r = math.Pi / (180.0 * 3600.0)
	// 角秒转度
	as2d = 1.0 / 3600.0
)

// ---- Frame Bias (IERS Conventions 2010) ----
// GCRF → J2000.0 mean equator/equinox 的常值偏移
// 3 个独立旋转角 (mas, 即 milliarcsecond)
const (
	frameBiasDX = -0.0146  // dα0, mas, 转为角秒 = -14.6 mas
	frameBiasDE = -0.016617 // ξ0,  mas, 转为角秒 = -16.617 mas
	frameBiasDP = -0.0068192 // η0,  mas, 转为角秒 =  -6.8192 mas
)

// frameBiasMatrix 计算 GCRF→J2000 框架偏差矩阵
// B = Rz(-dα0) · Ry(η0) · Rx(-ξ0)
func frameBiasMatrix() [3][3]float64 {
	// mas → rad
	dAlpha := frameBiasDX * 1e-3 * as2r // mas → rad
	xi := frameBiasDE * 1e-3 * as2r     // mas → rad
	eta := frameBiasDP * 1e-3 * as2r    // mas → rad

	sda, cda := math.Sin(dAlpha), math.Cos(dAlpha)
	sxi, cxi := math.Sin(xi), math.Cos(xi)
	seta, ceta := math.Sin(eta), math.Cos(eta)

	// B = Rz(-dα0) · Ry(η0) · RX(-ξ0)
	return [3][3]float64{
		{cda*ceta + sda*sxi*seta, -sda*ceta + cda*sxi*seta, cxi * seta},
		{sda * cxi, cda * cxi, sxi},
		{-cda*seta + sda*sxi*ceta, sda*seta + cda*sxi*ceta, cxi * ceta},
	}
}

// ---- IAU 2006 Precession (Classical 3-angle: ζ, z, θ) ----
// P(t) = Rz(−z_A) · Ry(θ_A) · Rz(−ζ_A)
// 参考: Vallado 4th ed. Table 3-5, IERS Conventions 2010
func precessionAngles(T float64) (zeta, z, theta float64) {
	zeta = (2.650545 +
		2306.083227*T +
		0.2988499*T*T +
		0.01801828*T*T*T -
		0.000005971*T*T*T*T -
		0.0000003173*T*T*T*T*T) * as2r

	theta = (2004.191903*T -
		0.4294934*T*T -
		0.04182264*T*T*T -
		0.000007089*T*T*T*T -
		0.0000001274*T*T*T*T*T) * as2r

	z = (-2.650545 +
		2306.077181*T +
		1.0927348*T*T +
		0.01826837*T*T*T -
		0.000028596*T*T*T*T -
		0.0000002904*T*T*T*T*T) * as2r

	return
}

func precessionMatrix(T float64) [3][3]float64 {
	zeta, z, theta := precessionAngles(T)

	// P = Rz(−z) · Ry(θ) · Rz(−ζ)
	return mul33(
		R3(-z),
		mul33(
			Ry(theta),
			R3(-zeta),
		),
	)
}

// ---- Nutation (IAU 2000B, 77 项) ----

// 基础自变量（角秒）
func fundArgs(T float64) (l, lp, F, D, Om float64) {
	l = (485868.249036 + 1717915923.2178*T + 31.8792*T*T + 0.051635*T*T*T - 0.00024470*T*T*T*T) * as2r
	lp = (1287104.79305 + 129596581.0481*T - 0.5532*T*T + 0.000136*T*T*T - 0.00001149*T*T*T*T) * as2r
	F = (335779.526232 + 1739527262.8478*T - 12.7512*T*T - 0.001037*T*T*T + 0.00000417*T*T*T*T) * as2r
	D = (1072260.70369 + 1602961601.2090*T - 6.3706*T*T + 0.006593*T*T*T - 0.00003169*T*T*T*T) * as2r
	Om = (450160.398036 - 6962890.5431*T + 7.4722*T*T + 0.007702*T*T*T - 0.00005939*T*T*T*T) * as2r
	return
}

// nut00bTerm 单个章动项的系数
// 系数来自 SOFA iau_nut00b，存储在内部整数单位 (×1e-8 得到角秒)
type nut00bTerm struct {
	nl, nlp, nF, nD, nOm int     // 自变量乘数
	sp, spt             float64 // Δψ 系数
	cp, cpt             float64 // Δε 系数
}

func initNut00b() []nut00bTerm {
	return []nut00bTerm{
		{0, 0, 0, 0, 1, -172064161, -174666, 92052331, 9086},
		{0, 0, 2, -2, 2, -13170906, -1675, 5730336, -3015},
		{0, 0, 2, 0, 2, -2276413, -234, 978459, -485},
		{0, 0, 0, 0, 2, 2074554, 207, -897492, 470},
		{0, 1, 0, 0, 0, 1475877, -3633, 73871, -184},
		{0, 1, 2, -2, 2, -516821, 1226, 224386, -677},
		{1, 0, 0, 0, 0, 711159, 73, -6750, 0},
		{0, 0, 2, 0, 1, -387298, -367, 200728, 18},
		{1, 0, 2, 0, 2, -301461, -36, 129025, -63},
		{0, -1, 2, -2, 2, 215829, -494, -95929, 299},
		{0, 0, 2, -2, 1, 128227, 137, -68982, -9},
		{-1, 0, 2, 0, 2, 123457, 11, -53311, 32},
		{-1, 0, 0, 2, 0, 156994, 10, -1235, 0},
		{1, 0, 0, 0, 1, 63110, 63, -33228, 0},
		{-1, 0, 0, 0, 1, -57976, -63, 31429, 0},
		{-1, 0, 2, 2, 2, -59641, -11, 25543, -11},
		{1, 0, 2, 0, 1, -51613, -42, 26366, 0},
		{-2, 0, 2, 0, 1, 45893, 50, -24236, -10},
		{0, 0, 0, 2, 0, 63384, 11, -1220, 0},
		{0, 0, 2, 2, 2, -38571, -1, 16452, -11},
		{0, -2, 2, -2, 2, 32481, 0, -13870, 0},
		{-2, 0, 0, 2, 0, -47722, 0, 477, 0},
		{2, 0, 2, 0, 2, -31046, -1, 13238, -11},
		{1, 0, 2, -2, 2, 28593, 0, -12338, 10},
		{-1, 0, 2, 0, 1, 20441, 21, -10758, 0},
		{2, 0, 0, 0, 0, 29243, 0, -609, 0},
		{0, 0, 2, 0, 0, 25887, 0, -550, 0},
		{0, 1, 0, 0, 1, -14053, -25, 7571, -2},
		{-1, 0, 0, 2, 1, 15164, 10, -8005, 0},
		{0, 2, 2, -2, 2, -15794, 72, 6850, -42},
		{0, 0, -2, 2, 0, 21783, 0, 167, 0},
		{1, 0, 0, -2, 1, -12873, -10, 6953, 0},
		{0, -1, 0, 0, 1, -12654, 11, 6415, 0},
		{-1, 0, 2, 2, 1, -10204, 0, 5222, 0},
		{0, 2, 0, 0, 0, 16707, -85, 168, -1},
		{1, 0, 2, 2, 2, -7691, 0, 3268, 0},
		{-2, 0, 2, 0, 0, -11024, 0, 104, 0},
		{0, 1, 2, 0, 2, 7566, -21, -3250, 0},
		{0, 0, 2, 2, 1, -6637, -11, 3353, 0},
		{0, -1, 2, 0, 2, -7141, 21, 3070, 0},
		{0, 0, 0, 2, 1, -6302, -11, 3272, 0},
		{1, 0, 2, -2, 1, 5800, 10, -3045, 0},
		{2, 0, 2, -2, 2, 6443, 0, -2768, 0},
		{-2, 0, 0, 2, 1, -5774, -11, 3041, 0},
		{2, 0, 2, 0, 1, -5350, 0, 2695, 0},
		{0, -1, 2, -2, 1, -4752, -11, 2719, 0},
		{0, 0, 0, -2, 1, -4940, -11, 2720, 0},
		{-1, -1, 0, 2, 0, 7350, 0, -51, 0},
		{2, 0, 0, -2, 1, 4065, 0, -2206, 0},
		{1, 0, 0, 2, 0, 6579, 0, -199, 0},
		{0, 1, 2, -2, 1, 3579, 0, -1900, 0},
		{1, -1, 0, 0, 0, 4725, 0, -41, 0},
		{-2, 0, 2, 0, 2, -3075, 0, 1313, 0},
		{3, 0, 2, 0, 2, -2904, 0, 1233, 0},
		{0, -1, 0, 2, 0, 4348, 0, -81, 0},
		{1, -1, 2, 0, 2, -2878, 0, 1232, 0},
		{0, 0, 0, 1, 0, -4230, 0, -20, 0},
		{-1, -1, 2, 2, 2, -2819, 0, 1207, 0},
		{-1, 0, 2, 0, 0, -4056, 0, 40, 0},
		{0, -1, 2, 2, 2, -2647, 0, 1129, 0},
		{-2, 0, 0, 0, 1, -2294, 0, 1266, 0},
		{1, 1, 2, 0, 2, 2481, 0, -1062, 0},
		{2, 0, 0, 0, 1, 2179, 0, -1129, 0},
		{-1, 1, 0, 1, 0, 3276, 0, -9, 0},
		{1, 1, 0, 0, 0, -3389, 0, 35, 0},
		{1, 0, 2, 0, 0, 3339, 0, -107, 0},
		{-1, 0, 2, -2, 1, -1987, 0, 1073, 0},
		{1, 0, 0, 0, 2, -1981, 0, 854, 0},
		{-1, 0, 0, 1, 0, 4026, 0, -55, 0},
		{0, 0, 2, 1, 2, 1660, 0, -710, 0},
		{-1, 0, 2, 4, 2, -1521, 0, 647, 0},
		{-1, 1, 0, 1, 1, 1314, 0, -700, 0},
		{0, -2, 2, -2, 1, -1283, 0, 672, 0},
		{1, 0, 2, 2, 1, -1331, 0, 663, 0},
		{-2, 0, 2, 2, 2, 1383, 0, -594, 0},
		{-1, 0, 0, 0, 2, 1405, 0, -610, 0},
		{1, 1, 2, -2, 2, 1290, 0, -556, 0},
	}
}

// nutation 计算章动 (Δψ, Δε)，单位 rad
func nutation(T float64) (dpsi, deps float64) {
	l, lp, F, D, Om := fundArgs(T)
	terms := initNut00b()

	dpsiSum := 0.0
	depsSum := 0.0

	for _, t := range terms {
		arg := float64(t.nl)*l + float64(t.nlp)*lp + float64(t.nF)*F +
			float64(t.nD)*D + float64(t.nOm)*Om

		// Δψ: (sp + spt × T) × sin(arg) × 0.1 mas
		dpsiSum += (t.sp + t.spt*T) * math.Sin(arg)
		// Δε: (cp + cpt × T) × cos(arg) × 0.1 mas
		depsSum += (t.cp + t.cpt*T) * math.Cos(arg)
	}

	// 系数存储在 0.1 µas × 1000 单位 → 乘以 1e-8 得到角秒
	dpsi = dpsiSum * 1e-8 * as2r
	deps = depsSum * 1e-8 * as2r
	return
}

// meanObliquity IAU 2006 平黄赤交角 (rad)
func meanObliquity(T float64) float64 {
	const eps0 = 84381.406 // "
	eps := (eps0 +
		(-46.836769)*T +
		(-0.0001831)*T*T +
		(0.00200340)*T*T*T +
		(-5.76e-7)*T*T*T*T +
		(-4.34e-8)*T*T*T*T*T) * as2r
	return eps
}

// nutationMatrix 章动矩阵 N(t)
// N = Rx(-ε) · Rz(-Δψ) · Rx(ε_mean)
func nutationMatrix(T float64) [3][3]float64 {
	dpsi, deps := nutation(T)
	epsMean := meanObliquity(T)
	epsTrue := epsMean + deps

	sep, cep := math.Sin(epsTrue), math.Cos(epsTrue)
	sd, cd := math.Sin(dpsi), math.Cos(dpsi)
	sem, cem := math.Sin(epsMean), math.Cos(epsMean)

	// N = Rx(-ε_true) · Rz(-Δψ) · Rx(ε_mean)
	return [3][3]float64{
		{cd, -sd * cem, -sd * sem},
		{sd * cep, cd*cep*cem + sep*sem, cd*cep*sem - sep*cem},
		{sd * sep, cd*sep*cem - cep*sem, cd*sep*sem + cep*cem},
	}
}

// ---- GAST (格林威治视恒星时) ----
// GAST = GMST + 赤经章动 (equation of equinoxes)
// equation of equinoxes = Δψ · cos(ε_true)
func GAST(jd float64) float64 {
	T := (jd - 2451545.0) / 36525.0
	gmst := greenwichsrt(jd)
	dpsi, deps := nutation(T)
	epsTrue := meanObliquity(T) + deps
	eeq := dpsi * math.Cos(epsTrue) // 赤经章动 (分点差)
	return math.Mod(gmst+eeq, tau)
}

// ---- 完整变换链 ----
// ECI(GCRF) → ECEF(ITRF): r_ecef = Rz(GAST) × N × P × B × r_eci
// ECEF(ITRF) → ECI(GCRF): r_eci = B^T × P^T × N^T × Rz(-GAST) × r_ecef

// GCRF2ITRF 计算完整 GCRF → ITRF 旋转矩阵（不含极移）
// M = Rz(GAST) · N · P · B
func GCRF2ITRF(jd float64) [3][3]float64 {
	T := (jd - 2451545.0) / 36525.0
	B := frameBiasMatrix()
	P := precessionMatrix(T)
	N := nutationMatrix(T)
	gast := GAST(jd)

	return mul33(R3(gast), mul33(N, mul33(P, B)))
}

// ---- 辅助矩阵运算 ----

// Rx 绕 X 轴旋转 angle 弧度的旋转矩阵
func Rx(angle float64) [3][3]float64 {
	c, s := math.Cos(angle), math.Sin(angle)
	return [3][3]float64{
		{1, 0, 0},
		{0, c, s},
		{0, -s, c},
	}
}

// Ry 绕 Y 轴旋转 angle 弧度的旋转矩阵
func Ry(angle float64) [3][3]float64 {
	c, s := math.Cos(angle), math.Sin(angle)
	return [3][3]float64{
		{c, 0, -s},
		{0, 1, 0},
		{s, 0, c},
	}
}

// mul33 两个 3×3 矩阵相乘
func mul33(a, b [3][3]float64) [3][3]float64 {
	return [3][3]float64{
		{
			a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0],
			a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1],
			a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2],
		},
		{
			a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0],
			a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1],
			a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2],
		},
		{
			a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0],
			a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1],
			a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2],
		},
	}
}
