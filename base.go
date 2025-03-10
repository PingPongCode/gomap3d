package gomap3d

import (
	"math"
	"time"
)

// gomap3d基础算子集合

// ENU2AER 将ENU坐标转换为方位角、仰角和斜距
func ENU2AER(e, n, u float64) (az, el, srange float64) {
	if math.Abs(e) < 1e-3 {
		e = 0
	}
	if math.Abs(n) < 1e-3 {
		n = 0
	}
	if math.Abs(u) < 1e-3 {
		u = 0
	}
	r := math.Hypot(e, n)
	srange = math.Hypot(r, u)
	el = math.Atan2(u, r) * 180 / math.Pi
	az = math.Mod(math.Atan2(e, n)*180/math.Pi+360, 360)
	return
}

// AER2ENU 将方位角、仰角和斜距转换为ENU坐标
func AER2ENU(az, el, srange float64) (e, n, u float64) {
	azRad := az * math.Pi / 180
	elRad := el * math.Pi / 180

	r := srange * math.Cos(elRad)
	e = r * math.Sin(azRad)
	n = r * math.Cos(azRad)
	u = srange * math.Sin(elRad)
	return
}

// Geodetic2ECEF 将地理坐标转换为ECEF坐标
func Geodetic2ECEF(lat, lon, alt float64, ell *Ellipsoid) (x, y, z float64) {
	latRad := lat * math.Pi / 180
	lonRad := lon * math.Pi / 180
	n := math.Pow(ell.SemimajorAxis, 2) / math.Hypot(ell.SemimajorAxis*math.Cos(latRad), ell.SemiminorAxis*math.Sin(latRad))
	x = (n + alt) * math.Cos(latRad) * math.Cos(lonRad)
	y = (n + alt) * math.Cos(latRad) * math.Sin(lonRad)
	z = (n*math.Pow((ell.SemiminorAxis/ell.SemimajorAxis), 2) + alt) * math.Sin(latRad)
	return
}

// ECEF2Geodetic 将ECEF坐标转换为地理坐标
func ECEF2Geodetic(x, y, z float64, ell *Ellipsoid) (lat, lon, alt float64) {
	// 实现You, Rey-Jer. (2000)算法
	// 此处简化为迭代法实现
	a := ell.SemimajorAxis
	b := ell.SemiminorAxis
	e2 := (a*a - b*b) / (a * a)
	eps := 1e-12

	p := math.Sqrt(x*x + y*y)
	lon = math.Atan2(y, x)

	// 初始估计
	lat = math.Atan2(z, p*(1-e2))

	for {
		n := a / math.Sqrt(1-e2*math.Pow(math.Sin(lat), 2))
		h := p/math.Cos(lat) - n

		latNew := math.Atan(z / p * 1 / (1 - e2*n/(n+h)))

		if math.Abs(latNew-lat) < eps {
			lat = latNew
			alt = h
			break
		}
		lat = latNew
	}

	lat = lat * 180 / math.Pi
	lon = lon * 180 / math.Pi
	return
}

// ECEF2ENU 将ECEF坐标转换为ENU坐标
func ECEF2ENU(x, y, z, lat0, lon0, h0 float64, ell *Ellipsoid) (e, n, u float64) {
	// 转换为本地ENU坐标
	x0, y0, z0 := Geodetic2ECEF(lat0, lon0, h0, ell)
	dx := x - x0
	dy := y - y0
	dz := z - z0

	lat0Rad := lat0 * math.Pi / 180
	lon0Rad := lon0 * math.Pi / 180

	// 旋转矩阵
	e = -math.Sin(lon0Rad)*dx + math.Cos(lon0Rad)*dy
	n = -math.Sin(lat0Rad)*math.Cos(lon0Rad)*dx -
		math.Sin(lat0Rad)*math.Sin(lon0Rad)*dy +
		math.Cos(lat0Rad)*dz
	u = math.Cos(lat0Rad)*math.Cos(lon0Rad)*dx +
		math.Cos(lat0Rad)*math.Sin(lon0Rad)*dy +
		math.Sin(lat0Rad)*dz
	return
}

// ENU2ECEF 将ENU坐标转换为ECEF坐标
func ENU2ECEF(e, n, u, lat0, lon0, h0 float64, ell *Ellipsoid) (x, y, z float64) {
	lat0Rad := lat0 * math.Pi / 180
	lon0Rad := lon0 * math.Pi / 180

	// 旋转矩阵转置
	dx := -math.Sin(lon0Rad)*e -
		math.Sin(lat0Rad)*math.Cos(lon0Rad)*n +
		math.Cos(lat0Rad)*math.Cos(lon0Rad)*u
	dy := math.Cos(lon0Rad)*e -
		math.Sin(lat0Rad)*math.Sin(lon0Rad)*n +
		math.Cos(lat0Rad)*math.Sin(lon0Rad)*u
	dz := math.Cos(lat0Rad)*n +
		math.Sin(lat0Rad)*u

	x0, y0, z0 := Geodetic2ECEF(lat0, lon0, h0, ell)
	return x0 + dx, y0 + dy, z0 + dz
}

const tau = 2 * math.Pi

// 天文计算相关函数
// juliandate 计算给定时间的儒略日
func juliandate(t time.Time) float64 {
	year := t.Year()
	month := t.Month()
	// 处理月份调整
	if month < time.March {
		year--
		month += 12
	}
	A := int(year / 100)
	B := 2 - A + int(A/4)
	C := ((t.Second()+t.Nanosecond()/1e9)/60 + t.Hour()) / 24
	result := float64(int(365.25*float64(year+4716)+float64(int(30.6001*float64(month+1))))) + float64(t.Day()) + float64(B) - 1524.5 + float64(C)
	return result

}

// greenwichsrt 计算格林威治恒星时（弧度）
func greenwichsrt(jd float64) float64 {
	tUT1 := (jd - 2451545.0) / 36525.0

	gmstSec := 67310.54841 +
		(876600*3600+8640184.812866)*tUT1 +
		0.093104*math.Pow(tUT1, 2) -
		6.2e-6*math.Pow(tUT1, 3)

	gmstRad := gmstSec * tau / 86400.0
	return math.Mod(gmstRad, tau)
}

// 生成绕Z轴旋转x弧度的矩阵R3
func R3(angle float64) [3][3]float64 {
	cos := math.Cos(angle)
	sin := math.Sin(angle)
	return [3][3]float64{
		{cos, sin, 0},
		{-sin, cos, 0},
		{0, 0, 1},
	}
}

// 矩阵转置
func transpose(m [3][3]float64) [3][3]float64 {
	return [3][3]float64{
		{m[0][0], m[1][0], m[2][0]},
		{m[0][1], m[1][1], m[2][1]},
		{m[0][2], m[1][2], m[2][2]},
	}
}

// 乘法运算，矩阵和向量
func multiplyMatrixVector(m [3][3]float64, v [3]float64) [3]float64 {
	return [3]float64{
		m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2],
		m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2],
		m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2],
	}
}

// ECI2ECEF 将ECI坐标转换为ECEF坐标
func ECI2ECEF(x, y, z float64, t time.Time) (xEcef, yEcef, zEcef float64) {
	jd := juliandate(t)
	gst := greenwichsrt(jd)
	rotMat := R3(gst)
	eciVec := [3]float64{x, y, z}
	ecefVec := multiplyMatrixVector(rotMat, eciVec)
	return ecefVec[0], ecefVec[1], ecefVec[2]
}

// ECEF2ECI 将ECEF坐标转换为ECI坐标
func ECEF2ECI(x, y, z float64, t time.Time) (xEci, yEci, zEci float64) {
	jd := juliandate(t)
	gst := greenwichsrt(jd)
	rotMat := R3(gst)
	transposed := transpose(rotMat)
	ecefVec := [3]float64{x, y, z}
	eciVec := multiplyMatrixVector(transposed, ecefVec)
	return eciVec[0], eciVec[1], eciVec[2]
}
