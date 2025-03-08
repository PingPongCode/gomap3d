package gomap3d

import (
	"math"
	"time"
)

// gomap3d基础算子集合

// ENU2AER 将ENU坐标转换为方位角、仰角和斜距
func ENU2AER(e, n, u float64) (az, el, srange float64) {
	r := math.Hypot(e, n)
	srange = math.Hypot(r, u)
	el = math.Atan2(u, r) / math.Pi * 180
	tau := 2 * math.Pi
	az = math.Mod(math.Atan2(e, n), tau) / math.Pi * 180
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

// 天文计算相关函数
// juliandate 计算给定时间的儒略日
func juliandate(t time.Time) float64 {
	j2000 := time.Date(2000, 1, 1, 12, 0, 0, 0, time.UTC)
	duration := t.Sub(j2000)
	days := duration.Seconds() / 86400.0
	return 2451545.0 + days
}

// greenwichsrt 计算格林威治恒星时（弧度）
func greenwichsrt(jd float64) float64 {
	T := (jd - 2451545.0) / 36525.0
	theta := 280.46061837 + 360.98564736629*(jd-2451545.0) + 0.000387933*T*T - (T*T*T)/38710000.0
	theta = math.Mod(theta, 360.0)
	if theta < 0 {
		theta += 360.0
	}
	return theta * math.Pi / 180.0
}

// rotationMatrix3 生成绕Z轴旋转x弧度的矩阵
func rotationMatrix3(x float64) [3][3]float64 {
	cosX := math.Cos(x)
	sinX := math.Sin(x)
	return [3][3]float64{
		{cosX, sinX, 0},
		{-sinX, cosX, 0},
		{0, 0, 1},
	}
}

// multiplyMatrixVector 矩阵乘以向量
func multiplyMatrixVector(matrix [3][3]float64, vector [3]float64) [3]float64 {
	x := matrix[0][0]*vector[0] + matrix[0][1]*vector[1] + matrix[0][2]*vector[2]
	y := matrix[1][0]*vector[0] + matrix[1][1]*vector[1] + matrix[1][2]*vector[2]
	z := matrix[2][0]*vector[0] + matrix[2][1]*vector[1] + matrix[2][2]*vector[2]
	return [3]float64{x, y, z}
}

// ECI2ECEF 将ECI坐标转换为ECEF坐标
func ECI2ECEF(x, y, z float64, t time.Time) (xEcef, yEcef, zEcef float64) {
	jd := juliandate(t)
	gst := greenwichsrt(jd)
	matrix := rotationMatrix3(gst)
	vec := [3]float64{x, y, z}
	result := multiplyMatrixVector(matrix, vec)
	return result[0], result[1], result[2]
}

// ECEF2ECI 将ECEF坐标转换为ECI坐标
func ECEF2ECI(x, y, z float64, t time.Time) (xEci, yEci, zEci float64) {
	jd := juliandate(t)
	gst := greenwichsrt(jd)
	matrix := rotationMatrix3(-gst) // 使用逆旋转
	vec := [3]float64{x, y, z}
	result := multiplyMatrixVector(matrix, vec)
	return result[0], result[1], result[2]
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
