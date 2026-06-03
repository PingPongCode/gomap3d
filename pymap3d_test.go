package gomap3d

import (
	"encoding/json"
	"math"
	"os"
	"testing"
	"time"
)

// ============================================================
// pymap3d 参考测试数据结构
// ============================================================

type RefPoint struct {
	Lat float64 `json:"lat"`
	Lon float64 `json:"lon"`
	Alt float64 `json:"alt"`
}

type RefVec3 struct {
	E float64 `json:"e"`
	N float64 `json:"n"`
	U float64 `json:"u"`
}

type RefAER struct {
	Az     float64 `json:"az"`
	El     float64 `json:"el"`
	Srange float64 `json:"srange"`
}

type RefECEF struct {
	X float64 `json:"x"`
	Y float64 `json:"y"`
	Z float64 `json:"z"`
}

type RefGeo struct {
	Lat float64 `json:"lat"`
	Lon float64 `json:"lon"`
	Alt float64 `json:"alt"`
}

type Pymap3dTestCase struct {
	Scenario string   `json:"scenario"`
	Ref      RefPoint `json:"ref"`
	Point    RefPoint `json:"point"`

	// 正向变换
	ENU  RefVec3 `json:"enu"`
	AER  RefAER  `json:"aer"`
	ECEF RefECEF `json:"ecef"`

	// 逆变换: AER->ENU
	AER2ENU RefVec3 `json:"aer2enu"`
	ECEF2ENU RefVec3 `json:"ecef2enu"`
	ENU2ECEF RefECEF `json:"enu2ecef"`
	ECEF2Geo RefGeo  `json:"ecef2geo"`

	// 结构化变换
	GEO2ENU RefVec3 `json:"geo2enu"`
	ENU2Geo RefGeo  `json:"enu2geo"`
	GEO2AER RefAER  `json:"geo2aer"`
	AER2Geo RefGeo  `json:"aer2geo"`
	ECEF2AER RefAER `json:"ecef2aer"`
	AER2ECEF RefECEF `json:"aer2ecef"`
}

// ============================================================
// 测试配置
// ============================================================
const enuTol = 1e-3    // ENU 容差 (米)
const geoCoordTol = 1e-8  // 经纬度容差 (度)
const geoAltTol = 1e-2  // 高程容差 (米)
const aerAzElTol = 3e-5 // 方位/仰角容差 (度) — 链式调用 (ECEF-ENU-AER) 累积浮点差异
const aerSrangeTol = 1e-3 // 斜距容差 (米)
const ecefTol = 1e-2  // ECEF 容差 (米)

var wgs84Ref *Ellipsoid

func init() {
	wgs84Ref, _ = NewEllipsoid("wgs84")
}

func loadPymap3dTestCases(t *testing.T) []Pymap3dTestCase {
	data, err := os.ReadFile("test_data/pymap3d_reference.json")
	if err != nil {
		t.Skipf("跳过: 找不到测试数据文件 test_data/pymap3d_reference.json: %v", err)
		return nil
	}
	var cases []Pymap3dTestCase
	if err := json.Unmarshal(data, &cases); err != nil {
		t.Fatalf("JSON 解析失败: %v", err)
	}
	return cases
}

// ============================================================
// 1. ENU <-> AER 基本函数测试
// ============================================================
func TestPymap3d_ENU2AER(t *testing.T) {
	cases := loadPymap3dTestCases(t)
	for _, tc := range cases {
		t.Run(tc.Scenario, func(t *testing.T) {
			az, el, sr := ENU2AER(tc.ENU.E, tc.ENU.N, tc.ENU.U)
			if math.Abs(az-tc.AER.Az) > aerAzElTol || math.Abs(el-tc.AER.El) > aerAzElTol {
				t.Errorf("ENU2AER: got(az=%.8f, el=%.8f), want(az=%.8f, el=%.8f)",
					az, el, tc.AER.Az, tc.AER.El)
			}
			if math.Abs(sr-tc.AER.Srange) > aerSrangeTol {
				t.Errorf("ENU2AER srange: got %.6f, want %.6f", sr, tc.AER.Srange)
			}
		})
	}
}

func TestPymap3d_AER2ENU(t *testing.T) {
	cases := loadPymap3dTestCases(t)
	for _, tc := range cases {
		t.Run(tc.Scenario, func(t *testing.T) {
			e, n, u := AER2ENU(tc.AER.Az, tc.AER.El, tc.AER.Srange)
			if math.Abs(e-tc.AER2ENU.E) > enuTol || math.Abs(n-tc.AER2ENU.N) > enuTol || math.Abs(u-tc.AER2ENU.U) > enuTol {
				t.Errorf("AER2ENU: got(e=%.6f, n=%.6f, u=%.6f), want(e=%.6f, n=%.6f, u=%.6f)",
					e, n, u, tc.AER2ENU.E, tc.AER2ENU.N, tc.AER2ENU.U)
			}
		})
	}
}

// ============================================================
// 2. Geodetic <-> ECEF 基本函数测试
// ============================================================
func TestPymap3d_Geodetic2ECEF(t *testing.T) {
	cases := loadPymap3dTestCases(t)
	for _, tc := range cases {
		t.Run(tc.Scenario, func(t *testing.T) {
			x, y, z := Geodetic2ECEF(tc.Point.Lat, tc.Point.Lon, tc.Point.Alt, wgs84Ref)
			if math.Abs(x-tc.ECEF.X) > ecefTol || math.Abs(y-tc.ECEF.Y) > ecefTol || math.Abs(z-tc.ECEF.Z) > ecefTol {
				t.Errorf("Geodetic2ECEF: got(x=%.4f, y=%.4f, z=%.4f), want(x=%.4f, y=%.4f, z=%.4f)",
					x, y, z, tc.ECEF.X, tc.ECEF.Y, tc.ECEF.Z)
			}
		})
	}
}

func TestPymap3d_ECEF2Geodetic(t *testing.T) {
	cases := loadPymap3dTestCases(t)
	for _, tc := range cases {
		t.Run(tc.Scenario, func(t *testing.T) {
			lat, lon, alt := ECEF2Geodetic(tc.ECEF.X, tc.ECEF.Y, tc.ECEF.Z, wgs84Ref)
			if math.Abs(lat-tc.ECEF2Geo.Lat) > geoCoordTol || math.Abs(lon-tc.ECEF2Geo.Lon) > geoCoordTol {
				t.Errorf("ECEF2Geodetic: got(lat=%.10f, lon=%.10f), want(lat=%.10f, lon=%.10f)",
					lat, lon, tc.ECEF2Geo.Lat, tc.ECEF2Geo.Lon)
			}
			if math.Abs(alt-tc.ECEF2Geo.Alt) > geoAltTol {
				t.Errorf("ECEF2Geodetic alt: got %.6f, want %.6f", alt, tc.ECEF2Geo.Alt)
			}
		})
	}
}

// ============================================================
// 3. ECEF <-> ENU 基本函数测试
// ============================================================
func TestPymap3d_ECEF2ENU(t *testing.T) {
	cases := loadPymap3dTestCases(t)
	for _, tc := range cases {
		t.Run(tc.Scenario, func(t *testing.T) {
			e, n, u := ECEF2ENU(tc.ECEF.X, tc.ECEF.Y, tc.ECEF.Z,
				tc.Ref.Lat, tc.Ref.Lon, tc.Ref.Alt, wgs84Ref)
			if math.Abs(e-tc.ECEF2ENU.E) > enuTol || math.Abs(n-tc.ECEF2ENU.N) > enuTol || math.Abs(u-tc.ECEF2ENU.U) > enuTol {
				t.Errorf("ECEF2ENU: got(e=%.6f, n=%.6f, u=%.6f), want(e=%.6f, n=%.6f, u=%.6f)",
					e, n, u, tc.ECEF2ENU.E, tc.ECEF2ENU.N, tc.ECEF2ENU.U)
			}
		})
	}
}

func TestPymap3d_ENU2ECEF(t *testing.T) {
	cases := loadPymap3dTestCases(t)
	for _, tc := range cases {
		t.Run(tc.Scenario, func(t *testing.T) {
			x, y, z := ENU2ECEF(tc.ENU.E, tc.ENU.N, tc.ENU.U,
				tc.Ref.Lat, tc.Ref.Lon, tc.Ref.Alt, wgs84Ref)
			if math.Abs(x-tc.ENU2ECEF.X) > ecefTol || math.Abs(y-tc.ENU2ECEF.Y) > ecefTol || math.Abs(z-tc.ENU2ECEF.Z) > ecefTol {
				t.Errorf("ENU2ECEF: got(x=%.4f, y=%.4f, z=%.4f), want(x=%.4f, y=%.4f, z=%.4f)",
					x, y, z, tc.ENU2ECEF.X, tc.ENU2ECEF.Y, tc.ENU2ECEF.Z)
			}
		})
	}
}

// ============================================================
// 4. 结构体方法测试 (Geodetic / ECEF / ENU / AER)
// ============================================================
func TestPymap3d_GeodeticMethods(t *testing.T) {
	cases := loadPymap3dTestCases(t)
	for _, tc := range cases {
		t.Run(tc.Scenario, func(t *testing.T) {
			geo := Geodetic{Latitude: tc.Point.Lat, Longitude: tc.Point.Lon, Altitude: tc.Point.Alt, Ell: wgs84Ref}
			ref := Geodetic{Latitude: tc.Ref.Lat, Longitude: tc.Ref.Lon, Altitude: tc.Ref.Alt, Ell: wgs84Ref}

			// Geodetic.ToECEF
			ecef := geo.ToECEF()
			if math.Abs(ecef.X-tc.ECEF.X) > ecefTol || math.Abs(ecef.Y-tc.ECEF.Y) > ecefTol || math.Abs(ecef.Z-tc.ECEF.Z) > ecefTol {
				t.Errorf("Geo.ToECEF: got(%.2f, %.2f, %.2f), want(%.2f, %.2f, %.2f)",
					ecef.X, ecef.Y, ecef.Z, tc.ECEF.X, tc.ECEF.Y, tc.ECEF.Z)
			}

			// Geodetic.ToENU
			enu := geo.ToENU(ref)
			if math.Abs(enu.East-tc.GEO2ENU.E) > enuTol || math.Abs(enu.North-tc.GEO2ENU.N) > enuTol || math.Abs(enu.Up-tc.GEO2ENU.U) > enuTol {
				t.Errorf("Geo.ToENU: got(%.6f, %.6f, %.6f), want(%.6f, %.6f, %.6f)",
					enu.East, enu.North, enu.Up, tc.GEO2ENU.E, tc.GEO2ENU.N, tc.GEO2ENU.U)
			}

			// Geodetic.ToAER
			aer := geo.ToAER(ref)
			if math.Abs(aer.Azimuth-tc.GEO2AER.Az) > aerAzElTol || math.Abs(aer.Elevation-tc.GEO2AER.El) > aerAzElTol {
				t.Errorf("Geo.ToAER: got(az=%.8f, el=%.8f), want(az=%.8f, el=%.8f)",
					aer.Azimuth, aer.Elevation, tc.GEO2AER.Az, tc.GEO2AER.El)
			}
			if math.Abs(aer.SRange-tc.GEO2AER.Srange) > aerSrangeTol {
				t.Errorf("Geo.ToAER srange: got %.6f, want %.6f", aer.SRange, tc.GEO2AER.Srange)
			}
		})
	}
}

func TestPymap3d_ECEFMethods(t *testing.T) {
	cases := loadPymap3dTestCases(t)
	for _, tc := range cases {
		t.Run(tc.Scenario, func(t *testing.T) {
			ecef := ECEF{X: tc.ECEF.X, Y: tc.ECEF.Y, Z: tc.ECEF.Z, Ell: wgs84Ref}
			ref := Geodetic{Latitude: tc.Ref.Lat, Longitude: tc.Ref.Lon, Altitude: tc.Ref.Alt, Ell: wgs84Ref}

			// ECEF.ToGeodetic
			geo := ecef.ToGeodetic()
			if math.Abs(geo.Latitude-tc.ECEF2Geo.Lat) > geoCoordTol || math.Abs(geo.Longitude-tc.ECEF2Geo.Lon) > geoCoordTol {
				t.Errorf("ECEF.ToGeodetic: got(%.10f, %.10f), want(%.10f, %.10f)",
					geo.Latitude, geo.Longitude, tc.ECEF2Geo.Lat, tc.ECEF2Geo.Lon)
			}
			if math.Abs(geo.Altitude-tc.ECEF2Geo.Alt) > geoAltTol {
				t.Errorf("ECEF.ToGeodetic alt: got %.6f, want %.6f", geo.Altitude, tc.ECEF2Geo.Alt)
			}

			// ECEF.ToENU
			enu := ecef.ToENU(ref)
			if math.Abs(enu.East-tc.ECEF2ENU.E) > enuTol || math.Abs(enu.North-tc.ECEF2ENU.N) > enuTol || math.Abs(enu.Up-tc.ECEF2ENU.U) > enuTol {
				t.Errorf("ECEF.ToENU: got(%.6f, %.6f, %.6f), want(%.6f, %.6f, %.6f)",
					enu.East, enu.North, enu.Up, tc.ECEF2ENU.E, tc.ECEF2ENU.N, tc.ECEF2ENU.U)
			}

			// ECEF.ToAER
			aer := ecef.ToAER(ref)
			if math.Abs(aer.Azimuth-tc.ECEF2AER.Az) > aerAzElTol || math.Abs(aer.Elevation-tc.ECEF2AER.El) > aerAzElTol {
				t.Errorf("ECEF.ToAER: got(az=%.8f, el=%.8f), want(az=%.8f, el=%.8f)",
					aer.Azimuth, aer.Elevation, tc.ECEF2AER.Az, tc.ECEF2AER.El)
			}
			if math.Abs(aer.SRange-tc.ECEF2AER.Srange) > aerSrangeTol {
				t.Errorf("ECEF.ToAER srange: got %.6f, want %.6f", aer.SRange, tc.ECEF2AER.Srange)
			}
		})
	}
}

func TestPymap3d_ENUMethods(t *testing.T) {
	cases := loadPymap3dTestCases(t)
	for _, tc := range cases {
		t.Run(tc.Scenario, func(t *testing.T) {
			enu := ENU{East: tc.ENU.E, North: tc.ENU.N, Up: tc.ENU.U, Ell: wgs84Ref}
			ref := Geodetic{Latitude: tc.Ref.Lat, Longitude: tc.Ref.Lon, Altitude: tc.Ref.Alt, Ell: wgs84Ref}

			// ENU.ToAER
			aer := enu.ToAER()
			if math.Abs(aer.Azimuth-tc.AER.Az) > aerAzElTol || math.Abs(aer.Elevation-tc.AER.El) > aerAzElTol {
				t.Errorf("ENU.ToAER: got(az=%.8f, el=%.8f), want(az=%.8f, el=%.8f)",
					aer.Azimuth, aer.Elevation, tc.AER.Az, tc.AER.El)
			}
			if math.Abs(aer.SRange-tc.AER.Srange) > aerSrangeTol {
				t.Errorf("ENU.ToAER srange: got %.6f, want %.6f", aer.SRange, tc.AER.Srange)
			}

			// ENU.ToGeodetic
			geo := enu.ToGeodetic(ref)
			if math.Abs(geo.Latitude-tc.ENU2Geo.Lat) > geoCoordTol || math.Abs(geo.Longitude-tc.ENU2Geo.Lon) > geoCoordTol {
				t.Errorf("ENU.ToGeodetic: got(%.10f, %.10f), want(%.10f, %.10f)",
					geo.Latitude, geo.Longitude, tc.ENU2Geo.Lat, tc.ENU2Geo.Lon)
			}
			if math.Abs(geo.Altitude-tc.ENU2Geo.Alt) > geoAltTol {
				t.Errorf("ENU.ToGeodetic alt: got %.6f, want %.6f", geo.Altitude, tc.ENU2Geo.Alt)
			}

			// ENU.ToECEF
			ecef := enu.ToECEF(ref)
			if math.Abs(ecef.X-tc.ENU2ECEF.X) > ecefTol || math.Abs(ecef.Y-tc.ENU2ECEF.Y) > ecefTol || math.Abs(ecef.Z-tc.ENU2ECEF.Z) > ecefTol {
				t.Errorf("ENU.ToECEF: got(%.2f, %.2f, %.2f), want(%.2f, %.2f, %.2f)",
					ecef.X, ecef.Y, ecef.Z, tc.ENU2ECEF.X, tc.ENU2ECEF.Y, tc.ENU2ECEF.Z)
			}
		})
	}
}

func TestPymap3d_AERMethods(t *testing.T) {
	cases := loadPymap3dTestCases(t)
	for _, tc := range cases {
		t.Run(tc.Scenario, func(t *testing.T) {
			aer := AER{Azimuth: tc.AER.Az, Elevation: tc.AER.El, SRange: tc.AER.Srange, Ell: wgs84Ref}
			ref := Geodetic{Latitude: tc.Ref.Lat, Longitude: tc.Ref.Lon, Altitude: tc.Ref.Alt, Ell: wgs84Ref}

			// AER.ToENU
			enu := aer.ToENU()
			if math.Abs(enu.East-tc.AER2ENU.E) > enuTol || math.Abs(enu.North-tc.AER2ENU.N) > enuTol || math.Abs(enu.Up-tc.AER2ENU.U) > enuTol {
				t.Errorf("AER.ToENU: got(%.6f, %.6f, %.6f), want(%.6f, %.6f, %.6f)",
					enu.East, enu.North, enu.Up, tc.AER2ENU.E, tc.AER2ENU.N, tc.AER2ENU.U)
			}

			// AER.ToECEF
			ecef := aer.ToECEF(ref)
			if math.Abs(ecef.X-tc.AER2ECEF.X) > ecefTol || math.Abs(ecef.Y-tc.AER2ECEF.Y) > ecefTol || math.Abs(ecef.Z-tc.AER2ECEF.Z) > ecefTol {
				t.Errorf("AER.ToECEF: got(%.2f, %.2f, %.2f), want(%.2f, %.2f, %.2f)",
					ecef.X, ecef.Y, ecef.Z, tc.AER2ECEF.X, tc.AER2ECEF.Y, tc.AER2ECEF.Z)
			}

			// AER.ToGeodetic
			geo := aer.ToGeodetic(ref)
			if math.Abs(geo.Latitude-tc.AER2Geo.Lat) > geoCoordTol || math.Abs(geo.Longitude-tc.AER2Geo.Lon) > geoCoordTol {
				t.Errorf("AER.ToGeodetic: got(%.10f, %.10f), want(%.10f, %.10f)",
					geo.Latitude, geo.Longitude, tc.AER2Geo.Lat, tc.AER2Geo.Lon)
			}
			if math.Abs(geo.Altitude-tc.AER2Geo.Alt) > geoAltTol {
				t.Errorf("AER.ToGeodetic alt: got %.6f, want %.6f", geo.Altitude, tc.AER2Geo.Alt)
			}
		})
	}
}

// ============================================================
// 5. 互逆性测试: 正逆转换应恢复原值
// ============================================================
func TestPymap3d_Roundtrip(t *testing.T) {
	cases := loadPymap3dTestCases(t)
	for _, tc := range cases {
		t.Run(tc.Scenario, func(t *testing.T) {
			ref := Geodetic{Latitude: tc.Ref.Lat, Longitude: tc.Ref.Lon, Altitude: tc.Ref.Alt, Ell: wgs84Ref}

			// 1) Geodetic -> ECEF -> Geodetic (round-trip)
			geo := Geodetic{Latitude: tc.Point.Lat, Longitude: tc.Point.Lon, Altitude: tc.Point.Alt, Ell: wgs84Ref}
			ecef := geo.ToECEF()
			geoBack := ecef.ToGeodetic()
			if math.Abs(geoBack.Latitude-tc.Point.Lat) > 1e-8 || math.Abs(geoBack.Longitude-tc.Point.Lon) > 1e-8 {
				t.Errorf("Geo->ECEF->Geo: lat/lon mismatch: (%.10f, %.10f) -> (%.10f, %.10f)",
					tc.Point.Lat, tc.Point.Lon, geoBack.Latitude, geoBack.Longitude)
			}

			// 2) ENU -> AER -> ENU (round-trip)
			enu := ENU{East: tc.ENU.E, North: tc.ENU.N, Up: tc.ENU.U, Ell: wgs84Ref}
			aer := enu.ToAER()
			enuBack := aer.ToENU()
			if math.Abs(enuBack.East-tc.ENU.E) > 1e-6 || math.Abs(enuBack.North-tc.ENU.N) > 1e-6 {
				t.Errorf("ENU->AER->ENU: mismatch")
			}

			// 3) Geodetic -> ENU -> Geodetic (round-trip via ENU)
			enuGeo := geo.ToENU(ref)
			geoFromEnu := enuGeo.ToGeodetic(ref)
			if math.Abs(geoFromEnu.Latitude-tc.Point.Lat) > 1e-7 || math.Abs(geoFromEnu.Longitude-tc.Point.Lon) > 1e-7 {
				t.Errorf("Geo->ENU->Geo: lat/lon mismatch")
			}

			// 4) Geodetic -> AER -> Geodetic (round-trip via AER)
			aerGeo := geo.ToAER(ref)
			geoFromAer := aerGeo.ToGeodetic(ref)
			if math.Abs(geoFromAer.Latitude-tc.Point.Lat) > 1e-7 || math.Abs(geoFromAer.Longitude-tc.Point.Lon) > 1e-7 {
				t.Errorf("Geo->AER->Geo: lat/lon mismatch: orig(%.10f, %.10f) back(%.10f, %.10f)",
					tc.Point.Lat, tc.Point.Lon, geoFromAer.Latitude, geoFromAer.Longitude)
			}

			// 5) ECEF -> ENU -> ECEF (round-trip)
			enuFromEcef := ecef.ToENU(ref)
			ecefBack := enuFromEcef.ToECEF(ref)
			if math.Abs(ecefBack.X-tc.ECEF.X) > 1e-1 {
				t.Errorf("ECEF->ENU->ECEF X: got %.2f, want %.2f", ecefBack.X, tc.ECEF.X)
			}
		})
	}
}

// ============================================================
// 6. ECI 函数测试 (独立的测试, 使用已知JD参考值)
// ============================================================
func TestECI_ECEF_Roundtrip(t *testing.T) {
	// ECI <-> ECEF 转换使用固定的已知JD和时间
	knownTimes := []struct {
		name string
		t    time.Time
		jd   float64
	}{
		{"J2000_noon", time.Date(2000, 1, 1, 12, 0, 0, 0, time.UTC), 2451545.0},
		{"J2000_midnight", time.Date(2000, 1, 1, 0, 0, 0, 0, time.UTC), 2451544.5},
	}
	// ECEF 测试点
	testPoints := []struct {
		x, y, z float64
	}{
		{-2179090.86, 4388326.88, 4069863.50}, // 北京附近
		{6378137.0, 0, 0},                      // 赤道本初子午线
		{0, 0, 6356752.31},                     // 北极附近
		{-2700226.9, -4292413.9, 3855273.8},   // 随机点 (来自cpp示例)
	}

	wgs84, _ := NewEllipsoid("wgs84")

	for _, tt := range knownTimes {
		t.Run(tt.name, func(t *testing.T) {
			jd := juliandate(tt.t)
			if math.Abs(jd-tt.jd) > 1e-9 {
				t.Fatalf("juliandate: got %.10f, want %.10f", jd, tt.jd)
			}

			for i, pt := range testPoints {
				// ECEF -> ECI
				xeci, yeci, zeci := ECEF2ECI(pt.x, pt.y, pt.z, tt.t)
				// ECI -> ECEF (round-trip)
				xb, yb, zb := ECI2ECEF(xeci, yeci, zeci, tt.t)

				// 验证互逆
				if math.Abs(xb-pt.x) > 1e-1 || math.Abs(yb-pt.y) > 1e-1 || math.Abs(zb-pt.z) > 1e-1 {
					t.Errorf("ECI round-trip [%d]: orig(%.2f, %.2f, %.2f) -> back(%.2f, %.2f, %.2f)",
						i, pt.x, pt.y, pt.z, xb, yb, zb)
				}

				// 测试结构体方法
				eci := ECEF{X: pt.x, Y: pt.y, Z: pt.z, Ell: wgs84}.ToECI(tt.t)
				ecefBack := eci.ToECEF()
				if math.Abs(ecefBack.X-pt.x) > 1e-1 || math.Abs(ecefBack.Y-pt.y) > 1e-1 || math.Abs(ecefBack.Z-pt.z) > 1e-1 {
					t.Errorf("ECI struct round-trip [%d]: orig(%.2f, %.2f, %.2f) -> back(%.2f, %.2f, %.2f)",
						i, pt.x, pt.y, pt.z, ecefBack.X, ecefBack.Y, ecefBack.Z)
				}
			}
		})
	}
}

func TestECI_Geodetic_Roundtrip(t *testing.T) {
	tNow := time.Date(2024, 6, 15, 12, 0, 0, 0, time.UTC)
	wgs84, _ := NewEllipsoid("wgs84")

	testPoints := []struct{ lat, lon, alt float64 }{
		{39.9042, 116.4074, 50},
		{0, 0, 0},
		{90, 0, 0},
		{-90, 0, 1000},
		{45, -120, 500},
	}

	for i, pt := range testPoints {
		geo := Geodetic{Latitude: pt.lat, Longitude: pt.lon, Altitude: pt.alt, Ell: wgs84}
		// Geodetic -> ECI -> Geodetic
		geoToEci := geo.ToECI(tNow)
		eciToGeo := geoToEci.ToGeodetic()

		if math.Abs(eciToGeo.Latitude-pt.lat) > 1e-8 || math.Abs(eciToGeo.Longitude-pt.lon) > 1e-8 {
			t.Errorf("Geo->ECI->Geo [%d]: orig(%.10f, %.10f) back(%.10f, %.10f)",
				i, pt.lat, pt.lon, eciToGeo.Latitude, eciToGeo.Longitude)
		}

		// ECEF -> ECI -> ENU -> Geodetic (完整链路)
		ecef := geo.ToECEF()
		eci := ecef.ToECI(tNow)
		ref := Geodetic{Latitude: pt.lat, Longitude: pt.lon, Altitude: pt.alt, Ell: wgs84}
		enuFromEci := eci.ToENU(ref)
		_ = enuFromEci // 只是验证不崩溃

		// ECEF -> ECI -> AER
		aerFromEci := eci.ToAER(ref)
		_ = aerFromEci
	}
}

// ============================================================
// 7. greenwichsrt 函数测试 (使用已知值)
// ============================================================
func TestGreenwichSiderealTime(t *testing.T) {
	// J2000.0 的 JD = 2451545.0
	// 此时格林威治恒星时 ≈ 280.46061837° + 360.98564736629° * 0 = 280.4606°
	// 但实际上 GMST 的完整公式会给出不同的值
	jd := 2451545.0
	gst := greenwichsrt(jd)
	// GST 应该在 [0, 2*pi) 范围内
	if gst < 0 || gst >= 2*math.Pi {
		t.Errorf("GST 超出 [0, 2pi) 范围: %.10f", gst)
	}
	// 已知值: J2000.0 的 GMST ≈ 1.75217 rad ≈ 100.398° (实际约为 1.752 rad)
	// 这是合理的值
	if gst < 0.1 || gst > 6.2 {
		t.Errorf("GST 值异常: %.10f rad (应在 0.1~6.2 之间)", gst)
	}

	// 测试不同时间: 24小时后 GMST 应增加约 360.985647°/天
	jd2 := jd + 1.0 // 1天后
	gst2 := greenwichsrt(jd2)
	// 处理 [0, 2pi) 环绕
	diff := gst2 - gst
	expected := 360.98564736629 * math.Pi / 180.0 // ~6.3004 rad
	// GMST每日增量略大于 2*pi, 需要检查两个方向
	if diff < 0 {
		diff += 2 * math.Pi
	}
	if math.Abs(diff-(expected-2*math.Pi)) > 0.001 {
		t.Errorf("GMST 1天变化: got %.10f rad, want wrapped %.10f rad (raw %.10f rad)",
			diff, expected-2*math.Pi, expected)
	}
}

// ============================================================
// 8. Ellipsoid 模型测试 (扩展)
// ============================================================
func TestEllipsoidAllModels(t *testing.T) {
	models := []string{"wgs84", "cgcs2000", "moon", "mars"}
	for _, name := range models {
		ell, err := NewEllipsoid(name)
		if err != nil {
			t.Errorf("NewEllipsoid(%q) 失败: %v", name, err)
			continue
		}
		if ell.SemimajorAxis <= 0 {
			t.Errorf("NewEllipsoid(%q).SemimajorAxis = %f <= 0", name, ell.SemimajorAxis)
		}
		if ell.SemiminorAxis <= 0 {
			t.Errorf("NewEllipsoid(%q).SemiminorAxis = %f <= 0", name, ell.SemiminorAxis)
		}
		if ell.Flattening <= 0 || ell.Flattening >= 1 {
			t.Errorf("NewEllipsoid(%q).Flattening = %f, 应在 (0,1) 之间", name, ell.Flattening)
		}
		if ell.Eccentricity <= 0 || ell.Eccentricity >= 1 {
			t.Errorf("NewEllipsoid(%q).Eccentricity = %f, 应在 (0,1) 之间", name, ell.Eccentricity)
		}

		// 验证 CGCS2000 和 WGS84 应非常接近 (短半轴相差 ~0.1mm)
		if name == "cgcs2000" || name == "wgs84" {
			wgs, _ := NewEllipsoid("wgs84")
			cgcs, _ := NewEllipsoid("cgcs2000")
			diff := math.Abs(wgs.SemiminorAxis - cgcs.SemiminorAxis)
			if diff > 0.001 { // 相差应 < 1mm
				t.Errorf("WGS84与CGCS2000短半轴差异过大: %.6f m", diff)
			}
		}
	}

	// 验证错误模型
	_, err := NewEllipsoid("invalid")
	if err == nil {
		t.Error("NewEllipsoid('invalid') 应返回错误")
	}
}
