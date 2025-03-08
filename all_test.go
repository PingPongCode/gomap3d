package gomap3d

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
	"testing"
	"time"
)

const epsilon = 1e-3 // 浮点比较精度
type ObservationData struct {
	Observation []float64
	Data        [][]float64
}

type TestData struct {
	DataColumns        []string
	ObservationColumns []string
	Data               []ObservationData
}

func TestCoordinateTransform(t *testing.T) {
	ell, _ := NewEllipsoid("wgs84")
	testData := &TestData{}
	var fileContent, err = os.ReadFile("data/randomPoints.json")
	if err != nil {
		t.Error(err)
	}
	json.Unmarshal(fileContent, testData)
	if err != nil {
		t.Error(err)
	}

	for i, observationData := range testData.Data {
		var observation = observationData.Observation
		var points = observationData.Data
		var ref = Geodetic{
			Latitude:  observation[0],
			Longitude: observation[1],
			Altitude:  observation[2],
		}
		for j, point := range points {
			var timestamp int64 = int64(point[0])
			var tNow = time.Unix(timestamp, 0)
			azimuth, elevation, srange := point[1], point[2], point[3]
			east, north, up := point[4], point[5], point[6]
			xECEF, yECEF, zECEF := point[7], point[8], point[9]
			latitude, longitude, altitude := point[10], point[11], point[12]
			xECI, yECI, zECI := point[13], point[14], point[15]

			// 正运算
			var aer = AER{point[1], point[2], point[3], ell}
			var enu = aer.ToENU()
			var ecef = aer.ToECEF(ref)
			var geodetic = aer.ToGeodetic(ref)
			var eci = aer.ToECI(ref, tNow)
			// 逆运算
			var ecef2 = eci.ToECEF()
			var geodetic2 = eci.ToGeodetic()
			var enu2 = eci.ToENU(ref)
			var aer2 = eci.ToAER(ref)
			var needOut = false
			// 测试正逆运算之间是否可逆
			if (math.Abs(aer2.Azimuth-aer.Azimuth) > epsilon) || (math.Abs(aer2.Elevation-aer.Elevation) > epsilon) || (math.Abs(aer2.SRange-aer.SRange) > epsilon) {
				needOut = true
				t.Errorf("正逆运算之间不可逆")
			}
			// 测试各项值和真值之间误差是否小于预设的容差值epsilon
			// AER
			if (math.Abs(aer.Azimuth-azimuth) > epsilon) || (math.Abs(aer.Elevation-elevation) > epsilon) || (math.Abs(aer.SRange-srange) > epsilon) {
				needOut = true
				t.Errorf("AER 误差超出")
			}
			// ENU
			if (math.Abs(enu.East-east) > epsilon) || (math.Abs(enu.North-north) > epsilon) || (math.Abs(enu.Up-up) > epsilon) {
				needOut = true
				t.Errorf("ENU 误差超出")
			}
			// ECEF
			if (math.Abs(ecef.X-xECEF) > epsilon) || (math.Abs(ecef.Y-yECEF) > epsilon) || (math.Abs(ecef.Z-zECEF) > epsilon) {
				needOut = true
				t.Errorf("ECEF 误差超出")
			}
			// Geodetic
			if (math.Abs(geodetic.Latitude-latitude) > epsilon) || (math.Abs(geodetic.Longitude-longitude) > epsilon) || (math.Abs(geodetic.Altitude-altitude) > epsilon) {
				needOut = true
				t.Errorf("Geodetic 误差超出")
			}
			//ECI
			if (math.Abs(eci.X-xECI) > epsilon) || (math.Abs(eci.Y-yECI) > epsilon) || (math.Abs(eci.Z-zECI) > epsilon) {
				needOut = true
				t.Errorf("ECI 误差超出")
			}

			// AER逆
			if (math.Abs(aer2.Azimuth-azimuth) > epsilon) || (math.Abs(aer2.Elevation-elevation) > epsilon) || (math.Abs(aer2.SRange-srange) > epsilon) {
				needOut = true
				t.Errorf("AER逆 误差超出")
			}
			// ENU逆
			if (math.Abs(enu2.East-east) > epsilon) || (math.Abs(enu2.North-north) > epsilon) || (math.Abs(enu2.Up-up) > epsilon) {
				needOut = true
				t.Errorf("ENU逆 误差超出")
			}
			// ECEF逆
			if (math.Abs(ecef2.X-xECEF) > epsilon) || (math.Abs(ecef2.Y-yECEF) > epsilon) || (math.Abs(ecef2.Z-zECEF) > epsilon) {
				needOut = true
				t.Errorf("ECEF逆 误差超出")
			}
			// Geodetic逆
			if (math.Abs(geodetic2.Latitude-latitude) > epsilon) || (math.Abs(geodetic2.Longitude-longitude) > epsilon) || (math.Abs(geodetic2.Altitude-altitude) > epsilon) {
				needOut = true
				t.Errorf("Geodetic逆 误差超出")
			}
			if needOut {
				fmt.Print(
					"order ", i, "-", j, " timestamp ", timestamp, "\n",
					"azimuth: ", azimuth, aer.Azimuth, aer2.Azimuth, " elevation: ", elevation, aer.Elevation, aer2.Elevation, " srange: ", srange, aer.SRange, aer2.SRange, "\n",
					"east: ", east, enu.East, enu2.East, " north: ", north, enu.North, enu2.North, " up: ", up, enu.Up, enu2.Up, "\n",
					"xECEF: ", xECEF, ecef.X, ecef2.X, " yECEF: ", yECEF, ecef.Y, ecef2.Y, " zECEF: ", zECEF, ecef.Z, ecef.Z, "\n",
					"latitude: ", latitude, geodetic.Latitude, geodetic2.Latitude, " longitude: ", longitude, geodetic.Longitude, geodetic2.Longitude, " altitude: ", altitude, geodetic.Altitude, geodetic2.Altitude, "\n",
					"xECI: ", xECI, eci.X, " yECI: ", yECI, eci.Y, " zECI: ", zECI, eci.Z, "\n",
				)
				fmt.Print("---------------------------------------------", "\n")
			}
		}
	}

}

func TestEllipsoidModels(t *testing.T) {
	tests := []struct {
		model    string
		expected Ellipsoid
	}{
		{
			model: "wgs84",
			expected: Ellipsoid{
				SemimajorAxis: 6378137.0,
				SemiminorAxis: 6356752.31424518,
			},
		},
		{
			model: "mars",
			expected: Ellipsoid{
				SemimajorAxis: 3396190,
				SemiminorAxis: 3376097.80585952,
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.model, func(t *testing.T) {
			ell, err := NewEllipsoid(tt.model)
			if err != nil {
				t.Fatalf("创建椭球体失败: %v", err)
			}

			if math.Abs(ell.SemimajorAxis-tt.expected.SemimajorAxis) > 1e-6 {
				t.Errorf("长半轴不匹配: 期望 %.6f, 实际 %.6f",
					tt.expected.SemimajorAxis, ell.SemimajorAxis)
			}
			if math.Abs(ell.SemiminorAxis-tt.expected.SemiminorAxis) > 1e-6 {
				t.Errorf("短半轴不匹配: 期望 %.6f, 实际 %.6f",
					tt.expected.SemiminorAxis, ell.SemiminorAxis)
			}
		})
	}

	// 测试错误处理
	_, err := NewEllipsoid("invalid")
	if err == nil {
		t.Error("预期错误未触发")
	}
}
