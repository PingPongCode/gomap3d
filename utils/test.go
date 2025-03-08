package utils

import (
	"math"
	"testing"
	"time"

	"github.com/PingPongCode/gomap3d"
)

const epsilon = 1e-6 // 浮点比较精度

func TestAERToENU() {
	t := testing.T{}
	ell, _ := gomap3d.NewEllipsoid("wgs84")
	aer := gomap3d.AER{
		Azimuth:   45.0,
		Elevation: 30.0,
		SRange:    2000.0,
		Ell:       ell,
	}

	enu := aer.ToENU()

	expectedEast := 612.372436
	expectedNorth := 612.372436
	expectedUp := 1000.0

	if math.Abs(enu.East-expectedEast) > epsilon {
		t.Errorf("East分量错误: 期望 %.6f, 实际 %.6f", expectedEast, enu.East)
	}
	if math.Abs(enu.North-expectedNorth) > epsilon {
		t.Errorf("North分量错误: 期望 %.6f, 实际 %.6f", expectedNorth, enu.North)
	}
	if math.Abs(enu.Up-expectedUp) > epsilon {
		t.Errorf("Up分量错误: 期望 %.6f, 实际 %.6f", expectedUp, enu.Up)
	}
}

func TestGeodeticRoundTrip() {
	t := testing.T{}
	ell, _ := gomap3d.NewEllipsoid("cgcs2000")
	original := gomap3d.Geodetic{
		Latitude:  35.6895,
		Longitude: 139.6917,
		Altitude:  131.0,
		Ell:       ell,
	}

	// Geodetic -> ECEF -> Geodetic
	ecef := original.ToECEF()
	converted := ecef.ToGeodetic()

	if math.Abs(original.Latitude-converted.Latitude) > 1e-9 {
		t.Errorf("纬度不匹配: 原始 %.9f, 转换后 %.9f", original.Latitude, converted.Latitude)
	}
	if math.Abs(math.Mod(original.Longitude-converted.Longitude, 360)) > 1e-9 {
		t.Errorf("经度不匹配: 原始 %.9f, 转换后 %.9f", original.Longitude, converted.Longitude)
	}
	if math.Abs(original.Altitude-converted.Altitude) > 1e-3 {
		t.Errorf("高度不匹配: 原始 %.3f, 转换后 %.3f", original.Altitude, converted.Altitude)
	}
}

func TestECIECEFConversion() {
	t := testing.T{}
	ell, _ := gomap3d.NewEllipsoid("wgs84")
	tTime := time.Date(2023, 1, 1, 0, 0, 0, 0, time.UTC)

	// 初始ECEF坐标
	ecef := gomap3d.ECEF{
		X:   6378137.0,
		Y:   0.0,
		Z:   0.0,
		Ell: ell,
	}

	// ECEF -> ECI -> ECEF
	eci := ecef.ToECI(tTime)
	converted := eci.ToECEF()

	if math.Abs(ecef.X-converted.X) > 1e-3 {
		t.Errorf("X轴不匹配: 原始 %.3f, 转换后 %.3f", ecef.X, converted.X)
	}
	if math.Abs(ecef.Y-converted.Y) > 1e-3 {
		t.Errorf("Y轴不匹配: 原始 %.3f, 转换后 %.3f", ecef.Y, converted.Y)
	}
	if math.Abs(ecef.Z-converted.Z) > 1e-3 {
		t.Errorf("Z轴不匹配: 原始 %.3f, 转换后 %.3f", ecef.Z, converted.Z)
	}
}

func TestEllipsoidModels() {
	t := testing.T{}
	tests := []struct {
		model    string
		expected gomap3d.Ellipsoid
	}{
		{
			model: "wgs84",
			expected: gomap3d.Ellipsoid{
				SemimajorAxis: 6378137.0,
				SemiminorAxis: 6356752.31424518,
			},
		},
		{
			model: "mars",
			expected: gomap3d.Ellipsoid{
				SemimajorAxis: 3396190,
				SemiminorAxis: 3376097.80585952,
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.model, func(t *testing.T) {
			ell, err := gomap3d.NewEllipsoid(tt.model)
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
	_, err := gomap3d.NewEllipsoid("invalid")
	if err == nil {
		t.Error("预期错误未触发")
	}
}
