package gomap3d

import (
	"fmt"
	"math"
)

// Ellipsoid 表示地球椭球体参数
type Ellipsoid struct {
	Model           string
	Name            string
	SemimajorAxis   float64
	SemiminorAxis   float64
	Flattening      float64
	ThirdFlattening float64
	Eccentricity    float64
}

var models = map[string]struct {
	Name string
	A    float64
	B    float64
}{
	// 地球椭球体
	// CGCS2000 坐标系
	"cgcs2000": {"CGCS-2000 (2008) ", 6378137.0, 6356752.31414},
	// WGS84 坐标系
	"wgs84": {"WGS-84 (1984)", 6378137.0, 6356752.31424518},
	// 月球
	"moon": {"Moon", 1738100, 1736000.0},
	// 火星
	"mars": {"Mars", 3396190, 3376097.80585952},
}

// NewEllipsoid 通过名称创建椭球体
func NewEllipsoid(model string) (*Ellipsoid, error) {
	m, ok := models[model]
	if !ok {
		return nil, fmt.Errorf("unknown ellipsoid model: %s", model)
	}

	a := m.A
	b := m.B
	f := (a - b) / a
	thirdF := (a - b) / (a + b)
	e := math.Sqrt(2*f - f*f)

	return &Ellipsoid{
		Model:           model,
		Name:            m.Name,
		SemimajorAxis:   a,
		SemiminorAxis:   b,
		Flattening:      f,
		ThirdFlattening: thirdF,
		Eccentricity:    e,
	}, nil
}
