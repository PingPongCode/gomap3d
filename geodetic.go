package gomap3d

import "time"

// LLA 大地坐标系 (纬度, 经度, 高度) (°, °, m)
type Geodetic struct {
	Latitude  float64
	Longitude float64
	Altitude  float64
	Ell       *Ellipsoid
}

// 转地心地固坐标系
func (geo *Geodetic) ToECEF() ECEF {
	x, y, z := Geodetic2ECEF(geo.Latitude, geo.Longitude, geo.Altitude, geo.Ell)
	return ECEF{
		X:   x,
		Y:   y,
		Z:   z,
		Ell: geo.Ell,
	}
}

// 转东北天坐标系
func (geo *Geodetic) ToENU(ref Geodetic) ENU {
	x, y, z := Geodetic2ECEF(geo.Latitude, geo.Longitude, geo.Altitude, geo.Ell)
	east, north, up := ECEF2ENU(x, y, z, ref.Latitude, ref.Longitude, ref.Altitude, geo.Ell)
	return ENU{
		East:  east,
		North: north,
		Up:    up,
		Ell:   geo.Ell,
	}
}

// 转站心坐标系
func (geo *Geodetic) ToAER(ref Geodetic) AER {
	x, y, z := Geodetic2ECEF(geo.Latitude, geo.Longitude, geo.Altitude, geo.Ell)
	east, north, up := ECEF2ENU(x, y, z, ref.Latitude, ref.Longitude, ref.Altitude, geo.Ell)
	azimuth, elevation, srange := ENU2AER(east, north, up)
	return AER{
		Azimuth:   azimuth,
		Elevation: elevation,
		SRange:    srange,
		Ell:       geo.Ell,
	}
}

// 转地心惯性坐标系
func (geo *Geodetic) ToECI(t time.Time) ECI {
	x, y, z := Geodetic2ECEF(geo.Latitude, geo.Longitude, geo.Altitude, geo.Ell)
	xECI, yECI, zECI := ECEF2ECI(x, y, z, t)
	return ECI{
		X:   xECI,
		Y:   yECI,
		Z:   zECI,
		T:   t,
		Ell: geo.Ell,
	}
}
