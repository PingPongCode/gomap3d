package gomap3d

import "time"

// ENU 东北天坐标系 (m)
type ENU struct {
	East, North, Up float64
	Ell             *Ellipsoid
}

// 转站心坐标系
func (enu *ENU) ToAER() AER {
	azimuth, elevation, srange := ENU2AER(enu.East, enu.North, enu.Up)
	return AER{
		Azimuth:   azimuth,
		Elevation: elevation,
		SRange:    srange,
		Ell:       enu.Ell,
	}
}

//转大地坐标系
func (enu *ENU) ToGeodetic(ref Geodetic) Geodetic {
	x, y, z := ENU2ECEF(enu.East, enu.North, enu.Up, ref.Latitude, ref.Longitude, ref.Altitude, enu.Ell)
	latitude, longitude, altitude := ECEF2Geodetic(x, y, z, enu.Ell)
	return Geodetic{
		Latitude:  latitude,
		Longitude: longitude,
		Altitude:  altitude,
		Ell:       enu.Ell,
	}
}

// 转地心地固坐标系
func (enu *ENU) ToECEF(ref Geodetic) ECEF {
	x, y, z := ENU2ECEF(enu.East, enu.North, enu.Up, ref.Latitude, ref.Longitude, ref.Altitude, enu.Ell)
	return ECEF{
		X:   x,
		Y:   y,
		Z:   z,
		Ell: enu.Ell,
	}
}

// 转地心惯性坐标系
func (enu *ENU) ToECI(ref Geodetic, t time.Time) ECI {
	x, y, z := ENU2ECEF(enu.East, enu.North, enu.Up, ref.Latitude, ref.Longitude, ref.Altitude, enu.Ell)
	xECI, yECI, zECI := ECEF2ECI(x, y, z, t)
	return ECI{
		X:   xECI,
		Y:   yECI,
		Z:   zECI,
		T:   t,
		Ell: enu.Ell,
	}
}
