package gomap3d

import (
	"time"
)

// ECI 地心惯性坐标系 (m)
type ECI struct {
	X, Y, Z float64
	T       time.Time
	Ell     *Ellipsoid
}

// 转地心地固坐标系
func (eci *ECI) ToECEF() ECEF {
	x, y, z := ECI2ECEF(eci.X, eci.Y, eci.Z, eci.T)
	return ECEF{
		X:   x,
		Y:   y,
		Z:   z,
		Ell: eci.Ell,
	}
}

// 转大地坐标系
func (eci *ECI) ToGeodetic() Geodetic {
	x, y, z := ECI2ECEF(eci.X, eci.Y, eci.Z, eci.T)
	latitude, longitude, altitude := ECEF2Geodetic(x, y, z, eci.Ell)
	return Geodetic{
		Latitude:  latitude,
		Longitude: longitude,
		Altitude:  altitude,
		Ell:       eci.Ell,
	}
}

// 转东北天坐标系
func (eci *ECI) ToENU(ref Geodetic) ENU {
	x, y, z := ECI2ECEF(eci.X, eci.Y, eci.Z, eci.T)
	e, n, u := ECEF2ENU(
		x, y, z,
		ref.Latitude, ref.Longitude, ref.Altitude,
		eci.Ell)
	return ENU{
		East:  e,
		North: n,
		Up:    u,
		Ell:   eci.Ell,
	}
}

// 转站心坐标系
func (eci *ECI) ToAER(ref Geodetic) AER {
	x, y, z := ECI2ECEF(eci.X, eci.Y, eci.Z, eci.T)
	east, north, up := ECEF2ENU(
		x, y, z,
		ref.Latitude, ref.Longitude, ref.Altitude,
		eci.Ell)
	a, e, r := ENU2AER(east, north, up)
	return AER{
		Azimuth:   a,
		Elevation: e,
		SRange:    r,
		Ell:       eci.Ell,
	}
}
