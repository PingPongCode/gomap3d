package gomap3d

import "time"

// ECEF 地心地固坐标系 (m)
type ECEF struct {
	X, Y, Z float64
	Ell     *Ellipsoid
}

// 转大地坐标系
func (ecef ECEF) ToGeodetic() Geodetic {
	lat, lon, alt := ECEF2Geodetic(ecef.X, ecef.Y, ecef.Z, ecef.Ell)
	return Geodetic{
		Latitude:  lat,
		Longitude: lon,
		Altitude:  alt,
		Ell:       ecef.Ell,
	}
}

// 转地心惯性坐标系
func (ecef ECEF) ToECI(t time.Time) ECI {
	x, y, z := ECEF2ECI(ecef.X, ecef.Y, ecef.Z, t)
	return ECI{
		X:   x,
		Y:   y,
		Z:   z,
		T:   t,
		Ell: ecef.Ell,
	}

}

// 转东北天坐标系
func (ecef ECEF) ToENU(ref Geodetic) ENU {
	east, north, up := ECEF2ENU(
		ecef.X, ecef.Y, ecef.Z,
		ref.Latitude, ref.Longitude, ref.Altitude,
		ecef.Ell)
	return ENU{
		East:  east,
		North: north,
		Up:    up,
		Ell:   ecef.Ell,
	}
}

// 转站心坐标系
func (ecef ECEF) ToAER(ref Geodetic) AER {
	east, north, up := ECEF2ENU(
		ecef.X, ecef.Y, ecef.Z,
		ref.Latitude, ref.Longitude, ref.Altitude,
		ecef.Ell)
	srange, azimuth, elevation := ENU2AER(east, north, up)
	return AER{
		Azimuth:   azimuth,
		Elevation: elevation,
		SRange:    srange,
		Ell:       ecef.Ell,
	}
}
