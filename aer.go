package gomap3d

import "time"

// AER 站心坐标系 (距离,方位,俯仰) (m, °, °)
type AER struct {
	Azimuth   float64
	Elevation float64
	SRange    float64
	Ell       *Ellipsoid
}

// 转东北天坐标系
func (aer *AER) ToENU() ENU {
	east, north, up := AER2ENU(aer.Azimuth, aer.Elevation, aer.SRange)
	return ENU{
		East:  east,
		North: north,
		Up:    up,
		Ell:   aer.Ell,
	}
}

// 转地心地固坐标系
func (aer *AER) ToECEF(ref Geodetic) ECEF {
	east, north, up := AER2ENU(aer.Azimuth, aer.Elevation, aer.SRange)
	x, y, z := ENU2ECEF(east, north, up, ref.Latitude, ref.Longitude, ref.Altitude, aer.Ell)
	return ECEF{
		X:   x,
		Y:   y,
		Z:   z,
		Ell: aer.Ell,
	}
}

// 转大地坐标系
func (aer *AER) ToGeodetic(ref Geodetic) Geodetic {
	east, north, up := AER2ENU(aer.Azimuth, aer.Elevation, aer.SRange)
	x, y, z := ENU2ECEF(east, north, up, ref.Latitude, ref.Longitude, ref.Altitude, aer.Ell)
	latitude, longitude, altitude := ECEF2Geodetic(x, y, z, aer.Ell)
	return Geodetic{
		Latitude:  latitude,
		Longitude: longitude,
		Altitude:  altitude,
		Ell:       aer.Ell,
	}
}

// 转地心惯性坐标系
func (aer *AER) ToECI(ref Geodetic, t time.Time) ECI {
	east, north, up := AER2ENU(aer.Azimuth, aer.Elevation, aer.SRange)
	x, y, z := ENU2ECEF(east, north, up, ref.Latitude, ref.Longitude, ref.Altitude, aer.Ell)
	xECI, yECI, zECI := ECEF2ECI(x, y, z, t)
	return ECI{
		X:   xECI,
		Y:   yECI,
		Z:   zECI,
		Ell: aer.Ell,
	}

}
