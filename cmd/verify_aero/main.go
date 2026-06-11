package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
	"time"

	"github.com/PingPongCode/gomap3d"
)

type MaxCheckInput struct {
	ID                string  `json:"id"`
	Name              string  `json:"name"`
	RadarLatDeg       float64 `json:"radar_lat_deg"`
	RadarLonDeg       float64 `json:"radar_lon_deg"`
	RadarHM           float64 `json:"radar_h_m"`
	RM                float64 `json:"R_m"`
	AzDeg             float64 `json:"Az_deg"`
	ElDeg             float64 `json:"El_deg"`
	VxEci             float64 `json:"Vx_eci"`
	VyEci             float64 `json:"Vy_eci"`
	VzEci             float64 `json:"Vz_eci"`
	TDatenum          float64 `json:"T_datenum"`
	ExpectedIsSpace   bool    `json:"expected_isSpace"`
	ExpectedScenario  string  `json:"expected_scenario"`
}

type MaxCheckRef struct {
	ID           string  `json:"id"`
	Name         string  `json:"name"`
	IsSpaceTarget bool   `json:"isSpaceTarget"`
	DtLaunchS    *float64 `json:"dt_launch_s"`
	DtImpactS    *float64 `json:"dt_impact_s"`
	LaunchLat    *float64 `json:"launch_lat_deg"`
	LaunchLon    *float64 `json:"launch_lon_deg"`
	ImpactLat    *float64 `json:"impact_lat_deg"`
	ImpactLon    *float64 `json:"impact_lon_deg"`
}

func datenumToTime(dn float64) time.Time {
	jd := dn + 1721058.5
	dayFrac := jd - math.Floor(jd)
	jdInt := int64(math.Floor(jd))

	L := jdInt + 68569
	N := 4 * L / 146097
	L = L - (146097*N+3)/4
	I := 4000 * (L + 1) / 1461001
	L = L - 1461*I/4 + 31
	J := 80 * L / 2447
	day := int(L - 2447*J/80)
	L = J / 11
	month := int(J + 2 - 12*L)
	year := int(100*(N-49) + I + L)

	secs := dayFrac * 86400.0
	h := int(secs) / 3600
	m := (int(secs) % 3600) / 60
	s := secs - float64(h*3600+m*60)

	return time.Date(year, time.Month(month), day, h, m, int(s),
		int((s-float64(int(s)))*1e9), time.UTC)
}

func main() {
	// 读取测试输入
	data, _ := os.ReadFile("check_input_maxcheck/test_mc_04_shortrange.json")
	var in MaxCheckInput
	json.Unmarshal(data, &in)

	tUTC := datenumToTime(in.TDatenum)
	ell, _ := gomap3d.NewEllipsoid("wgs84")

	// 测站 ECEF
	sx, sy, sz := gomap3d.Geodetic2ECEF(in.RadarLatDeg, in.RadarLonDeg, in.RadarHM, ell)

	// ENU → ECEF
	azR, elR := in.AzDeg*math.Pi/180, in.ElDeg*math.Pi/180
	e := in.RM * math.Cos(elR) * math.Sin(azR)
	n := in.RM * math.Cos(elR) * math.Cos(azR)
	u := in.RM * math.Sin(elR)
	tx, ty, tz := gomap3d.ENU2ECEF(e, n, u, in.RadarLatDeg, in.RadarLonDeg, in.RadarHM, ell)

	tgx, tgy, tgz := sx+tx, sy+ty, sz+tz

	// ECEF → ECI (新完整链)
	xi, yi, zi := gomap3d.ECEF2ECI(tgx, tgy, tgz, tUTC)
	// ECEF velocity → ECI velocity
	vxi, vyi, vzi := gomap3d.ECEFVel2ECIVel(in.VxEci, in.VyEci, in.VzEci, tgx, tgy, tgz, tUTC)

	fmt.Printf("Go Aero ECI position: (%.1f, %.1f, %.1f) km\n", xi/1000, yi/1000, zi/1000)
	fmt.Printf("Go Aero ECI velocity: (%.1f, %.1f, %.1f) m/s\n", vxi, vyi, vzi)

	// 读取 MATLAB 参考
	refData, _ := os.ReadFile("check_output_maxcheck/test_mc_04_shortrange_out.json")
	var ref MaxCheckRef
	json.Unmarshal(refData, &ref)

	fmt.Printf("\nMATLAB maxcheck:\n")
	fmt.Printf("  isSpaceTarget: %v\n", ref.IsSpaceTarget)
	if ref.DtLaunchS != nil {
		fmt.Printf("  dt_launch: %.1f s\n", *ref.DtLaunchS)
	}
	if ref.DtImpactS != nil {
		fmt.Printf("  dt_impact: %.1f s\n", *ref.DtImpactS)
	}
	if ref.LaunchLat != nil {
		fmt.Printf("  launch: (%.4f, %.4f)\n", *ref.LaunchLat, *ref.LaunchLon)
	}
	if ref.ImpactLat != nil {
		fmt.Printf("  impact: (%.4f, %.4f)\n", *ref.ImpactLat, *ref.ImpactLon)
	}
}
