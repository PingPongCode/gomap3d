package gomap3d

import (
	"math"
	"testing"
	"time"
)

// 儒略日已知参考值
// 参考: Meeus "Astronomical Algorithms", USNO, 在线转换器
var knownJDs = []struct {
	name string
	time time.Time
	jd   float64
	tol  float64 // 容差(天)
}{
	{
		name: "J2000.0 历元 (2000-01-01 12:00:00)",
		time: time.Date(2000, 1, 1, 12, 0, 0, 0, time.UTC),
		jd:   2451545.0,
		tol:  1e-9,
	},
	{
		name: "Unix 历元 (1970-01-01 00:00:00)",
		time: time.Date(1970, 1, 1, 0, 0, 0, 0, time.UTC),
		jd:   2440587.5,
		tol:  1e-9,
	},
	{
		name: "2019-01-01 00:00:00",
		time: time.Date(2019, 1, 1, 0, 0, 0, 0, time.UTC),
		jd:   2458484.5,
		tol:  1e-9,
	},
	{
		name: "2024-01-01 00:00:00",
		time: time.Date(2024, 1, 1, 0, 0, 0, 0, time.UTC),
		jd:   2460310.5,
		tol:  1e-9,
	},
	{
		name: "2000-03-01 00:00:00 (月份边界, 不调整月份)",
		time: time.Date(2000, 3, 1, 0, 0, 0, 0, time.UTC),
		jd:   2451604.5,
		tol:  1e-9,
	},
	{
		name: "2020-02-29 00:00:00 (闰日)",
		time: time.Date(2020, 2, 29, 0, 0, 0, 0, time.UTC),
		jd:   2458908.5,
		tol:  1e-9,
	},
	{
		name: "1850-01-01 00:00:00 (19世纪)",
		time: time.Date(1850, 1, 1, 0, 0, 0, 0, time.UTC),
		jd:   2396758.5,
		tol:  1e-9,
	},
	{
		name: "2000-01-01 06:30:15.500 (含时分秒毫秒)",
		time: time.Date(2000, 1, 1, 6, 30, 15, 500000000, time.UTC),
		jd:   2451544.5 + (6.0+30.0/60.0+15.5/3600.0)/24.0,
		tol:  1e-9,
	},
	{
		name: "2001-01-01 00:00:00 (新世纪午夜)",
		time: time.Date(2001, 1, 1, 0, 0, 0, 0, time.UTC),
		jd:   2451910.5,
		tol:  1e-9,
	},
	{
		name: "2001-01-01 12:00:00 (新世纪正午)",
		time: time.Date(2001, 1, 1, 12, 0, 0, 0, time.UTC),
		jd:   2451911.0,
		tol:  1e-9,
	},
	{
		name: "1999-12-31 00:00:00 (千禧年前夕)",
		time: time.Date(1999, 12, 31, 0, 0, 0, 0, time.UTC),
		jd:   2451543.5,
		tol:  1e-9,
	},
	{
		name: "1992-04-10 00:00:00 (Meeus 例3.a)",
		time: time.Date(1992, 4, 10, 0, 0, 0, 0, time.UTC),
		jd:   2448722.5,
		tol:  1e-9,
	},
	{
		name: "1957-10-04 00:00:00 (Sputnik 1发射日)",
		time: time.Date(1957, 10, 4, 0, 0, 0, 0, time.UTC),
		jd:   2436115.5,
		tol:  1e-9,
	},
	{
		name: "2026-06-03 10:47:41 (当前测试时间)",
		time: time.Date(2026, 6, 3, 10, 47, 41, 0, time.UTC),
		jd:   2461194.5 + (10.0+47.0/60.0+41.0/3600.0)/24.0,
		tol:  1e-9,
	},
}

func TestJulianDateKnownValues(t *testing.T) {
	for _, tc := range knownJDs {
		t.Run(tc.name, func(t *testing.T) {
			got := juliandate(tc.time)
			diff := math.Abs(got - tc.jd)
			if diff > tc.tol {
				t.Errorf("juliandate(%v) = %.15f\n  期望 %.15f\n  差异 %.2e (容差 %.0e)",
					tc.time, got, tc.jd, diff, tc.tol)
			} else {
				t.Logf("✓ got %.15f, want %.15f, diff %.2e ✓", got, tc.jd, diff)
			}
		})
	}
}

// 验证正午午夜关系
func TestJulianDateNoonMidnight(t *testing.T) {
	tests := []struct {
		date time.Time
		jd   float64
		msg  string
	}{
		{time.Date(2000, 1, 1, 12, 0, 0, 0, time.UTC), 2451545.0, "J2000正午"},
		{time.Date(2000, 1, 1, 0, 0, 0, 0, time.UTC), 2451544.5, "J2000午夜"},
		{time.Date(2024, 1, 1, 12, 0, 0, 0, time.UTC), 2460311.0, "2024正午"},
		{time.Date(2024, 1, 1, 0, 0, 0, 0, time.UTC), 2460310.5, "2024午夜"},
	}
	for _, tc := range tests {
		got := juliandate(tc.date)
		diff := math.Abs(got - tc.jd)
		if diff > 1e-9 {
			t.Errorf("%s: got %.10f, want %.10f, diff %.2e", tc.msg, got, tc.jd, diff)
		} else {
			t.Logf("✓ %s = %.10f ✓", tc.msg, got)
		}
	}
}

// 测试时间单调性
func TestJulianDateMonotonicity(t *testing.T) {
	base := time.Date(2024, 1, 1, 0, 0, 0, 0, time.UTC)
	prev := juliandate(base)
	for i := 1; i <= 3600; i++ {
		next := base.Add(time.Duration(i) * time.Second)
		jd := juliandate(next)
		if jd <= prev {
			t.Errorf("JD 不单调: t=%v, prev=%.15f, curr=%.15f", next, prev, jd)
			return
		}
		prev = jd
	}
	t.Log("✓ 时间单调性测试通过 (3600个连续秒点) ✓")
}

// 测试秒级精度
func TestJulianDatePrecision(t *testing.T) {
	t1 := time.Date(2024, 6, 1, 12, 0, 0, 0, time.UTC)
	t2 := t1.Add(time.Second) // 1秒后

	jd1 := juliandate(t1)
	jd2 := juliandate(t2)
	diff := jd2 - jd1
	expected := 1.0 / 86400.0 // 1秒 = 1/86400 天
	if math.Abs(diff-expected) > 1e-10 {
		t.Errorf("1秒差异应 = %.15e, 实际 = %.15e", expected, diff)
	} else {
		t.Logf("✓ 1秒精度: %.15e (期望 %.15e) ✓", diff, expected)
	}

	// 测试毫秒级精度
	t3 := t1.Add(100 * time.Millisecond)
	jd3 := juliandate(t3)
	diffMs := jd3 - jd1
	expectedMs := 0.1 / 86400.0 // 100ms = 0.1/86400 天
	if math.Abs(diffMs-expectedMs) > 1e-9 {
		t.Errorf("100ms差异应 = %.15e, 实际 = %.15e", expectedMs, diffMs)
	} else {
		t.Logf("✓ 100ms精度: %.15e ✓", diffMs)
	}
}

// 测试跨世纪日期 (验证格里高利历规则)
func TestJulianDateCenturies(t *testing.T) {
	// 2100年2月28日到3月1日 (2100不是闰年, 2月只有28天)
	feb28_2100 := time.Date(2100, 2, 28, 0, 0, 0, 0, time.UTC)
	mar1_2100 := time.Date(2100, 3, 1, 0, 0, 0, 0, time.UTC)
	diff := juliandate(mar1_2100) - juliandate(feb28_2100)
	if math.Abs(diff-1.0) > 1e-9 {
		t.Errorf("2100年非闰年: 2月28→3月1日应为1天, 实际 = %.10f", diff)
	} else {
		t.Logf("✓ 2100年非闰年验证通过: 2月28→3月1日 = %.10f天 ✓", diff)
	}

	// 2000年2月28日到3月1日 (2000是闰年, 有2月29日)
	feb28_2000 := time.Date(2000, 2, 28, 0, 0, 0, 0, time.UTC)
	mar1_2000 := time.Date(2000, 3, 1, 0, 0, 0, 0, time.UTC)
	diff2 := juliandate(mar1_2000) - juliandate(feb28_2000)
	if math.Abs(diff2-2.0) > 1e-9 {
		t.Errorf("2000年闰年: 2月28→3月1日应为2天, 实际 = %.10f", diff2)
	} else {
		t.Logf("✓ 2000年闰年验证通过: 2月28→3月1日 = %.10f天 ✓", diff2)
	}
}
