/* velocity_c.h -- 纯 C 实现的速度坐标转换函数
 *
 * 无需链接 Go DLL，直接编译到 C/C++ 程序中。
 * 公式与 gomap3d/velocity.go 中的 Go 实现一致。
 */

#ifndef VELOCITY_C_H
#define VELOCITY_C_H

#include <math.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* 常数定义（兼容 -std=c99 不提供 M_PI 的情况） */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===== 地球常数 ===== */
#define WE_C 7.2921150e-5  /* 地球自转角速度 (rad/s) */
#define RE_C 6378137.0     /* 地球赤道半径 (m) */
#define F_C  (1.0 / 298.257223563) /* 地球扁率 */
#define TAU_C 6.283185307179586 /* 2*PI */

/* ===== 儒略日与格林威治恒星时 ===== */

/* juliandate_c 计算 Unix 时间戳对应的儒略日
 * 基于 Meeus "Astronomical Algorithms" 公式，与 gomap3d 的 Go 实现一致 */
static inline double juliandate_c(double timestamp_sec) {
    int64_t sec = (int64_t)timestamp_sec;
    double frac = timestamp_sec - (double)sec;
    // Unix 时间戳 0 对应 1970-01-01 00:00:00 UTC
    // 转换为儒略日：
    // 1970-01-01 的儒略日为 2440587.5
    return 2440587.5 + (double)sec / 86400.0 + frac / 86400.0;
}

/* greenwichsrt_c 计算格林威治恒星时（弧度） */
static inline double greenwichsrt_c(double jd) {
    double t = (jd - 2451545.0) / 36525.0;
    double gmst = 67310.54841 +
        (876600LL * 3600 + 8640184.812866) * t +
        0.093104 * t * t -
        6.2e-6 * t * t * t;
    double gmstRad = fmod(gmst * TAU_C / 86400.0, TAU_C);
    if (gmstRad < 0) gmstRad += TAU_C;
    return gmstRad;
}

/* ===== 辅助函数 ===== */

/* 绕 Z 轴旋转 angle 弧度的 3x3 旋转矩阵 */
static inline void R3_c(double angle, double m[3][3]) {
    double c = cos(angle);
    double s = sin(angle);
    m[0][0] = c;  m[0][1] = s;  m[0][2] = 0;
    m[1][0] = -s; m[1][1] = c;  m[1][2] = 0;
    m[2][0] = 0;  m[2][1] = 0;  m[2][2] = 1;
}

/* 矩阵转置 */
static inline void transpose3_c(double m[3][3], double out[3][3]) {
    out[0][0] = m[0][0]; out[0][1] = m[1][0]; out[0][2] = m[2][0];
    out[1][0] = m[0][1]; out[1][1] = m[1][1]; out[1][2] = m[2][1];
    out[2][0] = m[0][2]; out[2][1] = m[1][2]; out[2][2] = m[2][2];
}

/* 3x3 矩阵乘以 3 维向量 */
static inline void mulMV_c(double m[3][3], double v[3], double out[3]) {
    out[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
    out[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
    out[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}

/* ===== 大地坐标与 ECEF 转换 ===== */

/* geodetic2ecef_c 将经纬高转换为 ECEF 坐标 */
static inline void geodetic2ecef_c(double latDeg, double lonDeg, double alt,
                                    double *x, double *y, double *z) {
    double latRad = latDeg * M_PI / 180.0;
    double lonRad = lonDeg * M_PI / 180.0;
    double e2 = 2.0 * F_C - F_C * F_C;
    double n = RE_C / sqrt(1.0 - e2 * sin(latRad) * sin(latRad));
    *x = (n + alt) * cos(latRad) * cos(lonRad);
    *y = (n + alt) * cos(latRad) * sin(lonRad);
    *z = (n * (1.0 - e2) + alt) * sin(latRad);
}

/* ecef2enu_c 将 ECEF 坐标转换为 ENU 坐标（需要参考点） */
static inline void ecef2enu_c(double x, double y, double z,
                               double lat0Deg, double lon0Deg, double h0,
                               double *e, double *n, double *u) {
    double x0, y0, z0;
    geodetic2ecef_c(lat0Deg, lon0Deg, h0, &x0, &y0, &z0);
    double dx = x - x0, dy = y - y0, dz = z - z0;
    double lat0Rad = lat0Deg * M_PI / 180.0;
    double lon0Rad = lon0Deg * M_PI / 180.0;
    *e = -sin(lon0Rad)*dx + cos(lon0Rad)*dy;
    *n = -sin(lat0Rad)*cos(lon0Rad)*dx - sin(lat0Rad)*sin(lon0Rad)*dy + cos(lat0Rad)*dz;
    *u = cos(lat0Rad)*cos(lon0Rad)*dx + cos(lat0Rad)*sin(lon0Rad)*dy + sin(lat0Rad)*dz;
}

/* enu2ecef_c 将 ENU 坐标转换为 ECEF 坐标（需要参考点） */
static inline void enu2ecef_c(double e, double n, double u,
                               double lat0Deg, double lon0Deg, double h0,
                               double *x, double *y, double *z) {
    double lat0Rad = lat0Deg * M_PI / 180.0;
    double lon0Rad = lon0Deg * M_PI / 180.0;
    double dx = -sin(lon0Rad)*e - sin(lat0Rad)*cos(lon0Rad)*n + cos(lat0Rad)*cos(lon0Rad)*u;
    double dy = cos(lon0Rad)*e - sin(lat0Rad)*sin(lon0Rad)*n + cos(lat0Rad)*sin(lon0Rad)*u;
    double dz = cos(lat0Rad)*n + sin(lat0Rad)*u;
    double x0, y0, z0;
    geodetic2ecef_c(lat0Deg, lon0Deg, h0, &x0, &y0, &z0);
    *x = x0 + dx; *y = y0 + dy; *z = z0 + dz;
}

/* enu2aer_c 将 ENU 坐标转换为 AER（方位角、俯仰角、斜距） */
static inline void enu2aer_c(double e, double n, double u,
                              double *azDeg, double *elDeg, double *srange) {
    double r = sqrt(e*e + n*n);
    *srange = sqrt(r*r + u*u);
    *elDeg = atan2(u, r) * 180.0 / M_PI;
    *azDeg = fmod(atan2(e, n) * 180.0 / M_PI + 360.0, 360.0);
}

/* ===== 速度转换函数 ===== */

/* ecefVel2eciVel_c: ECEF 速度 → ECI 速度
 * v_eci = R3(-GST) * (v_ecef + we × r_ecef) */
static inline void ecefVel2eciVel_c(double vx, double vy, double vz,
                                     double x, double y, double z,
                                     double timestamp,
                                     double *vxEci, double *vyEci, double *vzEci) {
    double jd = juliandate_c(timestamp);
    double gst = greenwichsrt_c(jd);
    double rot[3][3], tRot[3][3];
    R3_c(gst, rot);
    transpose3_c(rot, tRot);

    /* we × r_ecef */
    double wxr[3] = {-WE_C * y, WE_C * x, 0.0};
    double v[3] = {vx + wxr[0], vy + wxr[1], vz + wxr[2]};
    double result[3];
    mulMV_c(tRot, v, result);
    *vxEci = result[0]; *vyEci = result[1]; *vzEci = result[2];
}

/* eciVel2ecefVel_c: ECI 速度 → ECEF 速度（逆转换）
 * v_ecef = R3(GST) * v_eci - we × r_ecef
 * 其中 r_ecef = R3(GST) * r_eci */
static inline void eciVel2ecefVel_c(double vx, double vy, double vz,
                                     double x, double y, double z,
                                     double timestamp,
                                     double *vxEcef, double *vyEcef, double *vzEcef) {
    double jd = juliandate_c(timestamp);
    double gst = greenwichsrt_c(jd);
    double rot[3][3];
    R3_c(gst, rot);

    /* r_ecef = R3(GST) * r_eci */
    double rEci[3] = {x, y, z};
    double rEcef[3];
    mulMV_c(rot, rEci, rEcef);

    /* v_rot = R3(GST) * v_eci */
    double vEci[3] = {vx, vy, vz};
    double vRot[3];
    mulMV_c(rot, vEci, vRot);

    /* v_ecef = v_rot - we × r_ecef */
    *vxEcef = vRot[0] + WE_C * rEcef[1];
    *vyEcef = vRot[1] - WE_C * rEcef[0];
    *vzEcef = vRot[2];
}

/* enuVel2ecefVel_c: ENU 速度 → ECEF 速度 */
static inline void enuVel2ecefVel_c(double eVel, double nVel, double uVel,
                                     double latDeg, double lonDeg,
                                     double *vx, double *vy, double *vz) {
    double x0, y0, z0, x1, y1, z1;
    geodetic2ecef_c(latDeg, lonDeg, 0, &x0, &y0, &z0);
    enu2ecef_c(eVel, nVel, uVel, latDeg, lonDeg, 0, &x1, &y1, &z1);
    *vx = x1 - x0; *vy = y1 - y0; *vz = z1 - z0;
}

/* ecefVel2enuVel_c: ECEF 速度 → ENU 速度 */
static inline void ecefVel2enuVel_c(double vx, double vy, double vz,
                                     double latDeg, double lonDeg,
                                     double *e, double *n, double *u) {
    double x0, y0, z0;
    geodetic2ecef_c(latDeg, lonDeg, 0, &x0, &y0, &z0);
    ecef2enu_c(x0+vx, y0+vy, z0+vz, latDeg, lonDeg, 0, e, n, u);
}

/* aerDeriv2enuVel_c: AER 变化率（RAE 导数）→ ENU 速度 */
static inline void aerDeriv2enuVel_c(double R, double azDeg, double elDeg,
                                      double dR, double dAzDeg, double dElDeg,
                                      double *e, double *n, double *u) {
    double azRad = azDeg * M_PI / 180.0;
    double elRad = elDeg * M_PI / 180.0;
    double dAzRad = dAzDeg * M_PI / 180.0;
    double dElRad = dElDeg * M_PI / 180.0;
    double cosEl = cos(elRad), sinEl = sin(elRad);
    double cosAz = cos(azRad), sinAz = sin(azRad);

    *e = dR*cosEl*sinAz - R*sinEl*dElRad*sinAz + R*cosEl*cosAz*dAzRad;
    *n = dR*cosEl*cosAz - R*sinEl*dElRad*cosAz - R*cosEl*sinAz*dAzRad;
    *u = dR*sinEl + R*cosEl*dElRad;
}

/* enuVel2aerDeriv_c: ENU 速度 → AER 变化率（RAE 导数，逆转换） */
static inline void enuVel2aerDeriv_c(double eVel, double nVel, double uVel,
                                      double R, double azDeg, double elDeg,
                                      double *dR, double *dAzDeg, double *dElDeg) {
    double azRad = azDeg * M_PI / 180.0;
    double elRad = elDeg * M_PI / 180.0;
    double cosEl = cos(elRad), sinEl = sin(elRad);
    double cosAz = cos(azRad), sinAz = sin(azRad);

    /* dR = velocity projection onto line-of-sight */
    *dR = eVel*cosEl*sinAz + nVel*cosEl*cosAz + uVel*sinEl;

    if (fabs(cosEl) > 1e-12) {
        *dElDeg = ((uVel - *dR * sinEl) / (R * cosEl)) * 180.0 / M_PI;
        *dAzDeg = ((eVel*cosAz - nVel*sinAz) / (R * cosEl)) * 180.0 / M_PI;
    } else {
        *dElDeg = 0;
        *dAzDeg = 0;
    }
}

/* aerDeriv2eciVel_c: AER 变化率 → ECI 速度（一站式） */
static inline void aerDeriv2eciVel_c(double R, double azDeg, double elDeg,
                                      double dR, double dAzDeg, double dElDeg,
                                      double latDeg, double lonDeg,
                                      double timestamp,
                                      double *vx, double *vy, double *vz) {
    double eVel, nVel, uVel;
    aerDeriv2enuVel_c(R, azDeg, elDeg, dR, dAzDeg, dElDeg, &eVel, &nVel, &uVel);

    double vxE, vyE, vzE;
    enuVel2ecefVel_c(eVel, nVel, uVel, latDeg, lonDeg, &vxE, &vyE, &vzE);

    /* 计算 ECEF 位置 */
    double azRad = azDeg * M_PI / 180.0;
    double elRad = elDeg * M_PI / 180.0;
    double px = R * cos(elRad) * sin(azRad);
    double py = R * cos(elRad) * cos(azRad);
    double pz = R * sin(elRad);
    double ex, ey, ez;
    enu2ecef_c(px, py, pz, latDeg, lonDeg, 0, &ex, &ey, &ez);

    ecefVel2eciVel_c(vxE, vyE, vzE, ex, ey, ez, timestamp, vx, vy, vz);
}

#ifdef __cplusplus
}
#endif

#endif /* VELOCITY_C_H */
