/* test_velocity_c.c -- 纯 C 速度转换测试（无需链接任何 DLL） */
#include <stdio.h>
#include <math.h>
#include "velocity_c.h"

static int failures = 0;
static const double EPS = 1e-3;

static void check(const char* name, double got, double expected) {
    if (fabs(got - expected) > EPS) {
        printf("  FAIL %s: got %.6f, expected %.6f\n", name, got, expected);
        failures++;
    }
}

int main() {
    double timestamp = 1780472536.0; // 2026-06-03 07:42:16 UTC

    printf("=== 纯 C 速度转换测试 ===\n");

    /* 1. ECEF->ECI->ECEF 往返 */
    printf("\n1. ECEFVel->ECIVel->ECEFVel roundtrip:\n");
    double vx1, vy1, vz1, vx2, vy2, vz2;
    double rx_ecef = 5000000, ry_ecef = 3000000, rz_ecef = 2000000;
    ecefVel2eciVel_c(1500, 2500, 800, rx_ecef, ry_ecef, rz_ecef, timestamp, &vx1, &vy1, &vz1);
    /* 逆向需要 ECI 位置：由 ECEF2ECI 转换得到 */
    double jd = juliandate_c(timestamp);
    double gst = greenwichsrt_c(jd);
    double rot[3][3];
    R3_c(gst, rot);
    double rEci_in[3] = {rx_ecef, ry_ecef, rz_ecef};
    double rEci_out[3];
    double tRot[3][3];
    transpose3_c(rot, tRot);
    mulMV_c(tRot, rEci_in, rEci_out);  // r_eci = R3(-GST) * r_ecef
    eciVel2ecefVel_c(vx1, vy1, vz1, rEci_out[0], rEci_out[1], rEci_out[2], timestamp, &vx2, &vy2, &vz2);
    check("ECEF->ECI->ECEF Vx", vx2, 1500.0);
    check("ECEF->ECI->ECEF Vy", vy2, 2500.0);
    check("ECEF->ECI->ECEF Vz", vz2, 800.0);
    if (failures == 0) printf("  PASS\n");

    /* 2. ENU->ECEF->ENU 往返 */
    printf("\n2. ENUVel->ECEFVel->ENUVel roundtrip:\n");
    double vx3, vy3, vz3, e2, n2, u2;
    enuVel2ecefVel_c(500, -300, 200, 40.0, 125.0, &vx3, &vy3, &vz3);
    ecefVel2enuVel_c(vx3, vy3, vz3, 40.0, 125.0, &e2, &n2, &u2);
    check("ENU->ECEF->ENU E", e2, 500.0);
    check("ENU->ECEF->ENU N", n2, -300.0);
    check("ENU->ECEF->ENU U", u2, 200.0);
    if (failures == 0) printf("  PASS\n");

    /* 3. ECI->ECEF->ECI 往返 */
    printf("\n3. ECIVel->ECEFVel->ECIVel roundtrip:\n");
    double vx4, vy4, vz4, vx5, vy5, vz5;
    double rx_eci = 6378137, ry_eci = 0, rz_eci = 0;
    eciVel2ecefVel_c(3000, -2000, 1000, rx_eci, ry_eci, rz_eci, timestamp, &vx4, &vy4, &vz4);
    /* 逆向需要 ECEF 位置：由 ECI2ECEF 转换得到 */
    double rEcef_in[3] = {rx_eci, ry_eci, rz_eci};
    double rEcef_out[3];
    R3_c(gst, rot);  // rot = R3(GST)
    mulMV_c(rot, rEcef_in, rEcef_out);  // r_ecef = R3(GST) * r_eci
    ecefVel2eciVel_c(vx4, vy4, vz4, rEcef_out[0], rEcef_out[1], rEcef_out[2], timestamp, &vx5, &vy5, &vz5);
    check("ECI->ECEF->ECI Vx", vx5, 3000.0);
    check("ECI->ECEF->ECI Vy", vy5, -2000.0);
    check("ECI->ECEF->ECI Vz", vz5, 1000.0);
    if (failures == 0) printf("  PASS\n");

    /* 4. AER导数->ENU->AER导数 往返 */
    printf("\n4. AERDeriv->ENUVel->AERDeriv roundtrip:\n");
    double e3, n3, u3, dR2, dAz2, dEl2;
    aerDeriv2enuVel_c(200000, 45, 30, 1000, 0, 0, &e3, &n3, &u3);
    enuVel2aerDeriv_c(e3, n3, u3, 200000, 45, 30, &dR2, &dAz2, &dEl2);
    check("AER->ENU->AER dR", dR2, 1000.0);
    check("AER->ENU->AER dAz", dAz2, 0.0);
    check("AER->ENU->AER dEl", dEl2, 0.0);
    if (failures == 0) printf("  PASS\n");

    /* 5. 一站式 AER->ECI */
    printf("\n5. AERDeriv->ECIVel (one-shot):\n");
    double vx6, vy6, vz6;
    aerDeriv2eciVel_c(1200000, 180, 45, 3000, 0, 0.01, 40.0, 116.4, timestamp, &vx6, &vy6, &vz6);
    double speed = sqrt(vx6*vx6 + vy6*vy6 + vz6*vz6);
    printf("  ECI speed: %.2f m/s\n", speed);
    if (speed > 0 && speed < 20000) {
        printf("  PASS\n");
    } else {
        printf("  FAIL: unreasonable speed\n");
        failures++;
    }

    printf("\n=== ");
    if (failures == 0) {
        printf("ALL PASSED!");
    } else {
        printf("%d TESTS FAILED!", failures);
    }
    printf(" ===\n");

    return failures > 0 ? 1 : 0;
}
