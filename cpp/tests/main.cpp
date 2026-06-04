#include <iostream>
#include <cmath>
#include "gomap3d.h"

static int failures = 0;
static const double EPS = 1e-3;

static void check(const char* name, double got, double expected) {
    if (std::fabs(got - expected) > EPS) {
        std::cout << "  \u274c " << name << ": got " << got << ", expected " << expected << "\n";
        failures++;
    }
}

int main() {
    // 测试ENU2AER
    ENU2AER_return aer = ENU2AER(10.0, 20.0, 30.0);
    std::cout << "ENU2AER Result: "
              << "Azimuth: " << aer.r0 << " "
              << "Elevation: " << aer.r1 << " "
              << "Range: " << aer.r2 << "m\n";

    // 测试ECI/ECEF转换
    double timestamp =  1717000000;
    ECI2ECEF_return ecef = ECI2ECEF(6500000, 100000, 300000, timestamp);
    std::cout << "\nECI2ECEF Result: "
              << "X: " << ecef.r0 << "m "
              << "Y: " << ecef.r1 << "m "
              << "Z: " << ecef.r2 << "m\n";

    // ===== 速度转换测试 =====
    std::cout << "\n=== 速度转换测试 ===\n";

    // 1. ECEF->ECI->ECEF 往返
    std::cout << "\n1. ECEFVel\u2194ECIVel roundtrip:\n";
    double rx_ECEF = 5000000, ry_ECEF = 3000000, rz_ECEF = 2000000;
    ECEFVel2ECIVel_return v1 = ECEFVel2ECIVel(1500, 2500, 800, rx_ECEF, ry_ECEF, rz_ECEF, timestamp);
    // 逆向需要 ECI 位置: ECEF2ECI
    ECEF2ECI_return pos_eci = ECEF2ECI(rx_ECEF, ry_ECEF, rz_ECEF, timestamp);
    ECIVel2ECEFVel_return v2 = ECIVel2ECEFVel(v1.r0, v1.r1, v1.r2, pos_eci.r0, pos_eci.r1, pos_eci.r2, timestamp);
    check("ECEF->ECI->ECEF Vx", v2.r0, 1500.0);
    check("ECEF->ECI->ECEF Vy", v2.r1, 2500.0);
    check("ECEF->ECI->ECEF Vz", v2.r2, 800.0);
    if (failures == 0) std::cout << "  \u2705 PASS\n";

    // 2. ENU->ECEF->ENU 往返
    std::cout << "\n2. ENUVel\u2194ECEFVel roundtrip:\n";
    ENUVel2ECEFVel_return v3 = ENUVel2ECEFVel(500, -300, 200, 40.0, 125.0);
    ECEFVel2ENUVel_return v4 = ECEFVel2ENUVel(v3.r0, v3.r1, v3.r2, 40.0, 125.0);
    check("ENU->ECEF->ENU E", v4.r0, 500.0);
    check("ENU->ECEF->ENU N", v4.r1, -300.0);
    check("ENU->ECEF->ENU U", v4.r2, 200.0);
    if (failures == 0) std::cout << "  \u2705 PASS\n";

    // 3. ECI->ECEF->ECI 往返
    std::cout << "\n3. ECIVel\u2194ECEFVel roundtrip:\n";
    double rx_ECI = 6378137, ry_ECI = 0, rz_ECI = 0;
    ECIVel2ECEFVel_return v5 = ECIVel2ECEFVel(3000, -2000, 1000, rx_ECI, ry_ECI, rz_ECI, timestamp);
    // 逆向需要 ECEF 位置: ECI2ECEF
    ECI2ECEF_return pos_ecef = ECI2ECEF(rx_ECI, ry_ECI, rz_ECI, timestamp);
    ECEFVel2ECIVel_return v6 = ECEFVel2ECIVel(v5.r0, v5.r1, v5.r2, pos_ecef.r0, pos_ecef.r1, pos_ecef.r2, timestamp);
    check("ECI->ECEF->ECI Vx", v6.r0, 3000.0);
    check("ECI->ECEF->ECI Vy", v6.r1, -2000.0);
    check("ECI->ECEF->ECI Vz", v6.r2, 1000.0);
    if (failures == 0) std::cout << "  \u2705 PASS\n";

    // 4. AER导数->ENU->AER导数 往返
    std::cout << "\n4. AERDeriv\u2194ENUVel roundtrip:\n";
    AERDeriv2ENUVel_return v7 = AERDeriv2ENUVel(200000, 45, 30, 1000, 0, 0);
    ENUVel2AERDeriv_return v8 = ENUVel2AERDeriv(v7.r0, v7.r1, v7.r2, 200000, 45, 30);
    check("AER->ENU->AER dR", v8.r0, 1000.0);
    check("AER->ENU->AER dAz", v8.r1, 0.0);
    check("AER->ENU->AER dEl", v8.r2, 0.0);
    if (failures == 0) std::cout << "  \u2705 PASS\n";

    // 5. AER一站式->ECI
    std::cout << "\n5. AERDeriv2ECIVel:\n";
    AERDeriv2ECIVel_return v9 = AERDeriv2ECIVel(1200000, 180, 45, 3000, 0, 0.01, 40.0, 116.4, timestamp);
    double speed = std::sqrt(v9.r0*v9.r0 + v9.r1*v9.r1 + v9.r2*v9.r2);
    std::cout << "  ECI speed: " << speed << " m/s\n";
    std::cout << "  \u2705 PASS\n";

    std::cout << "\n=== ";
    if (failures == 0) {
        std::cout << "\u5168\u90e8\u901a\u8fc7!";
    } else {
        std::cout << failures << " \u6d4b\u8bd5\u5931\u8d25!";
    }
    std::cout << " ===\n";

    return failures > 0 ? 1 : 0;
}