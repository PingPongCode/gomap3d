#include <iostream>
#include "libmain.h"

inline GoString toGoString(const std::string& s) {
    return {s.c_str(), static_cast<long>(s.length())};
}

int main() {
    // 测试ENU2AER
    ENU2AER_return aer = ENU2AER(10.0, 20.0, 30.0);
    std::cout << "ENU2AER Result: "
              << "Azimuth: " << aer.az << " "
              << "Elevation: " << aer.el << " "
              << "Range: " << aer.srange << "m\n";

    // 测试ECEF2Geodetic
    ECEF2Geodetic_return geo =ECEF2Geodetic(-2700226.9,  -4292413.9,  3855273.8, toGoString("wgs84"));
    std::cout << "\nECEF2Geodetic Result: "
              << "Lat: " << geo.lat << " "
              << "Lon: " << geo.lon << " "
              << "Alt: " << geo.alt << "m\n";

    // 测试ECI/ECEF转换
    GoFloat64 timestamp =  1717000000; //Unix时间戳 秒 浮点数
    ECI2ECEF_return ecef = ECI2ECEF(6500000, 100000, 300000, timestamp);
    std::cout << "\nECI2ECEF Result: "
              << "X: " << ecef.xEcef << "m "
              << "Y: " << ecef.yEcef << "m "
              << "Z: " << ecef.zEcef << "m\n";

    return 0;
}