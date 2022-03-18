//
// Created by meng on 2021/2/19.
//

#ifndef GPS_IMU_FUSION_GPS_DATA_H
#define GPS_IMU_FUSION_GPS_DATA_H

#include <eigen3/Eigen/Core>
#include "../3rd/GeographicLib/include/Geocentric/LocalCartesian.hpp"
class GPSData{
public:
    GPSData() = default;

    Eigen::Vector3d position_lla = Eigen::Vector3d::Zero();//LLA
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero(); // 速度 NED
    Eigen::Vector3d position_ned = Eigen::Vector3d::Zero(); // 位置 NED

    Eigen::Vector3d true_velocity = Eigen::Vector3d::Zero();
    Eigen::Vector3d true_position_lla = Eigen::Vector3d::Zero();
    double time = 0.0;
};

#endif //GPS_IMU_FUSION_GPS_DATA_H
