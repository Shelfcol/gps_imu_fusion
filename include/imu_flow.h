//
// Created by meng on 2021/2/19.
//

#ifndef GPS_IMU_FUSION_IMU_FLOW_H
#define GPS_IMU_FUSION_IMU_FLOW_H

#include "imu_data.h"

#include <vector>
#include <deque>

class IMUFlow{
public:
    IMUFlow() = default;

    // 一个用deque，一个用vector
    static bool ReadIMUData(const std::string& path,
                            std::vector<IMUData>& imu_data_buff,
                            const int skip_rows = 1);

    static bool ReadIMUData(const std::string& path,
                            std::deque<IMUData>& imu_data_buff,
                            const int skip_rows = 1);
};

#endif //GPS_IMU_FUSION_IMU_FLOW_H
