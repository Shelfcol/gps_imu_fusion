//
// Created by meng on 2021/2/26.
//

#ifndef GPS_IMU_FUSION_TOOL_H
#define GPS_IMU_FUSION_TOOL_H
#pragma once

#include <ctime>
#include <cstdlib>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <iostream>

//? 坐标系关系
// 右前天->前右地。这里取逆，则为 前右地->右前天。因为数据坐标系是NED，而GPS转换的是右前天，所以需要把IMU相关的数据转换到右前天。
inline void TransformCoordinate(Eigen::Vector3d& vec){
    double kDegree2Radian = M_PI / 180.0;

    Eigen::Quaterniond Q_b_w = Eigen::AngleAxisd(90 * kDegree2Radian, Eigen::Vector3d::UnitZ()) *
                               Eigen::AngleAxisd(0 * kDegree2Radian, Eigen::Vector3d::UnitY()) *
                               Eigen::AngleAxisd(180 * kDegree2Radian, Eigen::Vector3d::UnitX());

    vec = Q_b_w.inverse() * vec;
}

constexpr double kDegree2Radian = M_PI / 180.0;

inline Eigen::Matrix3d BuildSkewMatrix(const Eigen::Vector3d& vec){
    Eigen::Matrix3d matrix;
    matrix << 0.0,     -vec[2],   vec[1],
              vec[2],    0.0,     -vec[0],
              -vec[1],   vec[0],    0.0;

    return matrix;
}

inline double condition_number(const Eigen::MatrixXd &A) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd singularValues;
    singularValues.resize(svd.singularValues().rows(), 1);
    singularValues = svd.singularValues();
    double cond = singularValues(0, 0) / singularValues(singularValues.rows() - 1, 0);
    return std::abs(cond);
}


class TicToc
{
  public:
    TicToc()
    {
        tic();
    }

    TicToc( bool _disp )
    {
        disp_ = _disp;
        tic();
    }

    void tic()
    {
        start = std::chrono::system_clock::now();
    }

    double toc()
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        return elapsed_seconds.count() * 1000;
    }

    void toc( std::string _about_task )
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        double elapsed_ms = elapsed_seconds.count() * 1000;

        if( disp_ )
        {
          std::cout.precision(3); // 10 for sec, 3 for ms 
          std::cout << _about_task << ": " << elapsed_ms << " msec." << std::endl;
        }
    }
  private:
    std::chrono::time_point<std::chrono::system_clock> start, end;
    bool disp_ = true;
};


#endif //GPS_IMU_FUSION_TOOL_H
