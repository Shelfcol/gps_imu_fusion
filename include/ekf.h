//
// Created by meng on 2021/2/19.
//

#ifndef GPS_IMU_FUSION_EKF_H
#define GPS_IMU_FUSION_EKF_H

#include "imu_data.h"
#include "gps_data.h"
#include "filter_interface.h"

#include <deque>
#include <yaml-cpp/yaml.h>
#include <eigen3/Eigen/Dense>

class EKF :public FilterInterface {
public:
    EKF(const YAML::Node &node);

    /*!
     * 用于EKF滤波器的初始化，设置初始位姿，初始速度
     * @param curr_gps_data 与imu时间同步的gps数据
     * @param curr_imu_data 与gps时间同步的imu数据
     * @return
     */
    bool Init(const GPSData &curr_gps_data, const IMUData &curr_imu_data);

    /*!
     * 滤波器的预测，对应卡尔曼滤波器的前两个公式
     * @param curr_imu_data
     * @return
     */
    bool Predict(const IMUData &curr_imu_data);

    /*!
     * 滤波器的矫正，对应卡尔曼滤波器的后三个公式
     * @param curr_gps_data
     * @return
     */
    bool Correct(const GPSData &curr_gps_data);

    Eigen::Matrix4d GetPose() const;

    Eigen::Vector3d GetVelocity(){
        return velocity_;
    }

private:
    void SetCovarianceQ(double gyro_noise_cov, double accel_noise_cov);
    void SetCovarianceW(double gyro_noise_cov, double accel_noise_cov);

    void SetCovarianceR(double posi_noise_cov);

    void SetCovarianceP(double posi_noise, double velo_noise, double ori_noise,
                        double gyro_noise, double accel_noise);

    /*!
     * 通过IMU计算位姿和速度
     * @return
     */
    bool UpdateOdomEstimation();

    bool UpdateEkfState(const double t, const Eigen::Vector3d &accel, const Eigen::Vector3d& curr_angle_velocity,
                        const Eigen::Quaterniond& curr_ori );

    bool ComputeAngularDelta(Eigen::Vector3d &angular_delta);

    /*!
     * 通过IMU计算当前姿态
     * @param angular_delta
     * @param R_nm_nm_1
     * @param curr_R
     * @param last_R
     * @return
     */
    bool ComputeOrientation(const Eigen::Vector3d &angular_delta,
                            Eigen::Matrix3d &curr_R,
                            Eigen::Matrix3d &last_R);

    bool ComputeVelocity(Eigen::Vector3d &curr_vel,
                         Eigen::Vector3d &last_vel,
                         const Eigen::Matrix3d &curr_R,
                         const Eigen::Matrix3d last_R);

    Eigen::Vector3d GetUnbiasAccel(const Eigen::Vector3d &accel);

    /*!
     * 通过imu计算当前位移
     * @param curr_vel
     * @param last_vel
     * @return
     */
    bool ComputePosition(const Eigen::Vector3d& curr_vel, const Eigen::Vector3d& last_vel);

    void UpdateState();


private:
    static const unsigned int DIM_STATE = 16;
    static const unsigned int DIM_STATE_NOISE = 6; // 噪声只有6维，陀螺仪和加速度计的bias
    static const unsigned int DIM_MEASUREMENT = 3; // 观测向量只有3维
    static const unsigned int DIM_MEASUREMENT_NOISE = 3; // 观测噪声

    static const unsigned int INDEX_STATE_POSI = 0; // 位置
    static const unsigned int INDEX_STATE_VEL = 3; // 速度
    static const unsigned int INDEX_STATE_ORI = 6; // 角度
    static const unsigned int INDEX_STATE_GYRO_BIAS = 10; // 陀螺仪bias
    static const unsigned int INDEX_STATE_ACC_BIAS = 13; // 加速度计bias
    static const unsigned int INDEX_MEASUREMENT_POSI = 0;

    typedef typename Eigen::Matrix<double, DIM_STATE, 1> TypeVectorX; // 状态向量
    typedef typename Eigen::Matrix<double, DIM_MEASUREMENT, 1> TypeVectorY; // 观测向量 GPS的位置
    typedef typename Eigen::Matrix<double, DIM_STATE, DIM_STATE> TypeMatrixF; //15*15
    typedef typename Eigen::Matrix<double, DIM_STATE, DIM_STATE_NOISE> TypeMatrixB; //15*6
    typedef typename Eigen::Matrix<double, DIM_STATE_NOISE, 1> TypeMatrixW; //6*1 陀螺仪和加速度计噪声
    typedef typename Eigen::Matrix<double, DIM_STATE_NOISE, DIM_STATE_NOISE> TypeMatrixQ;
    typedef typename Eigen::Matrix<double, DIM_STATE, DIM_STATE> TypeMatrixP;
    typedef typename Eigen::Matrix<double, DIM_STATE, DIM_MEASUREMENT> TypeMatrixK;
    typedef typename Eigen::Matrix<double, DIM_MEASUREMENT_NOISE, DIM_MEASUREMENT_NOISE> TypeMatrixC;
    typedef typename Eigen::Matrix<double, DIM_MEASUREMENT, DIM_STATE> TypeMatrixG;
    typedef typename Eigen::Matrix<double, DIM_MEASUREMENT, DIM_MEASUREMENT> TypeMatrixR;

    TypeVectorX X_;
    TypeVectorY Y_;
    TypeMatrixF F_;
    TypeMatrixB B_;
    TypeMatrixW W_;
    TypeMatrixQ Q_;
    TypeMatrixP P_;
    TypeMatrixK K_;
    TypeMatrixC C_;
    TypeMatrixG G_;
    TypeMatrixC R_;

    TypeMatrixF Ft_;
    TypeVectorX gt_; // 加速度项

    Eigen::Vector3d init_velocity_ = Eigen::Vector3d::Zero();
    Eigen::Vector3d velocity_ = Eigen::Vector3d::Zero();
    Eigen::Matrix4d init_pose_ = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d pose_ = Eigen::Matrix4d::Identity();

    Eigen::Vector3d gyro_bias_ = Eigen::Vector3d::Zero();
    Eigen::Vector3d accel_bias_ = Eigen::Vector3d::Zero();

    Eigen::Vector3d g_;//重力加速度
    Eigen::Vector3d w_;//地球自传角速度

    GPSData curr_gps_data_;

    double L_ = 0.0;//纬度

    std::deque<IMUData> imu_data_buff_; // 只保存两个IMU数据

public:
    // void GetFGY(TypeMatrixF& F,TypeMatrixG& G, TypeVectorY & Y);
};

#endif //GPS_IMU_FUSION_EKF_H