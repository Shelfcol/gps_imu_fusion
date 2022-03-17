//
// Created by meng on 2021/2/19.
//
#include "ekf.h"
#include "../3rd/sophus/se3.hpp"

constexpr double kDegree2Radian = M_PI / 180.0;

EKF::EKF(const YAML::Node &node) {
    double gravity = node["earth"]["gravity"].as<double>();
    double earth_rotation_speed = node["earth"]["rotation_speed"].as<double>();
    double cov_prior_posi = node["covariance"]["prior"]["posi"].as<double>();
    double cov_prior_vel = node["covariance"]["prior"]["vel"].as<double>();
    double cov_prior_ori = node["covariance"]["prior"]["ori"].as<double>();
    double cov_prior_epsilon = node["covariance"]["prior"]["epsilon"].as<double>();
    double cov_prior_delta = node["covariance"]["prior"]["delta"].as<double>();
    double cov_measurement_posi = node["covariance"]["measurement"]["posi"].as<double>();
    double cov_process_gyro = node["covariance"]["process"]["gyro"].as<double>();
    double cov_process_accel = node["covariance"]["process"]["accel"].as<double>();
    L_ = node["earth"]["latitude"].as<double>();
    g_ = Eigen::Vector3d(0.0, 0.0, -gravity);


    SetCovarianceP(cov_prior_posi, cov_prior_vel, cov_prior_ori,
                   cov_prior_epsilon, cov_prior_delta);
    SetCovarianceR(cov_measurement_posi);
    SetCovarianceQ(cov_process_gyro, cov_process_accel);
    SetCovarianceW(cov_process_gyro, cov_process_accel);

    X_.setZero(); // 初始化为零矩阵
    X_.block<3,1>(INDEX_STATE_POSI,0) = init_pose_.block<3,1>(0,3);
    X_.block<3,1>(INDEX_STATE_VEL,0) = init_velocity_;

    F_.setZero(); // 初始化为零矩阵
    C_.setIdentity(); // 单位矩阵
    G_.block<3,3>(INDEX_MEASUREMENT_POSI, INDEX_MEASUREMENT_POSI) = Eigen::Matrix3d::Identity();

    F_.block<3,3>(INDEX_STATE_POSI, INDEX_STATE_VEL) = Eigen::Matrix3d::Identity(); 
    B_.setZero();
    gt_.setZero();
    gt_.block<3,1>(INDEX_STATE_VEL,0) = Eigen::Vector3d(0, 0, -gravity);
}

void EKF::SetCovarianceW(double gyro_noise, double accel_noise) {
    W_.setZero();
    W_.block<3,1>(0,0) = Eigen::Vector3d(gyro_noise,gyro_noise,gyro_noise);
    W_.block<3,1>(3,0) = Eigen::Vector3d(accel_noise,accel_noise,accel_noise);
}

void EKF::SetCovarianceQ(double gyro_noise, double accel_noise) {
    Q_.setZero();
    Q_.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * gyro_noise * gyro_noise; // 平方
    Q_.block<3,3>(3,3) = Eigen::Matrix3d::Identity() * accel_noise * accel_noise;
}

void EKF::SetCovarianceR(double posi_noise) {
    R_.setZero();
    R_ = Eigen::Matrix3d::Identity() * posi_noise * posi_noise;
}

// 设置P矩阵
void EKF::SetCovarianceP(double posi_noise, double velo_noise, double ori_noise,
                          double gyro_noise, double accel_noise) {
    P_.setZero();
    P_.block<3,3>(INDEX_STATE_POSI, INDEX_STATE_POSI) = Eigen::Matrix3d::Identity() * posi_noise;
    P_.block<3,3>(INDEX_STATE_VEL, INDEX_STATE_VEL) = Eigen::Matrix3d::Identity() * velo_noise;
    P_.block<3,3>(INDEX_STATE_ORI, INDEX_STATE_ORI) = Eigen::Matrix3d::Identity() * ori_noise;
    P_.block<3,3>(INDEX_STATE_GYRO_BIAS, INDEX_STATE_GYRO_BIAS) = Eigen::Matrix3d::Identity() * gyro_noise;
    P_.block<3,3>(INDEX_STATE_ACC_BIAS, INDEX_STATE_ACC_BIAS) = Eigen::Matrix3d::Identity() * accel_noise;
}

bool EKF::Init(const GPSData &curr_gps_data, const IMUData &curr_imu_data) {
    init_velocity_ = curr_gps_data.true_velocity; // 用真实速度初始化
    velocity_ = init_velocity_;
    // 前右地
    Eigen::Quaterniond Q = Eigen::AngleAxisd(90 * kDegree2Radian, Eigen::Vector3d::UnitZ()) *
                           Eigen::AngleAxisd(0 * kDegree2Radian, Eigen::Vector3d::UnitY()) *
                           Eigen::AngleAxisd(180 * kDegree2Radian, Eigen::Vector3d::UnitX());
    init_pose_.block<3,3>(0,0) = Q.toRotationMatrix();
    pose_ = init_pose_;

    imu_data_buff_.clear(); // 这时EKF中的imu数据
    imu_data_buff_.push_back(curr_imu_data);

    curr_gps_data_ = curr_gps_data;

    return true;
}

// void EKF::GetFGY(TypeMatrixF &F, TypeMatrixG &G, TypeVectorY &Y) {
//     F = Ft_;
//     G = G_;
//     Y = Y_;
// }

bool EKF::Correct(const GPSData &curr_gps_data) {
    curr_gps_data_ = curr_gps_data;

    Y_ = curr_gps_data.position_ned; //Y_measure

    K_ = P_ * G_.transpose() * (G_ * P_ * G_.transpose() + C_ * R_ * C_.transpose()).inverse(); // kalman增益

    P_ = (TypeMatrixP::Identity() - K_ * G_) * P_;
    X_ = X_ + K_ * (Y_ - G_ * X_);

    UpdateState();
    return true;
}

// IMU数据预测
bool EKF::Predict(const IMUData &curr_imu_data) {
    imu_data_buff_.push_back(curr_imu_data);

    UpdateOdomEstimation(); // 更新角度，速度，位置

    double delta_t = curr_imu_data.time - imu_data_buff_.front().time; // dt

    // Eigen::Vector3d curr_accel = pose_.block<3, 3>(0, 0)* curr_imu_data.linear_accel; // 导航坐标系下的加速度
    Eigen::Vector3d curr_accel = curr_imu_data.linear_accel; // 导航坐标系下的加速度
    // Eigen::Vector3d curr_angle_velocity = pose_.block<3,3>(0,0)*curr_imu_data.angle_velocity; // 当前角速度
    Eigen::Vector3d curr_angle_velocity = curr_imu_data.angle_velocity; // 当前角速度
    Eigen::Quaterniond curr_ori = Eigen::Quaterniond(pose_.block<3, 3>(0, 0));
    UpdateEkfState(delta_t, curr_accel,curr_angle_velocity,curr_ori);

    imu_data_buff_.pop_front();
    return true;
}

bool EKF::UpdateEkfState(const double t, const Eigen::Vector3d &accel, const Eigen::Vector3d& curr_angle_velocity,
                        const Eigen::Quaterniond& curr_ori ) {
    double q0 = curr_ori.w();
    double q1 = curr_ori.x();
    double q2 = curr_ori.y();
    double q3 = curr_ori.z();
    double FVq0 = 2 * Eigen::Vector3d(q0,-q3,q2).transpose()*accel;
    double FVq1 = 2 * Eigen::Vector3d(q1,q2,q3).transpose()*accel;
    double FVq2 = 2 * Eigen::Vector3d(-q2,q1,q0).transpose()*accel;
    double FVq3 = 2 * Eigen::Vector3d(-q3,-q0,q1).transpose()*accel;
    Eigen::Matrix<double,3,4> FVq = (Eigen::Matrix<double,3,4>()<<FVq0,FVq1,FVq2,FVq3,
                                                                    -FVq3,-FVq2,FVq1,FVq0,
                                                                    FVq2,-FVq3,-FVq0,FVq1).finished();

    F_.block<3,4>(INDEX_STATE_VEL, INDEX_STATE_ORI) = FVq;
    F_.block<3,3>(INDEX_STATE_VEL, INDEX_STATE_ACC_BIAS) = pose_.block<3,3>(0,0);

    Eigen::Vector3d w = curr_angle_velocity;
    Eigen::Matrix<double,4,4> Fqq = 0.5* (Eigen::Matrix<double,4,4>()<<0,-w.x(),-w.y(),-w.z(),
                                                                        w.x(),0,w.z(),-w.y(),
                                                                        w.y(),-w.z(),0,w.x(),
                                                                        w.z(),w.y(),-w.x(),0).finished();
    F_.block<4,4>(INDEX_STATE_ORI,INDEX_STATE_ORI) = Fqq;

    Eigen::Matrix<double,4,3> Fqkesi  = 0.5 * (Eigen::Matrix<double,4,3>()<<-q1,-q2,-q3,
                                                                        q0,-q3,q2,
                                                                        q3,q0,-q1,
                                                                        -q2,q1,q0).finished();

    F_.block<4,3>(INDEX_STATE_ORI,INDEX_STATE_GYRO_BIAS) = Fqkesi;

    B_.block<3,3>(INDEX_STATE_VEL, 3) = pose_.block<3,3>(0,0);
    B_.block<4,3>(INDEX_STATE_ORI, 0) = Fqkesi;



    TypeMatrixF Fk = TypeMatrixF::Identity() + F_ * t;
    TypeMatrixB Bk = B_ * t;

    Ft_ = F_ * t;

    X_ = Fk * X_ + Bk * W_;// + gt_*t ; 
    P_ = Fk * P_ * Fk.transpose() + Bk * Q_ * Bk.transpose();

    return true;
}

bool EKF::UpdateOdomEstimation() {
    Eigen::Vector3d angular_delta;
    ComputeAngularDelta(angular_delta); // 平均角速度求转动过的角度，以此求delta_R


    Eigen::Matrix3d curr_R, last_R;
    ComputeOrientation(angular_delta, curr_R, last_R);

    Eigen::Vector3d curr_vel, last_vel;
    ComputeVelocity(curr_vel, last_vel, curr_R, last_R);

    ComputePosition(curr_vel, last_vel);

    return true;
}

bool EKF::ComputeAngularDelta(Eigen::Vector3d &angular_delta) {
    IMUData curr_imu_data = imu_data_buff_.at(1);
    IMUData last_imu_data = imu_data_buff_.at(0);

    double delta_t = curr_imu_data.time - last_imu_data.time;

    if (delta_t <= 0){
        return false;
    }

    Eigen::Vector3d curr_angular_vel = curr_imu_data.angle_velocity;

    Eigen::Vector3d last_angular_vel = last_imu_data.angle_velocity;

    // 直接使用last_R来对gyro_bias进行旋转
    Eigen::Matrix3d last_R = pose_.block<3, 3>(0, 0);

    Eigen::Vector3d curr_unbias_angular_vel = curr_angular_vel;
    Eigen::Vector3d last_unbias_angular_vel = last_angular_vel;

    angular_delta = 0.5 * (curr_unbias_angular_vel + last_unbias_angular_vel) * delta_t; // 中值

    return true;
}



bool EKF::ComputeOrientation(const Eigen::Vector3d &angular_delta,
                              Eigen::Matrix3d &curr_R,
                              Eigen::Matrix3d &last_R) {
    Eigen::AngleAxisd angle_axisd(angular_delta.norm(), angular_delta.normalized()); // 轴角公式，前一个为转动角度，后一个为向量。角度转旋转矩阵
    last_R = pose_.block<3, 3>(0, 0);

    curr_R = pose_.block<3, 3>(0, 0) * angle_axisd.toRotationMatrix(); // R*delta_R

    pose_.block<3, 3>(0, 0) = curr_R;

    return true;
}

// 使用去除重力影响和加速度bias的平均加速度计算速度
bool EKF::ComputeVelocity(Eigen::Vector3d &curr_vel, Eigen::Vector3d& last_vel,
                                             const Eigen::Matrix3d &curr_R,
                                             const Eigen::Matrix3d last_R) {
    IMUData curr_imu_data = imu_data_buff_.at(1);
    IMUData last_imu_data = imu_data_buff_.at(0);
    double delta_t = curr_imu_data.time - last_imu_data.time;
    if (delta_t <=0 ){
        return false;
    }

    Eigen::Vector3d curr_accel = curr_imu_data.linear_accel; // body系
    Eigen::Vector3d curr_unbias_accel = GetUnbiasAccel(curr_R * curr_accel);

    Eigen::Vector3d last_accel = last_imu_data.linear_accel;
    Eigen::Vector3d last_unbias_accel = GetUnbiasAccel(last_R * last_accel); // 减去重力影响

    last_vel = velocity_;

    velocity_ += delta_t * 0.5 * (curr_unbias_accel + last_unbias_accel);
    curr_vel = velocity_;

    return true;
}

Eigen::Vector3d EKF::GetUnbiasAccel(const Eigen::Vector3d &accel) {
//    return accel - accel_bias_ + g_; // z方向精度提高很多
    return accel + g_;
}

bool EKF::ComputePosition(const Eigen::Vector3d& curr_vel, const Eigen::Vector3d& last_vel){
    double delta_t = imu_data_buff_.at(1).time - imu_data_buff_.at(0).time;

    pose_.block<3,1>(0,3) += 0.5 * delta_t * (curr_vel + last_vel);

    return true;
}


void EKF::UpdateState() {
    pose_.block<3,1>(0,3) =  X_.block<3,1>(INDEX_STATE_POSI, 0);

    velocity_ =  X_.block<3,1>(INDEX_STATE_VEL, 0);
    Eigen::Quaterniond q;
    q.w()=X_(INDEX_STATE_ORI,0);
    q.x()=X_(INDEX_STATE_ORI+1,0);
    q.y()=X_(INDEX_STATE_ORI+2,0);
    q.z()=X_(INDEX_STATE_ORI+3,0);

    // 修改旋转矩阵
    pose_.block<3,3>(0,0) = q.toRotationMatrix(); // 固定坐标系更新，左乘
    gyro_bias_ = X_.block<3,1>(INDEX_STATE_GYRO_BIAS, 0);
    accel_bias_ = X_.block<3,1>(INDEX_STATE_ACC_BIAS, 0);
}

Eigen::Matrix4d EKF::GetPose() const {
    return pose_;
}