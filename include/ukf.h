//
// Created by meng on 2021/2/19.
//

#ifndef GPS_IMU_FUSION_UKF_H
#define GPS_IMU_FUSION_UKF_H
#include <omp.h>
#include "imu_data.h"
#include "gps_data.h"
#include "filter_interface.h"
#include "tool.h"
#include <deque>
#include <yaml-cpp/yaml.h>
#include <eigen3/Eigen/Dense>

// quaternion kinematics for error state kalman filter
class UKF :public FilterInterface {
public:
    UKF(const YAML::Node &node);

    /*!
     * 用于ESKF滤波器的初始化，设置初始位姿，初始速度
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
        Eigen::Vector3d vel = X_.block<3,1>(INDEX_STATE_VEL,0);
        return vel;
    }

private:

    // 设置初始P矩阵
    void SetCovarianceP(double posi_noise, double velo_noise, double ori_noise,
                        double accel_noise, double gyro_noise);
    // 观测噪声
    void SetCovarianceR(double measurement_noise);
    
    // 加速度计，陀螺仪噪声，加速度计，陀螺仪零偏噪声
    void SetCovarianceQ(const double cov_acc_delta, const double cov_gyro_delta, 
                        const double cov_accel_bias_delta, const double cov_gyro_bias_delta);

    void calcWeights();
    Eigen::MatrixXd   GetRootSquare(const Eigen::MatrixXd& M, const bool is_SVD);
    // 使用状态量和方差生成sigma点。这里可能是增广矩阵和非增广矩阵，行数可能不包含噪声项
    Eigen::MatrixXd GenerateSigmaMatrix(const Eigen::MatrixXd& X, const Eigen::MatrixXd& P);
    /*!
     * 通过IMU计算位姿和速度
     * @return
     */
    bool GenerateSigmaPoints(); // 生成sigma点，注意里面的P_aug_中的Q部分需要重置为Q，因为这部分本来是不会改变的
    bool PredictSigmaPoints(const Eigen::Vector3d& a_m, const Eigen::Vector3d& w_m, const double dt); // 通过X_sig_aug_和状态矩阵传递方程预测sigma点
    bool PredictStateMeanAndCovariance(); // 通过预测的sigma点和对应的权重预测均值和方差，同时更新X_和P_，便于观测数据来

    // 更新
    bool PredictMeasurementAndNoise(); // 根据观测矩阵和预测的均值求解预测的观测值和预测的观测噪声
    bool UpdateState(const Eigen::Vector3d measurement);


private:
    static const unsigned int DIM_STATE = 16; // 状态向量 p,v,q,a_bias,w_bias
    static const unsigned int DIM_IMU_NOISE = 12; // 噪声只有12维， a_n,w_n,a_w,w_w: 加速度角速度噪声，加速度角速度bias噪声
    static const unsigned int DIM_AUG_STATE {DIM_STATE+DIM_IMU_NOISE}; // 增广状态量维度，因为生成sigma点时噪声也是参与非线性传播
    static const unsigned int DIM_SIGMA {2*DIM_AUG_STATE+1}; // sigma矩阵的维度
    static const unsigned int DIM_MEASUREMENT = 3; // 观测向量只有3维
    static const unsigned int DIM_MEASUREMENT_NOISE = 3; // 观测噪声

    static const unsigned int INDEX_STATE_POSI = 0; // 位置
    static const unsigned int INDEX_STATE_VEL = 3; // 速度
    static const unsigned int INDEX_STATE_ORI = 6; // 角度
    static const unsigned int INDEX_STATE_ACC_BIAS = 10; // 加速度计bias
    static const unsigned int INDEX_STATE_GYRO_BIAS = 13; // 陀螺仪bias
    static const unsigned int INDEX_STATE_ACC_NOISE = 16; // 加速度计噪声
    static const unsigned int INDEX_STATE_GYRO_NOISE = 19; // 陀螺仪噪声
    static const unsigned int INDEX_STATE_ACC_BIAS_NOISE = 22; // 加速度计零偏噪声
    static const unsigned int INDEX_STATE_GYRO_BIAS_NOISE = 25; // 陀螺仪零偏噪声
    static const unsigned int INDEX_MEASUREMENT_POSI = 0;

    typedef typename Eigen::Matrix<double, DIM_STATE, 1> TypeVectorX; // 状态向量
    typedef typename Eigen::Matrix<double, DIM_AUG_STATE, 1> TypeVectorXAug; // 增广状态向量
    typedef typename Eigen::Matrix<double, DIM_AUG_STATE, DIM_SIGMA> TypeMatrixXSigAug; // sigma点形成的矩阵
    typedef typename Eigen::Matrix<double, DIM_STATE, DIM_SIGMA> TypeMatrixXSigAugWithoutNoise; // sigma点形成的矩阵，不包含增广的噪声项
    typedef typename Eigen::Matrix<double, DIM_SIGMA, 1> TypeVectorWeight; // sigma点权重矩阵
    typedef typename Eigen::Matrix<double, DIM_STATE, DIM_STATE> TypeMatrixP;
    typedef typename Eigen::Matrix<double, DIM_AUG_STATE, DIM_AUG_STATE> TypeMatrixPAug; //增广的P矩阵
    typedef typename Eigen::Matrix<double, DIM_AUG_STATE, DIM_SIGMA> TypeMatrixPSigAug; //增广的Sigma点得到的P矩阵
    typedef typename Eigen::Matrix<double, DIM_STATE, DIM_SIGMA> TypeMatrixPSigAugWithoutNoise; //Sigma点得到的P矩阵，不包含增广的噪声项

    typedef typename Eigen::Matrix<double, DIM_MEASUREMENT, 1> TypeVectorZ; // 观测向量 GPS的位置
    typedef typename Eigen::Matrix<double, DIM_MEASUREMENT, DIM_SIGMA> TypeMatrixZSig; // 预测的观测矩阵

    typedef typename Eigen::Matrix<double, DIM_MEASUREMENT, DIM_MEASUREMENT> TypeMatrixR; // 观测噪声
    typedef typename Eigen::Matrix<double, DIM_IMU_NOISE, DIM_IMU_NOISE> TypeMatrixQ;
    typedef typename Eigen::Matrix<double, DIM_STATE, DIM_MEASUREMENT> TypeMatrixK;
    typedef typename Eigen::Matrix<double, DIM_STATE, DIM_MEASUREMENT> TypeMatrixT;
    typedef typename Eigen::Matrix<double, DIM_MEASUREMENT, DIM_MEASUREMENT> TypeMatrixS;

    TypeVectorX X_; // 状态量
    TypeVectorXAug X_aug_; // 增广的状态量
    TypeMatrixXSigAugWithoutNoise X_sig_; // 不包含增广的噪声项状态量，用于观测更新
    TypeMatrixXSigAug X_sig_aug_; // 增广的sigma点组成的矩阵

    TypeVectorWeight weights_m_; // sigma点权重矩阵
    TypeVectorWeight weights_c_; // sigma点方差矩阵

    TypeMatrixP P_; // 方差
    TypeMatrixPAug P_aug_; // 增广的噪声矩阵

    TypeVectorZ Z_; // 预测的测量矩阵
    TypeMatrixZSig Z_sig_; // 预测的状态量得到的观测的sigma矩阵


    TypeMatrixR R_; // 观测噪声 
    TypeMatrixQ Q_; // 噪声矩阵， 噪声只有12维， a_n,w_n,a_w,w_w: 加速度角速度噪声，加速度角速度bias噪声
    TypeMatrixK K_; // 卡尔曼增益
    TypeMatrixT T_; // cross-correlation
    TypeMatrixS S_; // 预测的观测噪声

    int n_a; // 增广状态量维度
    int sigma_points_num_; // sigma点个数
    double scale; // 生成sigma点时的尺度信息

    // // 状态量
    // Eigen::Vector3d pose_ = Eigen::Vector3d::Zero();
    // Eigen::Vector3d velocity_ = Eigen::Vector3d::Zero();
    // Eigen::Quaterniond q_;
    // Eigen::Vector3d accel_bias_ = Eigen::Vector3d::Zero();
    // Eigen::Vector3d gyro_bias_ = Eigen::Vector3d::Zero();
    Eigen::Vector3d g_;//重力加速度

    std::deque<IMUData> imu_data_buff_; // 只保存两个IMU数据

public:
    // void GetFGY(TypeMatrixF& F,TypeMatrixG& G, TypeVectorY & Y);
};

#endif //GPS_IMU_FUSION_ESKF_H