//
// Created by meng on 2021/2/19.
//
#include "ukf.h"
#include "../3rd/sophus/se3.hpp"
// https://blog.csdn.net/l2014010671/article/details/93305871
// https://zhuanlan.zhihu.com/p/482392082
UKF::UKF(const YAML::Node &node) {

    // 判断使用哪一组数据
    std::string  data_path = node["data_path"].as<std::string>();
    std::string cov_node_string;
    if(data_path == "/data/raw_data") {
        cov_node_string = "covariance";
    } else if(data_path == "/data/raw_data1"){
        cov_node_string = "covariance1";
    } else {
        printf("no corres covariance");
        exit(0);
    }

    double gravity = node["earth"]["gravity"].as<double>();

    // 这里直接给的平方项
    double cov_prior_posi = node["UKF"][cov_node_string]["prior"]["posi2"].as<double>();
    double cov_prior_vel = node["UKF"][cov_node_string]["prior"]["vel2"].as<double>();
    double cov_prior_ori = node["UKF"][cov_node_string]["prior"]["ori2"].as<double>();
    double cov_prior_accel = node["UKF"][cov_node_string]["prior"]["accel_delta2"].as<double>();
    double cov_prior_gyro = node["UKF"][cov_node_string]["prior"]["gyro_delta2"].as<double>();
    double cov_prior_g = node["UKF"][cov_node_string]["prior"]["g_delta2"].as<double>();

    // 陀螺仪和加速度计的噪声和零偏噪声方差，这里给
    double cov_acc_delta = node["UKF"][cov_node_string]["IMU_noise"]["acc_delta2"].as<double>();
    double cov_gyro_delta = node["UKF"][cov_node_string]["IMU_noise"]["gyro_delta2"].as<double>();
    double cov_accel_bias_delta = node["UKF"][cov_node_string]["IMU_noise"]["accel_bias_delta2"].as<double>();
    double cov_gyro_bias_delta = node["UKF"][cov_node_string]["IMU_noise"]["gyro_bias_delta2"].as<double>();

    double measurement_posi = node["UKF"][cov_node_string]["measurement"]["posi"].as<double>();
    g_ = Eigen::Vector3d(0.0, 0.0, -gravity);

    n_a = DIM_AUG_STATE; // 增广状态量维度
    sigma_points_num_ = 2*n_a+1;
    SetCovarianceP(cov_prior_posi, cov_prior_vel, cov_prior_ori,
                   cov_prior_accel, cov_prior_gyro);
    SetCovarianceQ(cov_acc_delta, cov_gyro_delta,cov_accel_bias_delta,cov_gyro_bias_delta);
    SetCovarianceR(measurement_posi);
    calcWeights();
}

void UKF::calcWeights()
{   
    double alpha = 0.1;  // 1 × 10−4 ≤ α ≤ 1
    double beta = 2.;
    double kappa = 0;  // 0. or 3 - N;
    double alpha2 = alpha * alpha;

    double lambda = alpha2 * (n_a + kappa) - n_a;
    scale = sqrt(lambda+n_a);
    // 每个sigma点的权重
    weights_m_[0]=lambda/(lambda+n_a); 
    weights_c_[0]=lambda/(lambda+n_a)+(1 - alpha2 + beta); //! 保证P正定
    for(int i=1;i<weights_c_.size();++i)
    {
        weights_m_[i] = 0.5/(lambda+n_a);
        weights_c_[i] = 0.5/(lambda+n_a);
    }
}

// 设置P矩阵
void UKF::SetCovarianceP(double posi_noise, double velo_noise, double ori_noise,
                          double accel_noise, double gyro_noise) {
    P_.setZero();
    P_.block<3,3>(INDEX_STATE_POSI, INDEX_STATE_POSI) = Eigen::Matrix3d::Identity() * posi_noise;
    P_.block<3,3>(INDEX_STATE_VEL, INDEX_STATE_VEL) = Eigen::Matrix3d::Identity() * velo_noise;
    P_.block<3,3>(INDEX_STATE_ORI, INDEX_STATE_ORI) = Eigen::Matrix3d::Identity() * ori_noise;
    P_.block<3,3>(INDEX_STATE_ACC_BIAS, INDEX_STATE_ACC_BIAS) = Eigen::Matrix3d::Identity() * accel_noise;
    P_.block<3,3>(INDEX_STATE_GYRO_BIAS, INDEX_STATE_GYRO_BIAS) = Eigen::Matrix3d::Identity() * gyro_noise;
}

void UKF::SetCovarianceQ(const double cov_acc_delta, const double cov_gyro_delta, const double cov_accel_bias_delta, const double cov_gyro_bias_delta) {
    Q_.setZero();
    Q_.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * cov_acc_delta; // 平方
    Q_.block<3,3>(3,3) = Eigen::Matrix3d::Identity() * cov_gyro_delta;
    Q_.block<3,3>(6,6) = Eigen::Matrix3d::Identity() * cov_accel_bias_delta;
    Q_.block<3,3>(9,9) = Eigen::Matrix3d::Identity() * cov_gyro_bias_delta;
}

void UKF::SetCovarianceR(double measurement_noise) {
    R_ = Eigen::Matrix3d::Identity()*measurement_noise*measurement_noise;
}   

bool UKF::Init(const GPSData &curr_gps_data, const IMUData &curr_imu_data) {
    // 前右地
    Eigen::Quaterniond Q = Eigen::AngleAxisd(90 * kDegree2Radian, Eigen::Vector3d::UnitZ()) *
                           Eigen::AngleAxisd(0 * kDegree2Radian, Eigen::Vector3d::UnitY()) *
                           Eigen::AngleAxisd(180 * kDegree2Radian, Eigen::Vector3d::UnitX());
    // 状态量的初始化
    X_.setZero();
    X_.block<3,1>(INDEX_STATE_POSI,0) = curr_gps_data.position_ned;;
    X_.block<3,1>(INDEX_STATE_VEL,0)  = curr_gps_data.true_velocity;
    X_(INDEX_STATE_ORI,  0) = Q.w();
    X_(INDEX_STATE_ORI+1,0) = Q.x();
    X_(INDEX_STATE_ORI+2,0) = Q.y();
    X_(INDEX_STATE_ORI+3,0) = Q.z();

    imu_data_buff_.clear(); // 这时ESKF中的imu数据
    imu_data_buff_.push_back(curr_imu_data);
    return true;
}

Eigen::MatrixXd UKF::GetRootSquare(const Eigen::MatrixXd& M, const bool is_SVD)
{
    Eigen::MatrixXd  A;
    if(!is_SVD){ // 0.2ms
        A = M.llt().matrixL();
    }
    else{
        // for not-PSD cov matrix  4ms
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const auto &S = svd.matrixU().inverse() * M * svd.matrixV().transpose().inverse();  // S = U^-1 * A * VT * -1
        const Eigen::MatrixXd &Smat = S.matrix();
        const Eigen::MatrixXd &S_sqr = Smat.llt().matrixL();
        A = svd.matrixU().matrix() * S_sqr;  // U * sqrt(S)
    }
    return A;
}

// 使用状态量和方差生成sigma点。这里可能是增广矩阵和非增广矩阵，行数可能不包含噪声项
Eigen::MatrixXd UKF::GenerateSigmaMatrix(const Eigen::MatrixXd& X, const Eigen::MatrixXd& P)
{
    bool is_SVD_P = true;
    TicToc root_tic;
    Eigen::MatrixXd A = GetRootSquare(P, is_SVD_P);
    // root_tic.toc("GetRootSquare");
    int rows = X.rows(); // 状态矩阵维度
    Eigen::MatrixXd  SigMat(rows, 2*n_a+1);
    SigMat.col(0) = X;
    TicToc tic;
    for(int i=0;i<n_a;++i)
    {
        SigMat.col(i+1)     = X - scale * A.col(i);      
        SigMat.col(i+1+n_a) = X + scale * A.col(i);      
    }
    // tic.toc("SigMat");
    return SigMat;
}


/*!
    * 通过IMU计算位姿和速度
    * @return
    */
bool UKF::GenerateSigmaPoints() // 生成sigma点，注意里面的P_aug_中的Q部分需要重置为Q，因为这部分本来是不会改变的
{
    // 重置P_aug_中的Q，因为传感器的噪声可以看作固定
    P_aug_.setZero();
    P_aug_.topLeftCorner(DIM_STATE,DIM_STATE) = P_;
    P_aug_.bottomRightCorner(DIM_IMU_NOISE,DIM_IMU_NOISE) = Q_;

    X_aug_.setZero(); // 噪声项对应的状态量初始为0
    X_aug_.topLeftCorner(DIM_STATE,1) = X_;
    TicToc tic;
    X_sig_aug_ = GenerateSigmaMatrix(X_aug_, P_aug_);
    // tic.toc("GenerateSigmaMatrix");
    return true;
}

// 不带增广的噪声项 ，求解X_sig_aug_pred_
bool UKF::PredictSigmaPoints(const Eigen::Vector3d& a_m, const Eigen::Vector3d& w_m, const double dt) // 通过X_sig_aug_和状态矩阵传递方程预测sigma点
{
    omp_set_num_threads(20); //设置线程的个数
	#pragma omp parallel for // 加了之后耗时由3ms变为0.6ms
    for(int i=0; i<sigma_points_num_; ++i)
    {
        Eigen::Vector3d p           = X_sig_aug_.col(i).block<3,1>(INDEX_STATE_POSI,0);
        Eigen::Vector3d v           = X_sig_aug_.col(i).block<3,1>(INDEX_STATE_VEL,0);
        Eigen::Quaterniond q;
        q.w()                       = X_sig_aug_.col(i)(INDEX_STATE_ORI  ,0);
        q.x()                       = X_sig_aug_.col(i)(INDEX_STATE_ORI+1,0);
        q.y()                       = X_sig_aug_.col(i)(INDEX_STATE_ORI+2,0);
        q.z()                       = X_sig_aug_.col(i)(INDEX_STATE_ORI+3,0);
        q.normalize();
        Eigen::Matrix3d R           = q.toRotationMatrix();
        Eigen::Vector3d a_b         = X_sig_aug_.col(i).block<3,1>(INDEX_STATE_ACC_BIAS,0);
        Eigen::Vector3d w_b         = X_sig_aug_.col(i).block<3,1>(INDEX_STATE_GYRO_BIAS,0);
        Eigen::Vector3d acc_n       = X_sig_aug_.col(i).block<3,1>(INDEX_STATE_ACC_NOISE,0);
        Eigen::Vector3d gyro_n      = X_sig_aug_.col(i).block<3,1>(INDEX_STATE_GYRO_NOISE,0);
        Eigen::Vector3d acc_bias_n  = X_sig_aug_.col(i).block<3,1>(INDEX_STATE_ACC_BIAS_NOISE,0);
        Eigen::Vector3d gyro_bias_n = X_sig_aug_.col(i).block<3,1>(INDEX_STATE_GYRO_BIAS_NOISE,0);

        // 为了与Quternion kinematics for error state kalman filter中的变量保持一致
        auto a_n = acc_n;
        auto w_n = gyro_n;
        auto a_w = acc_bias_n;
        auto w_w = gyro_bias_n;

        //  process过程
        p += v*dt+0.5*(R*(a_m-a_b-a_n)+g_)*dt*dt;
        v += (R*(a_m-a_b-a_n)+g_)*dt;
        Eigen::Vector3d d_theta = (w_m-w_b-w_n)*dt;
        Eigen::AngleAxisd dq(d_theta.norm(),d_theta.normalized());
        q = q*Eigen::Quaterniond(dq);
        a_b += a_w*dt;
        w_b += w_w*dt;

        X_sig_.col(i).block<3,1>(INDEX_STATE_POSI,0) = p;
        X_sig_.col(i).block<3,1>(INDEX_STATE_VEL,0) = v;
        X_sig_.col(i)(INDEX_STATE_ORI,  0) = q.w();
        X_sig_.col(i)(INDEX_STATE_ORI+1,0) = q.x();
        X_sig_.col(i)(INDEX_STATE_ORI+2,0) = q.y();
        X_sig_.col(i)(INDEX_STATE_ORI+3,0) = q.z();
        X_sig_.col(i).block<3,1>(INDEX_STATE_ACC_BIAS,0) = a_b;
        X_sig_.col(i).block<3,1>(INDEX_STATE_GYRO_BIAS,0) = w_b;
    }

    return true;
}

// 更新X_, P_
bool UKF::PredictStateMeanAndCovariance() // 通过预测的sigma点和对应的权重预测均值和方差，同时更新X_和P_，便于观测数据来
{
    // 预测得到的sigma点均值
    X_ = X_sig_ * weights_m_;
    // 预测得到的P_
    P_.setZero();

	// omp_set_num_threads(20); //设置线程的个数
	// #pragma omp parallel for
    for(int i=0;i<sigma_points_num_;++i) // 1.33ms
    {
        TypeVectorX diff = X_sig_.col(i)-X_;
        P_ += weights_c_(i)*diff*diff.transpose();
    }
    return true;
}

bool UKF::PredictMeasurementAndNoise() // 根据观测矩阵和预测的均值求解预测的观测值和预测的观测噪声
{
    // 生成sigma点
    // Eigen::MatrixXd  X_sig = GenerateSigmaMatrix(X_,  P_); // 直接使用预测的矩阵
    omp_set_num_threads(20); //设置线程的个数
	#pragma omp parallel for
    for(int i=0; i<sigma_points_num_; ++i) 
    {
        Z_sig_.col(i)(0) = X_sig_.col(i)(0);
        Z_sig_.col(i)(1) = X_sig_.col(i)(1);
        Z_sig_.col(i)(2) = X_sig_.col(i)(2);
    }
    Z_ = Z_sig_ * weights_m_;

    S_ = R_;
    for(int i=0; i<sigma_points_num_; ++i)
    {
        Eigen::VectorXd diff = Z_sig_.col(i)-Z_;
        S_ += weights_c_(i)*diff*diff.transpose();
    }
    return true;
}

bool UKF::UpdateState(const Eigen::Vector3d measurement)
{
    T_.setZero();
    for(int i=0;i<sigma_points_num_;++i)
    {
        T_ += weights_c_(i)*(X_sig_.col(i)-X_)*(Z_sig_.col(i)-Z_).transpose();
    }
    K_ = T_*S_.inverse();
    X_ = X_ + K_*(measurement-Z_);
    P_ = P_ - K_*S_*K_.transpose();

    // P_ = 0.5 * (P_ + P_.transpose());

    // condition number
    // 解决：因观测误差较大，使P负定，致使后面P的Cholesky分解失败出现NaN，导致滤波器发散
    // {
    //   double cond_num = condition_number(P_);
    //   std::cout << "cond num of P: " << cond_num << std::endl;
    //   if (cond_num > 1e5) P_ = P_.diagonal().asDiagonal();
    // }

    return true;
}


// void UKF::GetFGY(TypeMatrixF &F, TypeMatrixG &G, TypeVectorY &Y) {
//     F = Ft_;
//     G = G_;
//     Y = Y_;
// }

// IMU数据预测
bool UKF::Predict(const IMUData &curr_imu_data) {
    imu_data_buff_.push_back(curr_imu_data);

    IMUData last_imu_data = imu_data_buff_.at(0);

    double delta_t = curr_imu_data.time - last_imu_data.time;
    Eigen::Vector3d a_m = curr_imu_data.linear_accel;
    Eigen::Vector3d w_m = curr_imu_data.angle_velocity;
    TicToc tic;
    tic.tic();
    GenerateSigmaPoints(); // 生成sigma点，注意里面的P_aug_中的Q部分需要重置为Q，因为这部分本来是不会改变的  4ms
    // tic.toc("GenerateSigmaPoints");
    tic.tic();
    PredictSigmaPoints( a_m, w_m, delta_t); // 通过X_sig_aug_和状态矩阵传递方程预测sigma点 3ms
    // tic.toc("PredictSigmaPoints");
    tic.tic();
    PredictStateMeanAndCovariance(); // 通过预测的sigma点和对应的权重预测均值和方差，同时更新X_和P_，便于观测数据来  1ms
    // tic.toc("PredictStateMeanAndCovariance");

    imu_data_buff_.pop_front();
    // pred_tic.toc("Predict");
    return true;
}

bool UKF::Correct(const GPSData &curr_gps_data)
{
    TicToc cor_tic;
    Eigen::Vector3d measurement = curr_gps_data.position_ned;
    PredictMeasurementAndNoise(); // 根据观测矩阵和预测的均值求解预测的观测值和预测的观测噪声
    UpdateState(measurement);
    // cor_tic.toc("correct");
    std::cout<<"a_b"<<X_.block<3,1>(INDEX_STATE_ACC_BIAS,0)<<std::endl;
    std::cout<<"w_b"<<X_.block<3,1>(INDEX_STATE_GYRO_BIAS,0)<<std::endl;
    return true;
}

Eigen::Matrix4d UKF::GetPose() const {
    Eigen::Quaterniond q;
    q.w()=X_(INDEX_STATE_ORI,0);
    q.x()=X_(INDEX_STATE_ORI+1,0);
    q.y()=X_(INDEX_STATE_ORI+2,0);
    q.z()=X_(INDEX_STATE_ORI+3,0);

    Eigen::Matrix4d pose;
    pose.block<3,3>(0,0) = q.toRotationMatrix();
    pose.block<3,1>(0,3) = X_.block<3,1>(INDEX_STATE_POSI,0);
    return pose;
}