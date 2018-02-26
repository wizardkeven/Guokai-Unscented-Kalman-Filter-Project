#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//globle constant definition
#define N_X 5
#define N_AUG 7
#define PARA_LAMBDA 3
#define AUG_SIG_NUM 2*N_AUG+1
#define NIS_RADAR_FILENAME "NIS_OF_RADAR.csv"
#define NIS_LIDAR_FILENAME "NIS_OF_LIDAR.csv"

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except duridng init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.4;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  n_x_ = N_X;

  n_aug_ = N_AUG;

  lambda_ = PARA_LAMBDA - n_aug_;

  //instanciate state vectors
  VectorXd x_ = VectorXd(n_x_);
  MatrixXd P_ = MatrixXd(n_x_ , n_x_);
  Xsig_pred_ = MatrixXd(n_x_, AUG_SIG_NUM);

  //initiate weights
  weights_ = VectorXd(AUG_SIG_NUM);
  // set weights
  double weight_0 = lambda_/(lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<AUG_SIG_NUM; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_ + lambda_);
    weights_(i) = weight;
  }

  lidarNisRecorder_.open(NIS_LIDAR_FILENAME, std::ofstream::out | std::ofstream::trunc);
  radarNisRecorder_.open(NIS_RADAR_FILENAME, std::ofstream::out | std::ofstream::trunc);
}

UKF::~UKF() {
  std::cout<<"Size of lidar vector: "<<lidarNIS_.size()<<"\n";
  if(!lidarNisRecorder_.is_open()){
    lidarNisRecorder_.open(NIS_LIDAR_FILENAME);
  }
  if(lidarNisRecorder_.is_open()){
    for(double ni : lidarNIS_){
      lidarNisRecorder_<<ni<<",";
    }
  }else std::cout << "Unable to open "<<NIS_LIDAR_FILENAME<<std::endl;
  lidarNisRecorder_.close();

  std::cout<<"Size of radar vector: "<<radarNIS_.size()<<std::endl;
  if(!radarNisRecorder_.is_open()){
    radarNisRecorder_.open(NIS_RADAR_FILENAME);
  }
  if(radarNisRecorder_.is_open()){
    for(double ni : radarNIS_){
      radarNisRecorder_<<ni<<",";
    }
  }else std::cout << "Unable to open "<<NIS_RADAR_FILENAME<<std::endl;
  radarNisRecorder_.close();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    //clear nis file 
    

    // first measurement
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      //initiate state covariance with identity matrix
      P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
      
      //Convert position data from polar coordinates to cartesian coordinates
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      //initialize state mean vector
      x_ << rho*cos(phi), rho*sin(phi), 0, 0 , 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      //initiate state covariance matrix with smaller px and py by variances
      P_ << 0.5, 0, 0, 0, 0,
            0, 0.5, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
      //set the state with the initial location and zero velocity
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0 , 0;
    }

    time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    std::cout<<"finish Initialization"<<std::endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  //unscented preprocess
  //TODO

  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {

    //measurement update for radar
    UpdateRadar(meas_package);

  } else {
    if(use_laser_){
      //measurement update for lidar
      UpdateLidar(meas_package);
    }
  }
  // print the output
  // cout << "x_ = \n" << ekf_.x_ << endl;
  // cout << "P_ = \n" << ekf_.P_ << endl;
  // cout <<endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //=============================Start of creating augmented sigma points==================================//
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, AUG_SIG_NUM);
 
  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  //=============================End of creating augmented sigma points==================================//

  //=============================Start of predicting sigma points==================================//
  //predict sigma points
  for (int i = 0; i< AUG_SIG_NUM; i++)
  {
    //extract values for better readability
    double p_x      = Xsig_aug(0,i);
    double p_y      = Xsig_aug(1,i);
    double v        = Xsig_aug(2,i);
    double yaw      = Xsig_aug(3,i);
    double yawd     = Xsig_aug(4,i);
    double nu_a     = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.000001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //=============================End of predicting sigma points==================================//

  //=============================Start of predicting mean and covariance==================================//

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < AUG_SIG_NUM; i++) {  //iterate over sigma points
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < AUG_SIG_NUM; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
  //=============================End of predicting mean and covariance==================================//
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  /**=============================Start of predicting Radar Measurement==================================*/
  
  //set measurement dimension, lidar can measure Cartesian cordinates x, y
  int n_z = 2;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd( n_z, AUG_SIG_NUM);
  Zsig = Xsig_pred_.topRows(2);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < AUG_SIG_NUM; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  //  //angle normalization
  // while (z_pred(3)> M_PI) z_pred(3)-=2.*M_PI;
  // while (z_pred(3)<-M_PI) z_pred(3)+=2.*M_PI;
  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < AUG_SIG_NUM; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // //angle normalization
    // while (z_diff(3)> M_PI) z_diff(3)-=2.*M_PI;
    // while (z_diff(3)<-M_PI) z_diff(3)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;
  /**=============================End of predicting Radar Measurement==================================*/

  /**=============================Start of update==================================*/
  //create example vector for incoming lidar measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
   //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < AUG_SIG_NUM; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    // while (z_diff(3)> M_PI) z_diff(3)-=2.*M_PI;
    // while (z_diff(3)<-M_PI) z_diff(3)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // //angle normalization
    // while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    // while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  // while (z_diff(3)> M_PI) z_diff(3)-=2.*M_PI;
  // while (z_diff(3)<-M_PI) z_diff(3)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ = P_ - K*S*K.transpose();

  /**=============================End of update==================================*/

  /**=============================Start of calculating NIS==================================*/
  VectorXd z_diffs = z - z_pred;
  double nis = z_diffs.transpose() * S.inverse() * z_diffs;

  lidarNIS_.push_back(nis);
  if(lidarNIS_.size() >= 20){
      std::cout<<"Size of lidar vector: "<<lidarNIS_.size()<<"\n";
    if(!lidarNisRecorder_.is_open()){
      lidarNisRecorder_.open(NIS_LIDAR_FILENAME,std::ios::app);
    }
    if(lidarNisRecorder_.is_open()){
      for(double ni : lidarNIS_){
        lidarNisRecorder_<<ni<<"\n";
      }
    }else std::cout << "Unable to open "<<NIS_LIDAR_FILENAME<<std::endl;
    lidarNisRecorder_.close();
    lidarNIS_.clear();
  }
  //write in to file
  // if (!lidarNisRecorder_.is_open())
  // {
  //   lidarNisRecorder_.open(NIS_LIDAR_FILENAME);
  // }
  // if(lidarNisRecorder_.is_open()){
  //   lidarNisRecorder_ << nis;
  //   lidarNisRecorder_ << ",";
  // }
  // else std::cout << "Unable to open "<<NIS_LIDAR_FILENAME<<std::endl;
  /**=============================End of calculating NIS==================================*/
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  /**=============================Start of predicting Radar Measurement==================================*/
  
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, AUG_SIG_NUM);
  //transform sigma points into measurement space
  for (int i = 0; i < AUG_SIG_NUM; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < AUG_SIG_NUM; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < AUG_SIG_NUM; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;
  /**=============================End of predicting Radar Measurement==================================*/

  /**=============================Start of update==================================*/
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
   //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < AUG_SIG_NUM; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ = P_ - K*S*K.transpose();

  /**=============================End of update==================================*/

  /**=============================Start of calculating NIS==================================*/
  VectorXd z_diffs = z - z_pred;
  double nis = z_diffs.transpose() * S.inverse() * z_diffs;
  radarNIS_.push_back(nis);
  if(radarNIS_.size() >= 20){
      std::cout<<"Size of radar vector: "<<radarNIS_.size()<<"\n";
    if(!radarNisRecorder_.is_open()){
      radarNisRecorder_.open(NIS_RADAR_FILENAME,std::ios::app);
    }
    if(radarNisRecorder_.is_open()){
      for(double ni : radarNIS_){
        radarNisRecorder_<<ni<<"\n";
      }
    }else std::cout << "Unable to open "<<NIS_RADAR_FILENAME<<std::endl;
    radarNisRecorder_.close();
    radarNIS_.clear();
  }
  //write into file
  // if (!radarNisRecorder_.is_open())
  // {
  //   radarNisRecorder_.open(NIS_RADAR_FILENAME);
  // }
  // if(radarNisRecorder_.is_open()){
  //   radarNisRecorder_ << nis;
  //   radarNisRecorder_ << ",";
  // }
  // else std::cout << "Unable to open "<<NIS_RADAR_FILENAME<<std::endl;
}
