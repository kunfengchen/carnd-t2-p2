#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;
  lambda_aug_ = 3 - n_aug_;

  // new values
  std_a_ = 0.2;
  std_yawdd_ = 0.2;


}

UKF::~UKF() {}

/**
 * REF:
 * https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/daf3dee8-7117-48e8-a27a-fc4769d2b954/concepts/a01da3d5-9c21-409a-b775-9b237987df46
 * @param Xsig_out
 */
void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {
  MatrixXd A = P.llt().matrixL();
  // *Xsig_out = MatrixXd(n_x_, 2 * n_x_ + 1);

  Xsig_out->col(0) = x_;

  for (int i =0; i < n_x_; i++) {
    Xsig_out->col(i+1) = x_  + sqrt(lambda_ +i n_x_) * A.col(i);
    Xsig_out->col(i+1+n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
  }
}

/**
 * REF:
 * https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/daf3dee8-7117-48e8-a27a-fc4769d2b954/concepts/07444c9b-b4be-4615-96e1-e1b221a9add6
 * @param Xsig_out
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  // *Xsig_out = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // augmented mean state
  Xsig_out->head(5) = x_;
  Xsig_out->(5) = 0;
  Xsig_out->(6) = 0;

  // augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;
  // square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // augmented sigma points
  Xsig_out->col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_out->col(i+1) = x_aug + sqrt(lambda_aug_+n_aug_) * A.col(i);
    Xsig_out->col(i+1+n_aug_) = x_aug - sqrt(lambda_aug_+n_aug_) * A.col(i);
  }
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
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
}
