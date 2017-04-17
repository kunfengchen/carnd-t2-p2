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
  // x_ = VectorXd(5);

  // initial covariance matrix
  // P_ = MatrixXd(5, 5);

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
  DONE:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = 5;
  n_aug_ = 7;
  // measurement dimension for radar: r, phi, and rdot
  n_z_ = 3;
  lambda_ = 3 - n_x_;
  lambda_aug_ = 3 - n_aug_;

  // new values
  /// std_a_ = 2.5; // 0.8 // from slack
  /// std_yawdd_ = 2.0; // 0.6 // from slack
  std_a_ = 0.2;
  std_yawdd_ = 0.2;
  std_radr_ = 0.3;
  std_radphi_ = 0.0175;
  std_radrd_ = 0.1;

  x_ = VectorXd(n_x_);
  P_ = MatrixXd(n_x_, n_x_);
  x_aug_ = VectorXd(n_aug_);
  P_aug_ = MatrixXd(n_aug_, n_aug_);
  X_aug_sig_ = MatrixXd(n_aug_, 2*n_aug_+1);
  X_aug_sig_pred_ = MatrixXd(n_x_, 2*n_aug_+1); // dimensions!!
  Z_aug_sig_ = MatrixXd(n_z_, 2*n_aug_+1);
  dt_ = 0;

  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_aug_/(lambda_aug_+n_aug_);
  double weight_others = 0.5/(lambda_aug_+n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2*n_aug_+1; i++) { // 2n+1 weights
    weights_(i) = weight_others;
  }

  // mean predicted measurement
  z_pred_ = VectorXd(n_z_);
  //measurement covariance matrix S
  S_ = MatrixXd(n_z_,n_z_);
}

UKF::~UKF() {}

/**
 * REF:
 * https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/daf3dee8-7117-48e8-a27a-fc4769d2b954/concepts/07444c9b-b4be-4615-96e1-e1b221a9add6
 */
void UKF::AugmentedSigmaPoints() {

  // augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  // augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(5, 5) = std_a_ * std_a_;
  P_aug_(6, 6) = std_yawdd_ * std_yawdd_;
  // square root matrix
  MatrixXd A = P_aug_.llt().matrixL();

  // augmented sigma points
  X_aug_sig_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; i++) {
    X_aug_sig_.col(i+1) = x_aug_ + sqrt(lambda_aug_+n_aug_) * A.col(i);
    X_aug_sig_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_aug_+n_aug_) * A.col(i);
  }
}

/**
 * REF:
 * https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/daf3dee8-7117-48e8-a27a-fc4769d2b954/concepts/19ab81df-0719-4e9c-8c04-236f7449261a
 */
void UKF::SigmaPointPrediction() {
  // predict
  double px, py, v, yaw, yawd, nu_a, nu_yawdd;
  // predicted stat values
  double px_pred, py_pred, v_pred, yaw_pred, yawd_pred;

  for (int i = 0; i<2*n_aug_+1; i++) {
    px = X_aug_sig_(0, i);
    py = X_aug_sig_(1, i);
    v = X_aug_sig_(2, i);
    yaw = X_aug_sig_(3, i);
    yawd = X_aug_sig_(4, i);
    nu_a = X_aug_sig_(5, i);
    nu_yawdd = X_aug_sig_(6, i);

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_pred = px + v/yawd * (sin(yaw+yawd*dt_)-sin(yaw));
      py_pred = py + v/yawd * (cos(yaw)-cos(yaw+yawd*dt_));
    } else {
      px_pred = px + v*dt_*cos(yaw);
      py_pred = py + v*dt_*sin(yaw);
    }

    v_pred = v;
    yaw_pred = yaw + yawd*dt_;
    yawd_pred = yawd;

    // add noise
    px_pred = px_pred + 0.5*nu_a*dt_*dt_*cos(yaw);
    py_pred = py_pred + 0.5*nu_a*dt_*dt_*sin(yaw);
    v_pred = v_pred + nu_a*dt_;

    yaw_pred = yaw_pred + 0.5*nu_yawdd*dt_*dt_;
    yawd_pred = yawd_pred+ nu_yawdd*dt_;

    // assign back to the right column
    X_aug_sig_pred_(0, i) = px_pred;
    X_aug_sig_pred_(1, i) = py_pred;
    X_aug_sig_pred_(2, i) = v_pred;
    X_aug_sig_pred_(3, i) = yaw_pred;
    X_aug_sig_pred_(4, i) = yawd_pred;
  }

}

/**
 * REF:
 * https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/daf3dee8-7117-48e8-a27a-fc4769d2b954/concepts/07b59fdd-adb3-479b-8566-336332cf5f09
 */
void UKF::PredictMeanAndCovariance() {

  // calculate
  x_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    x_ = x_+weights_(i)*X_aug_sig_pred_.col(i);
  }

  P_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    // state difference
    VectorXd x_diff = X_aug_sig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) >  M_PI) { x_diff(3)-=2.*M_PI; }
    while (x_diff(3) < -M_PI) { x_diff(3)+=2.*M_PI; }

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * REF:
 * https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/daf3dee8-7117-48e8-a27a-fc4769d2b954/concepts/d1a4ce03-73aa-4653-b216-1c6d965fc216
 */
void UKF::PredictRadarMeasurement() {
  double px, py, v, yaw, v1, v2;

  // transform sigma points into measuremnet space
  for (int i = 0; i < 2*n_aug_+1; i++) {
    px = X_aug_sig_pred_(0, i);
    py = X_aug_sig_pred_(1, i);
    v  = X_aug_sig_pred_(2, i);
    yaw = X_aug_sig_pred_(3, i);
    v1 = cos(yaw)*v;
    v2 = sin(yaw)*v;
    Z_aug_sig_(0, i) = sqrt(px*px+py*py);  // r
    Z_aug_sig_(1, i) = atan2(py, px);      // phi
    Z_aug_sig_(2, i) = (px*v1+py*v2) / sqrt(px*px+py*py);  //r_dot
  }

  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred_ = z_pred_ +weights_(i)*Z_aug_sig_.col(i);
  }

  S_.fill(0.0);
  for (int i = 0; i<2*n_aug_+1; i++) {
    // residual
    VectorXd z_diff = Z_aug_sig_.col(i) - z_pred_;

    // angle normalization
    while (z_diff(1) > M_PI) { z_diff(1)-=2.*M_PI; }
    while (z_diff(1) < -M_PI) { z_diff(1)+=2.*M_PI; }
    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z_, n_z_);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;

  S_ = S_ + R;
}

/**
 * REF:
 * https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/daf3dee8-7117-48e8-a27a-fc4769d2b954/concepts/e19bc36c-a671-4799-8c63-cec40544c2aa
 */
void UKF::UpdateState() {

  // TODO z
  // matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    // residaul
    VectorXd z_diff = Z_aug_sig_.col(i) - z_pred_;
    // angle normalization
    while (z_diff(1) > M_PI) { z_diff(1) -= 2.*M_PI; }
    while (z_diff(1) < -M_PI) { z_diff(1) += 2.*M_PI; }

    // state difference
    VectorXd x_diff = X_aug_sig_pred_.col(i) - x_;

    // angle normalization
    while (x_diff(3) > M_PI) { x_diff(3) -= 2.*M_PI; }
    while (x_diff(3) < -M_PI) { x_diff(3) += 2.*M_PI; }

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K
  MatrixXd K = Tc * S_.inverse();

  // residual
  VectorXd z_diff = z_ - z_pred_;

  // angle normalization
  while (z_diff(1) > M_PI) { z_diff(1) -= 2.*M_PI; }
  while (z_diff(1) < -M_PI) { z_diff(1) += 2.*M_PI; }

  // update state mean and covariance matrix

  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_*K.transpose();
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

  ///// Initialization
  if (!is_initialized_) {
    previous_timestampe_ = meas_package.timestamp_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       * Convert radar from polar to cartesian coordiates and initalize state.
       */
       double ro = meas_package.raw_measurements_[0];
       double theta = meas_package.raw_measurements_[1];
       double x = ro * cos(theta);
       double y = ro * sin(theta);
       x_ << x, y, 0.0, 0.0, 0.0, 0.0, 0.0; //TODO
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
       x_ << meas_package.raw_measurements_[0],
             meas_package.raw_measurements_[1],
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    }

    if (x_.norm() > 1e-4) {
      // done init. no need to predict or update
      is_initialized_ = true;
    } else {
      cout << "input is too small to initialied UKF." << endl;
    }
  }

  ///// Prediction
  // compute the time difference
  dt_ = (meas_package.timestamp_ - previous_timestampe_) / 1000000.0; // dt in seconds

  z_ = meas_package.raw_measurements_;

  Prediction(dt_);

  ///// Update
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar update
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    // Laser update
    UpdateLidar(meas_package);
  }

  UpdateState();

  cout << "x_ = " << endl << x_ << endl;
  cout << "P_ = " << endl << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  DONE:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  AugmentedSigmaPoints();
  SigmaPointPrediction();
  PredictMeanAndCovariance();
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
  DONE:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  PredictRadarMeasurement();
}
