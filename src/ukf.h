#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* Augmented mean stat vector;
  VectorXd x_aug_;

  ///* Augmented state covariance matrix
  MatrixXd P_aug_;

  ///* X augmented sigma points
  MatrixXd X_aug_sig_;

  ///* predicted augmented sigma points matrix
  MatrixXd X_aug_sig_pred_;

  ///* Measurement augmented sigma points
  MatrixXd Z_aug_sig_;


  ///* The mesasurement from sensor
  VectorXd z_laser_;
  VectorXd z_radar_;

  ///* mean predicted measurement
  VectorXd z_laser_pred_;
  VectorXd z_radar_pred_;

  ///* measurement covariance matrix S
  MatrixXd S_radar_;

  ///* laser measurement covariance matrix
  MatrixXd R_laser_;

  ///* laser measurement matrix
  MatrixXd H_laser_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Keep track of timestamp
  long long previous_timestampe_;

  ///* Keep track of times difference
  long long dt_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Measurement dimension
  int n_z_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* Augmented Sigma point spreading parameter
  double lambda_aug_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  ///* process mode
  enum ProcessMode {
    BOTH,
    LASER_ONLY,
    RADAR_ONLY
  } process_mode_ ;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * Generate augmented sigma Points
  */
  void AugmentedSigmaPoints();

  /**
   * Sigma point prediction
   */
  void SigmaPointPrediction();

  /**
   * Predict mean and covariance
   */
  void PredictMeanAndCovariance();

  /**
   * Precit Radar Measurement
   */
  void PredictRadarMeasurement();

  /**
   * Update radar state
   */
  void UpdateRadarState();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);
};

#endif /* UKF_H */
