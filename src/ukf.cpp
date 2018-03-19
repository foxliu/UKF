#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 1.5;

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
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // not initialized untilf first process measurement
  is_initialized_ = false;

  // set state dimentsion
  n_x_ = 5;

  // set augmented dimension
  n_aug_ = 7;

  // set lambda
  lambda_ = 3 - n_aug_;

  // sqrt lambda and aug
  sqrt_l_a = sqrt(lambda_ + n_aug_);

  // set the set the sigma points cols
  sigma_length = 2 * n_aug_ + 1;

  // init the time
  time_us_ = 0;

  // set the radar and ladar dimension
  n_z_radar_ = 3;
  n_z_ladar_ = 2;

  Xsig_pred_ = MatrixXd(n_x_, sigma_length);

  // the weights
  weights_ = VectorXd(sigma_length);
  double weight_start = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_start;

  double weight_i = 0.5 / (lambda_ + n_aug_);
  for (int i = 1; i < sigma_length; ++i)
  {
    weights_(i) = weight_i;
  }

  // the R for ladar and radar
  R_ladar_ = MatrixXd(n_z_ladar_, n_z_ladar_);
  R_ladar_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;
  
  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
              0, std_radphi_ * std_radphi_, 0,
              0, 0, std_radrd_ * std_radrd_;

}

UKF::~UKF() {}

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
  if(!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
      double ro, phi, ro_dot;
      ro = meas_package.raw_measurements_[0];
      phi = meas_package.raw_measurements_[1];
      ro_dot = meas_package.raw_measurements_[2];
      // x_ << ro * cos(phi), ro * sin(phi), 0, 0, 0;
      // x_ << ro * cos(phi), ro * sin(phi), 0, ro_dot * cos(phi), ro_dot * sin(phi);
      x_ << ro * cos(phi), ro * sin(phi), 4, ro_dot * cos(phi), ro_dot * sin(phi);

      ///* state covariance matrix

      // P_ << 1, 0, 0, 0, 0,
      //       0, 1, 0, 0, 0,
      //       0, 0, 1, 0, 0,
      //       0, 0, 0, 1, 0,
      //       0, 0, 0, 0, 1;

      P_ << std_radr_ * std_radr_, 0, 0, 0, 0,
          0, std_radr_ * std_radr_, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, std_radphi_ * std_radphi_, 0,
          0, 0, 0, 0, std_radphi_ * std_radphi_;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      ///* state covariance matrix
      // P_ << 1, 0, 0, 0, 0,
      //       0, 1, 0, 0, 0,
      //       0, 0, 1, 0, 0,
      //       0, 0, 0, 1, 0,
      //       0, 0, 0, 0, 1;

      P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
          0, std_laspy_ * std_laspy_, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_  = true;
    return;
  }
  double dt = double(meas_package.timestamp_ - time_us_) / 1000000;
  time_us_ = meas_package.timestamp_;

  if (dt > LEAST_DT) {
    Prediction(dt);
  }

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }
}

void UKF::AugmentedSigmaPoints(MatrixXd &Xsig_aug_) {
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, sigma_length);

  x_aug.head(5) = x_;
  x_aug[5] = 0;
  x_aug[6] = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_; 
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  MatrixXd A = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; ++i) {
    Xsig_aug.col(i+1) = x_aug + sqrt_l_a * A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt_l_a * A.col(i);
  }
  Xsig_aug_ = Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd Xsig_aug, double delta_t) {
  MatrixXd Xsig_pred = MatrixXd(n_x_, sigma_length);
  for (int i = 0; i < sigma_length; ++i) {
    double p_x, p_y, v, yaw, yawd, nu_a, nu_yawdd;
    p_x = Xsig_aug(0, i);
    p_y = Xsig_aug(1, i);
    v  = Xsig_aug(2, i);
    yaw = Xsig_aug(3, i);
    yawd = Xsig_aug(4, i);
    nu_a = Xsig_aug(5, i);
    nu_yawdd = Xsig_aug(6, i);

    double px_p, py_p;
    double yaw_p = yaw + yawd * delta_t;
    double sin_yaw = sin(yaw);
    double cos_yaw = cos(yaw);
    if (fabs(yawd) > LEAST_YAWD) {
      double r = v /yawd;
      px_p = p_x + r * (sin(yaw_p) - sin_yaw);
      py_p = p_y + r * (cos_yaw - cos(yaw_p));
    } else {
      px_p = p_x + v * delta_t * cos_yaw;
      py_p = p_y + v * delta_t * sin_yaw;
    }
    double v_p = v;
    double yawd_p = yawd;

    // add noise;
    double noise_delta = 0.5 * nu_a * delta_t * delta_t;
    px_p = px_p + noise_delta * cos_yaw;
    py_p = py_p + noise_delta * sin_yaw;
    v_p = v_p + nu_a * delta_t;
    yaw_p = yaw_p + noise_delta * nu_yawdd;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    Xsig_pred(0, i) = px_p;
    Xsig_pred(1, i) = py_p;
    Xsig_pred(2, i) = v_p;
    Xsig_pred(3, i) = yaw_p;
    Xsig_pred(4, i) = yawd_p;
  }
  Xsig_pred_ = Xsig_pred;
}

void UKF::PredictMeanAndCovariance()
{
  x_.fill(0.0);
  for (int i = 0; i < sigma_length; ++i)
  {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  P_.fill(0.0);
  for (int i = 0; i < sigma_length; ++i)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    NormalAngle(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  MatrixXd Xsig_aug_;
  this->AugmentedSigmaPoints(Xsig_aug_);
  this->SigmaPointPrediction(Xsig_aug_, delta_t);
  this->PredictMeanAndCovariance();
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
  MatrixXd Zsig;
  VectorXd z_pred;
  MatrixXd S;

  this->PredictMeasurement(z_pred, S, Zsig, n_z_ladar_);

  VectorXd z(n_z_ladar_);
  z << meas_package.raw_measurements_[0],
      meas_package.raw_measurements_[1];

  this->UpdateState(z, S, Zsig, z_pred, n_z_ladar_);
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
  MatrixXd Zsig;
  VectorXd z_pred;
  MatrixXd S;

  this->PredictMeasurement(z_pred, S, Zsig, n_z_radar_);

  VectorXd z(n_z_radar_);
  z << meas_package.raw_measurements_[0],
      meas_package.raw_measurements_[1],
      meas_package.raw_measurements_[2];
  
  this->UpdateState(z, S, Zsig, z_pred, n_z_radar_);
}

void UKF::NormalAngle(double &angle) {
  while (angle > M_PI) angle -= PI2;
  while (angle < -M_PI) angle += PI2;
}

void UKF::PredictMeasurement(VectorXd &z_out, MatrixXd &S_out, MatrixXd &Zsig_out, int n_z) {
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, sigma_length);
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  for (int i = 0; i < sigma_length; ++i) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    if (n_z == 3) {
      double v = Xsig_pred_(2, i);
      double phi = Xsig_pred_(3, i);
      double ro = sqrt(px * px + py * py);
      ro = RO_VALUE(ro);

      Zsig(0, i) = ro;
      Zsig(1, i) = ((px == 0 && py == 0) ? 0 : atan2(py, px));
      Zsig(2, i) = (px * cos(phi) * v + py * sin(phi) * v) / ro;
    } else if (n_z == 2) {
      Zsig(0, i) = px;
      Zsig(1, i) = py;
    }
  }
  z_pred.fill(0.0);
  for (int i = 0; i < sigma_length; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  S.fill(0.0);
  for (int i = 0; i < sigma_length; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    this->NormalAngle(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  if (n_z == 3) {
    S = S + R_radar_;
  } else if (n_z == 2) {
    S = S + R_ladar_;
  } else {
    throw length_error("Wrong z lengthe number");
  }
  S_out = S;
  z_out = z_pred;
  Zsig_out = Zsig;
}

void UKF::UpdateState(VectorXd z, MatrixXd S, MatrixXd Zsig, VectorXd z_pred, int n_z)
{
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Tc.fill(0.0);
  for (int i = 0; i < sigma_length; ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalAngle(x_diff(3));

    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (n_z == 3) {
      NormalAngle(z_diff(1));
    }
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  MatrixXd Si = S.inverse();
  MatrixXd K = Tc * Si;
  MatrixXd y = z - z_pred;
  if (n_z == 3) {
    NormalAngle(y(1));
  }
  x_ = x_ + K * y;
  P_ = P_ - K * S * K.transpose();
  
  if (n_z == 3) {
    NIS_Radar = y.transpose() * Si * y;
  } else if (n_z == 2) {
    NIS_Ladar = y.transpose() * Si * y;
  }
}