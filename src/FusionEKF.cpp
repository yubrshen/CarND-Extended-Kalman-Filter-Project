#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  H_laser_ <<
    1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0;

  // Hj will be computed when needed

  //measurement covariance matrix - laser
  R_laser_ <<
    0.0225, 0,
    0,      0.0225;

  //measurement covariance matrix - radar
  R_radar_ <<
    0.09, 0,      0,
    0,    0.0009, 0,
    0,    0,      0.09;

  /**
  DONE:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  // The following variables (objects)
  // local statically allocated in this constructor would live as long as the
  // object of the class lives.

  // initializing ekf_
  VectorXd x_in(4);
  MatrixXd P_in(4, 4);
	P_in <<
    1, 0, 0,    0,
    0, 1, 0,    0,
    0, 0, 1000, 0,
    0, 0, 0,    1000;

  MatrixXd F_in(4, 4); // To be updated in ProcessMeasurement by delta_t
  F_in <<              // identity as background
    1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0;

  MatrixXd Q_in(4, 4); // To be updated in ProcessMeasurement by delta_t
  Q_in <<              // all zeros as background
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0;

  MatrixXd H_in; // to be substituted according to measurement type
  MatrixXd R_in; // to be substituted according to measurement type

  ekf_.Init(x_in, P_in, F_in, H_in, R_in, Q_in);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

inline float SquaredDistance(const float& px, const float& py) {
  return px*px + py*py;
}

bool FusionEKF::GoodMeasurement(const MeasurementPackage &measurement_pack) {
  switch (measurement_pack.sensor_type_) {
  case MeasurementPackage::RADAR:
    if (0.0001 < measurement_pack.raw_measurements_[0]) {
      return true;
    }
    break;
  case MeasurementPackage::LASER:
    if (0.0001 < SquaredDistance(measurement_pack.raw_measurements_[0],
                                       measurement_pack.raw_measurements_[1])) {
      return true;
    }
    break;
  }
  return false;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  // pre-condition: the measurement is acceptable for processing

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
       DONE:
       * Initialize the state ekf_.x_ with the first measurement.
       * Create the covariance matrix. (the covariance P is initialized in the FusionEKF())
       * Remember: you'll need to convert radar from polar to cartesian coordinates.
       */
    // first measurement
    cout << "EKF: " << endl;
    previous_timestamp_ = measurement_pack.timestamp_;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
         Convert radar from polar to cartesian coordinates and initialize state.
      */
      ekf_.x_ <<
        measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]),
        measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]),
        measurement_pack.raw_measurements_[2] * cos(measurement_pack.raw_measurements_[1]),
        measurement_pack.raw_measurements_[2] * sin(measurement_pack.raw_measurements_[1]);
      // assume that the speed along rho can be also projected on x, and y axis.
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
         Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    // done initializing, no need to predict or update
    is_initialized_ = true;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   DONE:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix. (These constants are defined as FusionEKE members)
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  ekf_.Q_(0, 0) = dt_4/4*noise_ax; ekf_.Q_(0, 2) = dt_3/2*noise_ax;
  ekf_.Q_(1, 1) = dt_4/4*noise_ay; ekf_.Q_(1, 3) = dt_3/2*noise_ay;
  ekf_.Q_(2, 0) = dt_3/2*noise_ax; ekf_.Q_(2, 2) = dt_2*noise_ax;
  ekf_.Q_(3, 1) = dt_3/2*noise_ay; ekf_.Q_(3, 3) = dt_2*noise_ay;

  // ekf_.Q_ <<
  //   dt_4/4*noise_ax, 0,               dt_3/2*noise_ax, 0,
  //   0,               dt_4/4*noise_ay, 0,               dt_3/2*noise_ay,
  //   dt_3/2*noise_ax, 0,               dt_2*noise_ax,   0,
  //   0,               dt_3/2*noise_ay, 0,               dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   DONE:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
