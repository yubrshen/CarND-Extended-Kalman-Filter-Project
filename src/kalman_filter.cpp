#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  DONE:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::UpdateMeta(const VectorXd &z, const VectorXd &y) {
  /**
  DONE:
    * update the state by using Kalman Filter equations
  */
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  UpdateMeta(z, y);
}

VectorXd state_comparable_to_measurement(const VectorXd& state, const MatrixXd& H) {
  // convert state to the form and value comparable to external measurement in the case of RADAR measurement

  // The following code computes y by the accurate but nonlinear equation:
  VectorXd z_pred(3);
  float px = state[0];
  float py = state[1];
  float vx = state[2];
  float vy = state[3];

  float c1 = sqrt(px*px + py*py);
  float c2 = px*vx + py*vy;

  float c3 = 0; // the default used if c1 is too small. It should be OK,
  // as when c1 is small, the c3, the rate of change of rho (rho_dot)
  // can be assumed to be very small.
  if (0.0001 < fabs(c1)) {
    c3 = c2/c1;
  }
  z_pred <<
    c1, atan2(py, px), c3;
  return z_pred;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  DONE:
    * update the state by using Extended Kalman Filter equations
  */
  VectorXd z_pred = state_comparable_to_measurement(x_, H_);
  VectorXd y = z - z_pred;

  // Begin adjusting y[1]
  // make sure y[1], the angle is within [-pi, pi]
  float full_circle = 2*M_PI;
  while (y[1] < -M_PI) {
    y[1] += full_circle;
  }
  while (M_PI < y[1]) {
    y[1] -= full_circle;
  }
  // End adjusting y[1]

  UpdateMeta(z, y);
}
