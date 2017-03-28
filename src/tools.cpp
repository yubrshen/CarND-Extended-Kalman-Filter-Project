#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size

  if ((0 < estimations.size()) && (estimations.size() == ground_truth.size())) {
    VectorXd r, s;
    for(unsigned int i=0; i < estimations.size(); ++i){
      r = (estimations[ i ] - ground_truth[ i ]);
      s = r.array()*r.array();
      rmse += s;
    }
    rmse = rmse/estimations.size();
    rmse = rmse.array().sqrt();
  } else {
    cout << "Invalid estimation or ground_truth data" << endl;
  }
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
}
