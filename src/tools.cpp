#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

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
  MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if((fabs(c1) < 0.0001) || (fabs(c3) < 0.0001)){
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

  float c4 = (px*vy -py*vx)/c3;
  float c5 = px/c2;
  float c6 = py/c2;

	//compute the Jacobian matrix
	Hj <<
    c5,       c6,      0,  0,
    -(py/c1), (px/c1), 0,  0,
    py*(-c4), px*c4,   c5, c6;

	return Hj;
}
