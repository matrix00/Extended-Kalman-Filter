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
	// ... your code here
	
	if (estimations.size() == 0 || (estimations.size() != ground_truth.size()))
	    return rmse;

    
	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        	VectorXd error = (estimations[i] - ground_truth[i]);
        	VectorXd errorsq  = error.array()*error.array();
        	rmse = rmse + errorsq;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = (rmse.array()).sqrt(); 

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	
	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//check division by zero
	if (px ==0 || py == 0)
	    Hj << 0, 0, 0,
	          0, 0, 0,
	          0, 0, 0;
	    
	
	//compute the Jacobian matrix
	float pxy_sqrt = sqrt(px*px + py*py);
	float pxy = pxy_sqrt * pxy_sqrt;
	float pxy_32 = pxy * pxy_sqrt;
	
	Hj(0,0) = px/pxy_sqrt;
	Hj(0,1) = py/pxy_sqrt;

	Hj(1,0) = -py/pxy;
	Hj(1,1) = px/pxy;
    
	Hj(2,0) = py*(vx*py-vy*px)/pxy_32;
	Hj(2,1) = px*(vy*px-vx*py)/pxy_32;
    
    
	Hj(2,2) = px/pxy_sqrt;
	Hj(2,3) = py/pxy_sqrt;
    
    
	return Hj;
}

