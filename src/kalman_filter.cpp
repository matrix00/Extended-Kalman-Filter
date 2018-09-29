#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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

	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {

  cout << "Laser update start "<< endl;
	VectorXd z_pred = H_ * x_;
  cout << "Laser update H z_pred "<< H_ << z_pred<< endl;
	VectorXd y = z - z_pred;
  cout << "Laser update y "<< y << endl;
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

 // print the output
  cout << "Laser update KF x_ = " << x_ << endl;
  cout << "Laser Update KF P_ = " << P_ << endl;


}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */

  cout << "Radar update start "<< endl;
	
  	Tools tools;

	VectorXd h_x = tools.CalculatehxPrime(x_) ;

  cout << "Radar update hx_ = " << h_x << endl;
	MatrixXd Hj = tools.CalculateJacobian(x_);

  cout << "Radar update Hj = " << Hj << endl;
	VectorXd y = z - h_x ;
  cout << "Radar update y = " << y  << endl;
	MatrixXd Ht = Hj.transpose();
	MatrixXd S = Hj * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;
  cout << "Radar update K = " << K << endl;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj) * P_;
 
// print the output
  cout << "Radar update KF x_ = " << x_ << endl;
  cout << "Radar Update KF P_ = " << P_ << endl;
}
