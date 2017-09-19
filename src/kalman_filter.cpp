#include "kalman_filter.h"
#include "tools.h"

// STE
#include <iostream>
using namespace std;

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

  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
  
}

void KalmanFilter::Update(const VectorXd &z) {
  
  VectorXd z_pred = H_ * x_;  
  VectorXd y = z - z_pred; 
  MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;
	
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  
  Tools tools;
    
  //recover state parameters
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	
	//compute h(x')	
	VectorXd z_pred = VectorXd(3);
	z_pred << c2, atan2(py,px), (px*vx+py*vy)/c2;   // TODO Controllare risultati di atan
	
	cout << "RADAR: passato da qui 1" << endl;
  
  VectorXd y = z - z_pred;
  MatrixXd Hj_ = tools.CalculateJacobian(x_);
  MatrixXd Hjt = Hj_.transpose();
	MatrixXd S = Hj_ * P_ * Hjt + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Hjt;
	MatrixXd K = PHt * Si;
	
	x_ = x_ + (K * y);
	long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj_) * P_;
  
}
