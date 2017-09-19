#include "kalman_filter.h"
#include <iostream>

#define PI 3.14159265358979323846
#define VERYSMALL 0.0001

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
  Update_(y);
  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
      
  //recover state parameters
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);
	
	//compute h(x')	
	float rho = sqrt(px*px+py*py);
	float theta = atan2(py,px);
	float rho_dot = (fabs(rho) < VERYSMALL) ? 0 : (px*vx+py*vy)/rho;
	VectorXd z_pred = VectorXd(3);
	z_pred << rho, theta, rho_dot;
	VectorXd y = z - z_pred;
	
	// Correcting dtheta to avoid discontinuities
	if (y(1) > PI) y(1) -= 2 * PI;
  if (y(1) < -PI) y(1) += 2 * PI;
	Update_(y);
	
}

// Operations common to Update and UpdateEKF
void KalmanFilter::Update_(const VectorXd &y) {

  MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht * Si;
    
  x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
  
}

