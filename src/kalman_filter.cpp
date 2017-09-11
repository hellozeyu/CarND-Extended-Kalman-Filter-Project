#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;


KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

const double PI = 3.1415926;
void normalize_phi(double & phi)
{
    while(phi<-PI)
    {
        phi += 2*PI;
    }

    while(phi>PI)
    {
        phi -= 2*PI;
    }
}


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
    cout << x_ << '\n'<< endl;
    cout << "========" << endl;
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

    // KF Measurement update step
    VectorXd y = z - H_ * x_;
    normalize_phi(y(1));
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;

    MatrixXd I; // Identity matrix
    //new state
    x_ = x_ + (K * y);
    long x_size = x_.size();
    I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    double rho = sqrt(x_(0)*x_(0)+x_(1)*x_(1));
    double phi = atan2(x_(1), x_(0));
    double rho_dot = (x_(2)*x_(0) + x_(3)*x_(1))/rho;

    VectorXd z_pred(3);
    z_pred<<rho, phi, rho_dot;

    VectorXd y = z - z_pred;
    normalize_phi(y(1));
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
