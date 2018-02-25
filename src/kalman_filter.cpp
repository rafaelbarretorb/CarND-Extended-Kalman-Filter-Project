// Include tools
#include "tools.h"

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

void KalmanFilter::Predict() 
{
    /**
    TODO:
    * predict the state
    */
    // Predict state
    x_ = F_ * x_ ;

    // Predict Object Covariance Matrix
    P_ = F_ * P_ * (F_.transpose()) + Q_;
}

void KalmanFilter::Update(const VectorXd &z) 
{
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

    z_predicted = H_ * x_;

    // Residual/Error
    VectorXd y = z - z_predicted;

    MatrixXd Ht = H_.transpose();
    MatrixXd S_ = H_ * P_ * Ht + R_;

    MatrixXd Si = S_.inverse();
    MatrixXd K_ =  P_ * Ht * Si;

    // Update the state
    x_ = x_ + K_*y;

    MatrixXd I_ =  MatrixXd::Identity(x_.size(), x_.size());

    P_ = (I_ - K_*H_)*P_;
}

void KalmanFilter::UpdateEKF(VectorXd &z)
{
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

    // Convert predicted state to polar
    // h(x)
    z_predicted = VectorXd(3);
    z_predicted(0) = sqrt((pow(x_(0), 2.0) + pow(x_(1), 2.0)));
    z_predicted(1) = atan2(x_(1),x_(0));
    z_predicted(2) = (x_(0)*x_(2) + x_(1)*x_(3))/z_predicted(0);

    //  y = z - h(x)
    VectorXd y = z - z_predicted;

    // Normalized Angle phi ---> -pi < phi < pi
    y(1) = tools.AngleNorm(y(1));

    MatrixXd Ht = H_.transpose();

    MatrixXd S_ = H_ * P_ * Ht + R_;
    MatrixXd Si = S_.inverse();
    MatrixXd K_ =  P_ * Ht * Si;

    x_ = x_ + K_*y;

    MatrixXd I_ =  MatrixXd::Identity(x_.size(), x_.size());

    P_ = (I_ - K_*H_)*P_;
}
