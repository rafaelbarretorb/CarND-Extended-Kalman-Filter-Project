#include <iostream>
#include "tools.h"
#include <string>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;




Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth)
{
    /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	VectorXd residual;
	if (estimations.size() != ground_truth.size() || estimations.size()==0)
	{
		cout << "ERROR RMSE" << endl;
		cout << estimations.size() << endl;
		
		return rmse;
	}

	for(unsigned int i=0; i < estimations.size(); ++i)
	{
		VectorXd residual;
		residual = estimations[i] - ground_truth[i];
		residual = residual.array()*residual.array();
		rmse += residual;

	}
    
	//calculate the mean
	rmse = rmse/(estimations.size());

	//calculate the squared root
	rmse = (rmse.array()).sqrt();

	//return the result
	cout << "RMSE = ["<< rmse(0) << " , "<< rmse(1) << " , "<< rmse(2)<< " , " << rmse(3) << " ]"<< endl;
	cout << estimations.size() << endl;
	return rmse;

    
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) 
{
  /**
  TODO:
    * Calculate a Jacobian here.
  */

    MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float square_p = px*px + py*py; 

	if (square_p > 0.00001)
	{
		Jacobian_OK = true;
		Hj(0,0) = px/sqrt(square_p);
		Hj(0,1) = py/sqrt(square_p);;
		Hj(0,2) = 0;
		Hj(0,3) = 0;

		Hj(1,0) = - py/square_p;
		Hj(1,1) = px/square_p;
		Hj(1,2) = 0;
		Hj(1,3) = 0;

		Hj(2,0) = (py*(vx*py - vy*px))/pow(square_p,1.5);
		Hj(2,1) = (px*(vy*px - vx*py))/pow(square_p,1.5);
		Hj(2,2) = Hj(0,0);
		Hj(2,3) = Hj(0,1);
	}
	else
	{
		cout << "ERROR Jacobian" << endl;
		Jacobian_OK = false;
		return Hj;
	}
	 
	return Hj;

}


VectorXd Tools::Conversion_Polar2Cartesian(const VectorXd &z)
{
	VectorXd z_radar_cartesian;

	float phi;
	phi = Tools::AngleNorm(phi);

	z_radar_cartesian(0) = cos(phi)*z(0);
	z_radar_cartesian(1) = sin(phi)*z(0);      
	

	return z_radar_cartesian;
}

float Tools::AngleNorm(const float &phi)
{
	float phi_norm = phi;
	PI = 3.14159265358979;

	while(phi_norm > PI)
        phi_norm -= 2*PI;

    while(phi_norm <-PI)
        phi_norm += 2*PI;

	return phi_norm;
}
