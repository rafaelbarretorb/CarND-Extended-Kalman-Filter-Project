#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() 
{
	is_initialized_ = false;

	previous_timestamp_ = 0;

	// initializing matrices
	R_laser_ = MatrixXd(2, 2);	// measurement covariance matrix - laser
	R_radar_ = MatrixXd(3, 3);	// measurement covariance matrix - radar
	H_laser_ = MatrixXd(2, 4);
	H_laser_ << 1, 0, 0, 0,
		     	0, 1, 0, 0;

	Hj_      = MatrixXd(3, 4);

	//
	R_laser_ << 0.0225, 0,
		        0, 0.0225;

	
	R_radar_ << 0.09, 0, 0,
		        0, 0.0009, 0,
		        0, 0, 0.09;

	/**
	TODO:
	* Finish initializing the FusionEKF.
	* Set the process and measurement noises
	*/
	// 

	// State transition matrix
	F_ = MatrixXd(4, 4);
	F_ << 1, 0, 1, 0,
		  0, 1, 0, 1,
		  0, 0, 1, 0,
		  0, 0, 0, 1;	

    // P and Q are set from Lesson "Kalman Filter Equations in C++ Part 2"
	P_ = MatrixXd(4, 4);
	P_ << 1000, 0, 0, 0,
		  0, 1000, 0, 0,
		  0, 0, 1000, 0,
		  0, 0, 0, 1000;

	// Process covariance matrix
	Q_ = MatrixXd(4, 4);
	Q_ << 0, 0, 0, 0,
		  0, 0, 0, 0,
		  0, 0, 0, 0,
		  0, 0, 0, 0;

	// Measurements Vectors
	z_laser_ = VectorXd(2);
	z_radar_ = VectorXd(3);
	z_radar_cartesian = VectorXd(2);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) 
{
  /*****************************************************************************
   *  Initialization

	meas_px      
    meas_py                          
    timestamp
   ****************************************************************************/

	// LIDAR

	if (!is_initialized_) // if initialized is False
	{
		/**
		TODO:
		  * Initialize the state ekf_.x_ with the first measurement.
			- state vector x = (px py vx vy)


		  * Create the covariance matrix.
		  * Remember: you'll need to convert radar from polar to cartesian coordinates.
		*/
		// first measurement
		cout << "EKF: " << endl;
		ekf_.x_ = VectorXd(4);

		// velocity is zero at the beginning
		ekf_.x_ << 1, 1, 0, 0;

		current_timestamp_ = measurement_pack.timestamp_;

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
		{
			/**
		  	Convert radar from polar to cartesian coordinates and initialize state.
			*/

		  	// measure z RADAR
			z_radar_(0) = measurement_pack.raw_measurements_(0); // meas_rho
			z_radar_(1) = measurement_pack.raw_measurements_(1); // meas_phi
			z_radar_(2) = measurement_pack.raw_measurements_(2); // meas_rho_dot

			// Convert radar from polar to cartesian coordinates  
			z_radar_cartesian = tools.Conversion_Polar2Cartesian(z_radar_);         

			// Initialize state
			ekf_.x_(0) = z_radar_cartesian(0); // px
			ekf_.x_(1) = z_radar_cartesian(1); // py 

			// 
			ekf_.Init(ekf_.x_, P_ , F_, Hj_, R_radar_, Q_);
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
		{
		  /**
		  Initialize state.
		  */

			z_laser_(0) = measurement_pack.raw_measurements_(0); // measure py 
			z_laser_(1) = measurement_pack.raw_measurements_(1); // 

			ekf_.x_(0) = z_laser_(0); // px
			ekf_.x_(1) = z_laser_(1); // py
	
			ekf_.Init(ekf_.x_, P_ , F_, H_laser_, R_laser_, Q_);
		}

		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
	 * Update the state transition matrix F according to the new elapsed time.
	  - Time is measured in seconds.
	 * Update the process noise covariance matrix.
	 * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
	 
	previous_timestamp_ = current_timestamp_;
	current_timestamp_ = measurement_pack.timestamp_;

	//compute the time elapsed between the current and previous measurements
	dt = (current_timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds


	// 1. Update the state transition matrix F according to the new elapsed time
	ekf_.F_(0,2) = dt;
	ekf_.F_(1,3) = dt; 
	
	// 2. Update the process noise covariance matrix
	noise_ax = 9;
	noise_ay = 9;

	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ = Q_;

	ekf_.Q_(0,0) = pow(dt,4)*noise_ax/4;
	ekf_.Q_(0,2) = pow(dt,3)*noise_ax/2;

	ekf_.Q_(1,1) = pow(dt,4)*noise_ay/4;
	ekf_.Q_(1,3) = pow(dt,3)*noise_ay/2;

	ekf_.Q_(2,0) = pow(dt,3)*noise_ax/2;
	ekf_.Q_(2,2) = pow(dt,2)*noise_ax;

	ekf_.Q_(3,1) = pow(dt,3)*noise_ay/2;
	ekf_.Q_(3,3) = pow(dt,2)*noise_ay;



	// 3. Predict
	ekf_.Predict();

 	/*****************************************************************************
	*  Prediction
    ****************************************************************************/

	// 4. Update 
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
	{
  		// RADAR measure z 
		z_radar_(0) = measurement_pack.raw_measurements_(0); // meas_rho
		z_radar_(1) = measurement_pack.raw_measurements_(1); // meas_phi
		z_radar_(2) = measurement_pack.raw_measurements_(2); // meas_rho_dot

		ekf_.R_ = R_radar_;

		Hj_ = tools.CalculateJacobian(ekf_.x_);
		
		// Radar updates
		ekf_.H_ = Hj_;

		// Test the calculated Jacobian Matrix
		if (tools.Jacobian_OK == true)
			ekf_.UpdateEKF(z_radar_);
	} 
 	else
    {
    	// LASER measure z
    	z_laser_(0) = measurement_pack.raw_measurements_(0); // measure px
		z_laser_(1) = measurement_pack.raw_measurements_(1); // measure py


    	ekf_.R_ = R_laser_;
		ekf_.H_ = H_laser_;

	    // Laser update
		ekf_.Update(z_laser_);
    }

	// print the output
	cout << "x_ = " << ekf_.x_ << endl;
	cout << "P_ = " << ekf_.P_ << endl;
}
