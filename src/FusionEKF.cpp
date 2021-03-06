#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

#include "debug.h"

#define SMALLNUM 0.0001

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/***
 * Constructor.
 *
 * Initialize key matrices
 *
 *   -R_laser:  convariance for Laser measurement update
 *   -R_radar:  covariance for Radar measurement update
 *   -H_laser:  maps belief from laser into sensor measurement space
 *   -Hj:	Jacobian for mapping from Radar (polar) space to sensor measurement space
 *
 * Initialize the Kalman filter: Note the Kalman filter implements the math for
 * the prediction and update.   
 */
FusionEKF::FusionEKF() {
  if (DEBUG)
    cout << "FusionEKF::FusionEKF" << endl;

  is_initialized_ = false;

  previous_timestamp_ = 0.0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  if (DEBUG)
    cout << "FusionEKF::FusionEKF  step" << endl;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
 
  ekf_.x_ = VectorXd(4);
  ekf_.x_ << 0,0,0,0;
 
  //initialize uncertainty covariance, P
  ekf_.P_= MatrixXd(4,4);
  ekf_.P_ << 1, 0, 0, 0,
           0, 1, 0, 0,
           0, 0, 1000, 0,
	   0, 0, 0, 1000;

  if (DEBUG)
    cout << "FusionEKF::FusionEKF  step" << endl;

  //initiallize state transition matrix, F
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, 0, 0,
            0, 1, 0, 0,
	    0, 0, 1, 0,
	    0, 0, 0, 1;

  //initialize measurement function, H
  ekf_.H_ = MatrixXd(4,4);
  ekf_.H_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;

  if (DEBUG)
    cout << "FusionEKF::FusionEKF  step" << endl;

  //initialize process covariance, Q
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << 0, 0, 0, 0,
            0, 0, 0, 0,
	    0, 0, 0, 0,
	    0, 0, 0, 0;


  if (DEBUG)
    cout << "FusionEKF::FusionEKF  step" << endl;

  //create the identity matrix (used in equations)
  ekf_.I_ = MatrixXd(4,4);
  ekf_.I_ << 1, 0, 0, 0,
             0, 1, 0, 0,
	     0, 0, 1, 0,
	     0, 0, 0, 1;
             
  //H_laser projects belief about object's current state
  //into the measurement space of sensor.  For LIDAR it
  //discards velocity information.  Recall state vector
  // x = (px,py,vx,vy) so we extract only position info.

  H_laser_ << 1, 0, 0, 0,
             0, 1, 0, 0;


  //Hj is the jacobian that procects radar measurement
  //into the position space
  Hj_ << 0,0,0,0,
         0,0,0,0,
 	0,0,0,0;


  //set the process noise
  noise_ax = 9.0;
  noise_ay = 9.0;

  //set the measurement noise, R above
   
  if (DEBUG)
    cout << "FusionEKF::FusionEKF  DONE" << endl;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}


/***
 * FusionEKF::ProcessMeasurement
 *
 * const MeasurementPackage& measurement_pack-  measurement instance consisting
 *                                              of timestamp, type, and data
 */
void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  float rho= 0.0;
  float rhodot = 0.0;
  float phi = 0.0;

  float px = 0.0;
  float py = 0.0;
  float vx = 0.0;
  float vy = 0.0;

  if (DEBUG) {
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "type= RADAR "  << "rho= " << measurement_pack.raw_measurements_[0] << "  ";
      cout << "phi= " << measurement_pack.raw_measurements_[1] << " ";
      cout << "rhodot= " << measurement_pack.raw_measurements_[2] << endl;
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      cout << "type= LASER ";
      cout << "px = " << measurement_pack.raw_measurements_[0] << " ";
      cout << "py = " << measurement_pack.raw_measurements_[1] << endl;
    }
  }

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/

  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    //ekf_.x_ << 0.6, 0.6, 5.199937, 0;


    /***
     * If initial sensor measurement is RADAR, then convert from polar
     * coordinates to Cartesian coordinates because the state vector
     * is maintained in Cartesian coordinates.   Set velocity to zero.
     *
     * Otherwise if sensor measurement is LASER, then extract position
     * information and set velocity to zero.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
     /**
      * Initialize state from RADAR
      */
      rho = measurement_pack.raw_measurements_(0);
      phi = measurement_pack.raw_measurements_(1);
      rhodot = measurement_pack.raw_measurements_(3); 
      
      px = rho * cos(phi);
      py = rho * sin(phi);
      //vx = 5.199937;
      vx = 0;
      vy = 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
     /**
      * Initialize state from LASER
      */
      px = measurement_pack.raw_measurements_(0);
      py = measurement_pack.raw_measurements_(1);
      //vx = 5.199937;
      //vy = 0.00179685;
      vx = 0;
      vy = 0;
    }

   /***
    * Numerical precision adustment
    * Snap px to reasonable small number (1/10,000th currenty)
    * Snap py to reasonable small number
    */
    if (fabs(px) < SMALLNUM ) {
      px = SMALLNUM;
    }

    if (fabs(py) < SMALLNUM ) {
      py = SMALLNUM;
    }


    ekf_.x_ << px, py, vx, vy;

    previous_timestamp_ = measurement_pack.timestamp_;

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

  /***
   * Calculate time lapse
   */
  float dt = ( (float) (measurement_pack.timestamp_ - previous_timestamp_))/ 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;


  if (DEBUG) {
    float diff = (measurement_pack.timestamp_ - previous_timestamp_);

    cout << "FusionEKF::ProcessMeasurement:  timestamp= " << measurement_pack.timestamp_ << endl;
    cout << "FusionEKF::ProcessMeasurement:  dt= " << dt << endl;
    cout << "FusionEKF::ProcessMeasurement:  diff= " << diff << endl;
  }

  float dt_2= dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4= dt_3 * dt;

  /***
   * Populate the state transition matrix, F, accounting
   * for time lapse.
   */
  if (DEBUG)
    cout << "FusionEKF::ProcessMeasurement:  update state transition matrix F" << endl;

  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  if (DEBUG)
    cout << "ekf_.F_ = " << ekf_.F_ << endl;

  
 /***
  * Update the process covariance matrix, Q
  * 
  * Note: noise_ax and noise_ay were set in constructor (shoud it be updated?)
  */

  if (DEBUG)
    cout << "FusionEKF::ProcessMeasurement:  update process covariance, Q...." << endl;

  ekf_.Q_ <<  (dt_4/4.0)*noise_ax, 0, (dt_3/2.0)*noise_ax, 0,
  	      0, (dt_4/4.0)*noise_ay, 0.0, (dt_3/2.0)*noise_ay,
	      (dt_3/2.0)*noise_ax, 0.0, (dt_2*noise_ax), 0.0,
	      0.0, (dt_3/2.0)*noise_ay, 0.0, dt_2*noise_ay;

  if (DEBUG) {
    cout << "FusionEKF::ProcessMeasurement:  prediciton step, x= Fx + Q" << endl;
    cout << "FusionEKF::ProcessMeasurement:  ekf_.F_= " << ekf_.F_ << endl;
    cout << "FusionEKF::ProcessMeasurement:  ekf_.Q_= " << ekf_.Q_ << endl;
  }

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  if (DEBUG)
    cout << "FusionEKF::ProcessMeasurement:  update process noice covariance" << endl;

  /***
   * Math uses LaTex like notation...
   *
   * If the measurement update is from RADAR, calculate the Jacobian used for the
   * measurement update along with the uncertaintly matrix, R.  
   *
   * The current state variable, x=(px, py, vx, vy) is maintained since
   * this prediction is compared to linearized RADAR measuremnt z= (rho, phi, rhodot)
   * in calculating the Kalman gain.  This is done by a transform function h(x)
   * when computing the error term.   The Jacobian Hk is used for projection 
   * S= H P' H^{i} + R  and  K = P' H^t S^{-1} in computing kalman gain, K.
   *
   * If the measurement update is from LASER, use the transform matrix H and R_laser.
   * For this case, the Kalman gain will use the H matrix.
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates

    if (DEBUG)
      cout << "FusionEKF::ProcessMeasurement:  Radar measurment update..." << endl;

     Hj_ = tools.CalculateJacobian(ekf_.x_);

     //for a Jacobian that is degenerate, do not
     //perform the update

     if(Hj_.isZero(0)) {
        cout << "FusionEKF::ProcessMeasurement:  degenerate Jacobian Hj, skipping update..." << endl;
	return;
     } 

    /***
     * now upbade the EKF measurement matrix with
     * the linearizing Jacobian Hj_ and the measurement
     * uncertainty (covariance) matrix R_radar_
     */
     ekf_.H_ = Hj_;
     ekf_.R_ = R_radar_;

     if (DEBUG)
       cout << "FusionEKF::ProcessMeasurement:  performing EKF measurment update..." << endl;

     ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates

    if (DEBUG)
      cout << "FusionEKF::ProcessMeasurement:  Laser measurment update..." << endl;

    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;


    if (DEBUG)
      cout << "FusionEKF::ProcessMeasurement:  performing measurment update..." << endl;

    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;

}
