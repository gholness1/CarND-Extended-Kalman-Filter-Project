#include "kalman_filter.h"
#include "debug.h"
#include <iostream>
#include<math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  /**
    TODO:
    * predict the state
  */

  //Kalman Filter linear prediction using process uncertainty, Q_
  x_ = F_ *x_;

  MatrixXd Ft= F_.transpose();

  P_ = F_ * P_ * Ft + Q_;
}


/****
 * KalmanFilter::h:
 *
 * This function computes h(x) function used to convert cartesian
 * state 4-vector, x= (px,py,vx,vy) to polar coordinates h(x) = (rho,phi,rhodot)
 *
 * Eigen::VectorXd &x_prime-  state vector
 *
 * returns:
 *   VectorXd- 3-vector in polar coordinates, (rho,phi,rhodot)
 */

VectorXd KalmanFilter::h(const Eigen::VectorXd &x_prime) {
  VectorXd hx= VectorXd(3);

  float rho= 0.0;
  float phi= 0.0;
  float rhodot= 0.0;

  float px= x_prime[0];
  float py= x_prime[1];
  float vx= x_prime[2];
  float vy= x_prime[3];


  /**
   * Numerical precision adustment
   * Snap px to reasonable unit of 1/1000th
   * Snap py to reasonable unit of 1/1000th
   */
   if (fabs(px) < (1.0/1000.0) ) {
     px = 0.001;
   }

   if (fabs(py) < (1.0/1000.0) ) {
     py = 0.001;
   }

   rho = sqrt(px*px + py*py);

   if (rho < (1.0/1000.0)) {
     rho = 0.001;
   }

   float p_slope= py/px;
     
   //define zero slope to be 1/10,000th radians
   if ( abs(p_slope) < 1.0/1000) {
     phi= 0.0;
     rhodot= (px*vx + py*vy)/rho;
   } else {
     rho = sqrt(px*px + py*py);
     phi = atan2(py,px);
     rhodot= (px*vx + py*vy)/rho;
   } 

   #if 0
   if (phi > M_PI) {
     phi = phi - 2*M_PI;
   } else if (phi < -M_PI) {
     phi = phi + 2*M_PI;
   }
   #endif

   hx << rho, phi, rhodot;

   return hx;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  if (DEBUG) {
    cout << "KalmanFilter::Update on z= " << z << endl;
    cout << "KalmanFilter::Update x_= " << x_ << endl;
  }

  VectorXd y = z - H_ * x_;
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;
 
  x_ = x_ + (K *y);
  P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  if (DEBUG) {
    cout << "KalmanFilter::UpdateEKF on " << z << endl;
    cout << "KalmanFilter::UpdateEKF x_= " << x_ << endl;
  }

  VectorXd y = z - h(x_);

  float computed_phi = y(1);

  if (DEBUG)
    cout << "KalmanFilter::UpdateEKF:  phi, y(1), computed_phi= " << computed_phi << endl;

  /**
   * Transform angle, phi, within y=(rho,phi,rhodot) to be in closed interval [-PI, PI]
   *
   * Note: M_PI is from <math.h>
   */
  float quotient = 0.0;
  float div= 0.0;
  float remainder = 0.0;

  if (computed_phi > M_PI) {
    div = computed_phi/(2*M_PI);
    quotient = floor(computed_phi/(2*M_PI));
    remainder = computed_phi * (div - quotient);
    
    computed_phi = remainder;

  } else if (computed_phi < -M_PI) {
    div = -computed_phi/(2*M_PI);
    quotient = floor(computed_phi/(2*M_PI));
    remainder = computed_phi * (div - quotient);
    
    computed_phi = -remainder;
  }
 
  y(1) = computed_phi;
 
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;
      
  x_ = x_ + (K *y); 
  P_ = (I_ - K * H_) * P_;
}
