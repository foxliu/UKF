#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);

  rmse.fill(0.0);

  unsigned long est_size = estimations.size();
  unsigned long ground_size = ground_truth.size();

  if (est_size != ground_size || est_size == 0)
  {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  for (int i = 0; i < est_size; ++i)
  {
    VectorXd residual = estimations[i] - ground_truth[i];

    residual = residual.array() * residual.array();
    rmse += residual;
  }

  rmse = rmse / est_size;
  rmse = rmse.array().sqrt();
  return rmse;
}