#include "zero_horizon_disturbance_estimation_problem.hpp"

ZeroHorizonDisturbanceEstimationProblem::
ZeroHorizonDisturbanceEstimationProblem()
  : StateEstimationProblem(),
    dim_solution_(model_.dimConstraints()+model_.dimDisturbance()),
    lambda_vec_(linearalgebra::NewVector(model_.dimState())) {
}

ZeroHorizonDisturbanceEstimationProblem::
~ZeroHorizonDisturbanceEstimationProblem() {
  linearalgebra::DeleteVector(lambda_vec_);
}

void ZeroHorizonDisturbanceEstimationProblem::computeOptimalityResidual(
    const double time, const double* measured_control_input_vec, 
    const double* measured_output_vec, const double* estimated_state_vec, 
    const double* solution_vec, double* optimality_residual) {
  model_.phixFunc(time, estimated_state_vec, lambda_vec_);
  for (int i=0; i<dim_state_; ++i) {
    lambda_vec_[i] *= -1;
  }
  model_.hwFunc(
      time, estimated_state_vec, measured_control_input_vec, 
      measured_output_vec, solution_vec, &(solution_vec[dim_disturbance_]),
      lambda_vec_, optimality_residual);
  model_.constraintFunc(
      time, estimated_state_vec, measured_control_input_vec, 
      measured_output_vec, solution_vec, &(solution_vec[dim_disturbance_]),
      &(optimality_residual[dim_disturbance_]));
}

int ZeroHorizonDisturbanceEstimationProblem::dim_solution() const {
  return dim_solution_;
}