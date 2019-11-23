#include "newton_gmres_for_disturbance_estimation.hpp"

NewtonGMRESForDisturbanceEstimation::NewtonGMRESForDisturbanceEstimation(
    const double finite_difference_increment) 
  : problem_(), 
    dim_solution_(problem_.dim_solution()),
    finite_difference_increment_(finite_difference_increment),
    incremented_solution_vec_(linearalgebra::NewVector(problem_.dim_solution())),
    optimality_residual_(linearalgebra::NewVector(problem_.dim_solution())),
    optimality_residual_1_(linearalgebra::NewVector(problem_.dim_solution())) {
}

NewtonGMRESForDisturbanceEstimation::~NewtonGMRESForDisturbanceEstimation() {
  linearalgebra::DeleteVector(incremented_solution_vec_);
  linearalgebra::DeleteVector(optimality_residual_);
  linearalgebra::DeleteVector(optimality_residual_1_);
}

double NewtonGMRESForDisturbanceEstimation::errorNorm(
    const double time, const double* measured_control_input_vec,
    const double* measured_output_vec, const double* estimated_state_vec, 
    const double* solution_vec) {
  problem_.computeOptimalityResidual(time, measured_control_input_vec, 
                                     measured_output_vec, estimated_state_vec, 
                                     solution_vec, optimality_residual_);
  return std::sqrt(linearalgebra::SquaredNorm(dim_solution_, 
                                              optimality_residual_));
}


void NewtonGMRESForDisturbanceEstimation::bFunc(
    const double time, const double* measured_control_input_vec,
    const double* measured_output_vec, const double* estimated_state_vec, 
    const double* current_solution_vec, 
    const double* current_solution_update_vec, double* b_vec) {
  for (int i=0; i<dim_solution_; ++i) {
    incremented_solution_vec_[i] = current_solution_vec[i] 
        + finite_difference_increment_ * current_solution_update_vec[i];
  }
  problem_.computeOptimalityResidual(time, measured_control_input_vec, 
                                      measured_output_vec, estimated_state_vec,
                                      current_solution_vec, 
                                      optimality_residual_);
  problem_.computeOptimalityResidual(time, measured_control_input_vec, 
                                      measured_output_vec, estimated_state_vec,
                                      incremented_solution_vec_, 
                                      optimality_residual_1_);
  for (int i=0; i<dim_solution_; ++i) {
    b_vec[i] = - optimality_residual_[i] 
        - (optimality_residual_1_[i]-optimality_residual_[i]) 
        / finite_difference_increment_;
  }
}

void NewtonGMRESForDisturbanceEstimation::AxFunc(
    const double time, const double* measured_control_input_vec,
    const double* measured_output_vec, const double* estimated_state_vec, 
    const double* current_solution_vec, const double* direction_vec,
    double* ax_vec) {
  for (int i=0; i<dim_solution_; ++i) { 
    incremented_solution_vec_[i] = current_solution_vec[i] 
        + finite_difference_increment_ * direction_vec[i];
  }
  problem_.computeOptimalityResidual(time, measured_control_input_vec, 
                                      measured_output_vec, estimated_state_vec,
                                      incremented_solution_vec_, 
                                      optimality_residual_1_);
  for (int i=0; i<dim_solution_; ++i) {
    ax_vec[i] = (optimality_residual_1_[i]-optimality_residual_[i]) 
        / finite_difference_increment_;
  }
}

int NewtonGMRESForDisturbanceEstimation::dim_state() const {
  return problem_.dim_state();
}

int NewtonGMRESForDisturbanceEstimation::dim_control_input() const {
  return problem_.dim_control_input();
}

int NewtonGMRESForDisturbanceEstimation::dim_constraints() const {
  return problem_.dim_constraints();
}

int NewtonGMRESForDisturbanceEstimation::dim_output() const {
  return problem_.dim_output();
}

int NewtonGMRESForDisturbanceEstimation::dim_disturbance() const {
  return problem_.dim_disturbance();
}

int NewtonGMRESForDisturbanceEstimation::dim_solution() const {
  return problem_.dim_solution();
}