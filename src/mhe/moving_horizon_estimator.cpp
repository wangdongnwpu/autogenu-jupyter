#include "moving_horizon_estimator.hpp"

MovingHorizonEstimator::MovingHorizonEstimator(
    const double T_f, const double alpha, const int N,
    const double finite_difference_increment, const double zeta, 
    const int kmax, const double min_sampling_period)
  : continuation_problem_(T_f, alpha, N, finite_difference_increment, zeta),
    mfgmres_(continuation_problem_.dim_solution(), kmax), 
    control_input_data_(
        continuation_problem_.dim_control_input(), T_f/min_sampling_period), 
    output_data_(
        continuation_problem_.dim_output(), T_f/min_sampling_period), 
    horizon_(T_f, alpha),
    solution_initializer_(),
    dim_state_(continuation_problem_.dim_state()), 
    dim_control_input_(continuation_problem_.dim_control_input()), 
    dim_constraints_(continuation_problem_.dim_constraints()), 
    dim_disturbance_(continuation_problem_.dim_disturbance()), 
    dim_output_(continuation_problem_.dim_output()),
    N_(N),
    finite_difference_increment_(finite_difference_increment),
    solution_vec_(
        linearalgebra::NewVector(continuation_problem_.dim_solution())), 
    solution_update_vec_(
        linearalgebra::NewVector(continuation_problem_.dim_control_input())), 
    initial_solution_vec_(
        linearalgebra::NewVector(solution_initializer_.dim_solution())), 
    control_input_mat_(
        linearalgebra::NewMatrix(N+1, 
                                 continuation_problem_.dim_control_input())), 
    output_mat_(
        linearalgebra::NewMatrix(N+1, continuation_problem_.dim_output())), 
    incremented_control_input_mat_(
        linearalgebra::NewMatrix(N+1, 
                                 continuation_problem_.dim_control_input())),
    incremented_output_mat_(
        linearalgebra::NewMatrix(N+1, continuation_problem_.dim_output())) {
}

MovingHorizonEstimator::~MovingHorizonEstimator() {
  linearalgebra::DeleteVector(solution_vec_);
  linearalgebra::DeleteVector(solution_update_vec_);
  linearalgebra::DeleteVector(initial_solution_vec_);
  linearalgebra::DeleteMatrix(control_input_mat_);
  linearalgebra::DeleteMatrix(output_mat_);
  linearalgebra::DeleteMatrix(incremented_control_input_mat_);
  linearalgebra::DeleteMatrix(incremented_output_mat_);
}

void MovingHorizonEstimator::updateEstimation(const double time, 
                                              double* estimated_state_vec,
                                              const double sampling_period) {
  double horizon_length = horizon_.getLength(time);
  double incremented_time = time + finite_difference_increment_;
  double incremented_horizon_length = horizon_.getLength(incremented_time);

  control_input_data_.generateUniformlySpacedTimeSeriesData(
      time, horizon_length, N_, control_input_mat_);
  control_input_data_.generateUniformlySpacedTimeSeriesData(
      incremented_time, incremented_horizon_length, N_, 
      incremented_control_input_mat_);
  output_data_.generateUniformlySpacedTimeSeriesData(
      time, horizon_length, N_, output_mat_);
  output_data_.generateUniformlySpacedTimeSeriesData(
      incremented_time, incremented_horizon_length, N_, 
      incremented_output_mat_);
  mfgmres_.solveLinearProblem(continuation_problem_, time, control_input_data_, 
                              output_mat_, incremented_control_input_mat_, 
                              incremented_output_mat_, solution_vec_, 
                              solution_update_vec_);
  continuation_problem_.integrateSolution(solution_vec_, solution_update_vec_, 
                                          sampling_period);
  for (int i=0; i<dim_state_; ++i) {
    estimated_state_vec[i] 
        = solution_vec_[N_*(dim_constraints_+dim_disturbance_)+i];
  }
}

void MovingHorizonEstimator::getStateEstimation(
    double* estimated_state_vec) const {
  for (int i=0; i<dim_state_; ++i) {
    estimated_state_vec[i] 
        = solution_vec_[N_*(dim_constraints_+dim_disturbance_)+i];
  }
}

void MovingHorizonEstimator::setParametersForInitialization(
    const double* initial_guess_solution, 
    const double newton_residual_tolerance, const int max_newton_iteration) {
  solution_initializer_.setInitialGuessSolution(initial_guess_solution);
  solution_initializer_.setCriterionsOfNewtonTermination(
      newton_residual_tolerance, max_newton_iteration);
}

void MovingHorizonEstimator::initializeSolution(
    const double initial_time, const double* initial_control_input_vec,
    const double* initial_output_vec, 
    const double* initial_estimated_state_vec) {
  solution_initializer_.computeInitialSolution(
      initial_time, initial_control_input_vec, initial_output_vec, 
      initial_estimated_state_vec, initial_solution_vec_);
  for (int i=0; i<continuation_problem_.N(); ++i) {
    for (int j=0; j<dim_constraints_+dim_disturbance_; ++j) {
      solution_vec_[i*(dim_constraints_+dim_disturbance_)+j]
          = initial_solution_vec_[j];
    }
  }
  for (int i=0; i<dim_state_; ++i) {
    solution_vec_[continuation_problem_.N()*(dim_constraints_+dim_disturbance_)+i]
        = initial_estimated_state_vec[i];
  }
  control_input_data_.appendData(initial_time, initial_control_input_vec);
  output_data_.appendData(initial_time, initial_output_vec);
  continuation_problem_.resetHorizonLength(initial_time);
}

double MovingHorizonEstimator::getErrorNorm(const double time) {
  control_input_data_.generateUniformlySpacedTimeSeriesData(
      time, horizon_.getLength(time), N_, control_input_mat_);
  output_data_.generateUniformlySpacedTimeSeriesData(
      time, horizon_.getLength(time), N_, output_mat_);
  return continuation_problem_.computeErrorNorm(time, control_input_mat_, 
                                                output_mat_, solution_vec_);
}

void MovingHorizonEstimator::appendMeasuredControlInput(
    const double measured_time, const double* measured_control_input) {
  control_input_data_.appendData(measured_time, measured_control_input);
}

void MovingHorizonEstimator::appendMeasuredOutput(
    const double measured_time, const double* measured_output) {
  output_data_.appendData(measured_time, measured_output);
}