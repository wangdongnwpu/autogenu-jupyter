#include "state_estimate_continuation.hpp"

StateEstimateContinuation::StateEstimateContinuation(
    const double T_f, const double alpha, const int N,
    const double finite_difference_increment, const double zeta)
  : estimation_problem_(T_f, alpha, N),
    dim_state_(estimation_problem_.dim_state()),
    dim_control_input_(estimation_problem_.dim_control_input()), 
    dim_constraints_(estimation_problem_.dim_constraints()), 
    dim_disturbance_(estimation_problem_.dim_disturbance()), 
    dim_output_(estimation_problem_.dim_output()), 
    dim_solution_(estimation_problem_.dim_solution()),
    finite_difference_increment_(finite_difference_increment),
    zeta_(zeta),
    incremented_time_(0),
    incremented_solution_vec_(linearalgebra::NewVector(dim_solution_)),
    optimality_residual_(linearalgebra::NewVector(dim_solution_)),
    optimality_residual_1_(linearalgebra::NewVector(dim_solution_)),
    optimality_residual_2_(linearalgebra::NewVector(dim_solution_)) {
}

StateEstimateContinuation::StateEstimateContinuation(
    const double T_f, const double alpha, const int N, 
    const double initial_time, const double finite_difference_increment, 
    const double zeta)
  : estimation_problem_(T_f, alpha, N),
    dim_state_(estimation_problem_.dim_state()),
    dim_control_input_(estimation_problem_.dim_control_input()), 
    dim_constraints_(estimation_problem_.dim_constraints()), 
    dim_disturbance_(estimation_problem_.dim_disturbance()), 
    dim_output_(estimation_problem_.dim_output()), 
    dim_solution_(estimation_problem_.dim_solution()),
    finite_difference_increment_(finite_difference_increment),
    zeta_(zeta),
    incremented_time_(initial_time),
    incremented_solution_vec_(linearalgebra::NewVector(dim_solution_)),
    optimality_residual_(linearalgebra::NewVector(dim_solution_)),
    optimality_residual_1_(linearalgebra::NewVector(dim_solution_)),
    optimality_residual_2_(linearalgebra::NewVector(dim_solution_)) {
}

StateEstimateContinuation::~StateEstimateContinuation() {
  linearalgebra::DeleteVector(incremented_solution_vec_);
  linearalgebra::DeleteVector(optimality_residual_);
  linearalgebra::DeleteVector(optimality_residual_1_);
  linearalgebra::DeleteVector(optimality_residual_2_);
}

void StateEstimateContinuation::integrateSolution(
    double* solution_vec, const double* solution_update_vec, 
    const double integration_length) {
  for (int i=0; i<dim_solution_; ++i) {
    solution_vec[i] += integration_length * solution_update_vec[i];
  }
}

double StateEstimateContinuation::computeErrorNorm(
    const double time, double const* const* control_input_mat,
    double const* const* outpu_mat, const double* solution_vec) {
  estimation_problem_.computeOptimalityResidual(time, control_input_mat,
                                                outpu_mat, solution_vec, 
                                                optimality_residual_);
  return std::sqrt(linearalgebra::SquaredNorm(dim_solution_, 
                                              optimality_residual_));
}

void StateEstimateContinuation::resetHorizonLength(const double T_f, 
                                                   const double alpha, 
                                                   const double initial_time) {
  estimation_problem_.resetHorizonLength(T_f, alpha, initial_time);
}

void StateEstimateContinuation::resetHorizonLength(const double initial_time) {
  estimation_problem_.resetHorizonLength(initial_time);
}

void StateEstimateContinuation::bFunc(
    const double time, double const* const* control_input_mat,
    double const* const* output_mat, 
    double const* const* incremented_control_input_mat,
    double const* const* incremented_output_mat, 
    const double* current_solution_vec, 
    const double* current_solution_update_vec, double* b_vec) {
  incremented_time_ = time + finite_difference_increment_;
  for (int i=0; i<dim_solution_; ++i) {
    incremented_solution_vec_[i] = current_solution_vec[i] 
        + finite_difference_increment_ * solution_update_vec_[i];
  }
  estimation_problem_.computeOptimalityResidual(time, control_input_mat, 
                                                output_mat,
                                                current_solution_vec, 
                                                optimality_residual_);
  estimation_problem_.computeOptimalityResidual(incremented_time_, 
                                                incremented_control_input_mat,
                                                incremented_output_mat, 
                                                current_solution_vec, 
                                                optimality_residual_1_);
  estimation_problem_.computeOptimalityResidual(incremented_time_, 
                                                incremented_control_input_mat, 
                                                incremented_output_mat, 
                                                incremented_solution_vec_, 
                                                optimality_residual_2_);
  for (int i=0; i<dim_solution_; ++i) {
    b_vec[i] = (1/finite_difference_increment_-zeta_) * optimality_residual_[i] 
        - optimality_residual_2_[i] / finite_difference_increment_;
  }
}

void StateEstimateContinuation::AxFunc(
    const double time, double const* const* control_input_mat,
    double const* const* output_mat, 
    double const* const* incremented_control_input_mat,
    double const* const* incremented_output_mat, 
    const double* current_solution_vec, const double* direction_vec,
    double* ax_vec) {
  for (int i=0; i<dim_solution_; ++i) {
    incremented_solution_vec_[i] = current_solution_vec[i] 
        + finite_difference_increment_ * direction_vec[i];
  }
  estimation_problem_.computeOptimalityResidual(incremented_time_, 
                                                incremented_control_input_mat, 
                                                incremented_output_mat, 
                                                incremented_solution_vec_, 
                                                optimality_residual_2_);
  for (int i=0; i<dim_solution_; ++i) {
    ax_vec[i] = 
        (optimality_residual_2_[i]-optimality_residual_1_[i]) 
        / finite_difference_increment_;
  }
}

int StateEstimateContinuation::dim_state() const {
  return dim_state_;
}

int StateEstimateContinuation::dim_control_input() const {
  return dim_control_input_;
}

int StateEstimateContinuation::dim_constraints() const {
  return dim_constraints_;
}

int StateEstimateContinuation::dim_disturbance() const {
  return dim_disturbance_;
}

int StateEstimateContinuation::dim_output() const {
  return dim_output_;
}

int StateEstimateContinuation::dim_solution() const {
  return dim_solution_;
}

int StateEstimateContinuation::N() const {
  return estimation_problem_.N();
}