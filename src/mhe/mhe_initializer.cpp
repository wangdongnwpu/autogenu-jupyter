#include "mhe_initializer.hpp"

MHEInitializer::MHEInitializer(const double finite_difference_increment, 
                               const int kmax, 
                               const double newton_residual_tolerance, 
                               const int max_newton_iteration) 
  : newton_(finite_difference_increment),
    mfgmres_(newton_.dim_solution(), kmax),
    dim_constraints_(newton_.dim_constraints()),
    dim_disturbance_(newton_.dim_disturbance()),
    dim_solution_(newton_.dim_solution()),
    max_newton_iteration_(max_newton_iteration),
    newton_residual_tolerance_(newton_residual_tolerance),
    initial_guess_solution_vec_(
        linearalgebra::NewVector(newton_.dim_solution())),
    solution_update_vec_(
        linearalgebra::NewVector(newton_.dim_solution())) {
}


MHEInitializer::MHEInitializer(const double finite_difference_increment, 
                               const int kmax)
  : newton_(finite_difference_increment),
    mfgmres_(newton_.dim_solution(), kmax),
    dim_constraints_(newton_.dim_constraints()),
    dim_disturbance_(newton_.dim_disturbance()),
    dim_solution_(newton_.dim_solution()),
    max_newton_iteration_(50),
    newton_residual_tolerance_(1e-08),
    initial_guess_solution_vec_(
        linearalgebra::NewVector(newton_.dim_solution())),
    solution_update_vec_(
        linearalgebra::NewVector(newton_.dim_solution())) {
}

MHEInitializer::~MHEInitializer() {
  linearalgebra::DeleteVector(initial_guess_solution_vec_);
  linearalgebra::DeleteVector(solution_update_vec_);
}

void MHEInitializer::setCriterionsOfNewtonTermination(
    const double newton_residual_tolerance, const int max_newton_iteration) {
  newton_residual_tolerance_ = newton_residual_tolerance;
  max_newton_iteration_ = max_newton_iteration;
}

void MHEInitializer::setInitialGuessSolution(
    const double* initial_guess_solution_vec) {
  for (int i=0; i<dim_solution_; ++i) {
    initial_guess_solution_vec_[i] = initial_guess_solution_vec[i];
  }
}

void MHEInitializer::computeInitialSolution(
    const double initial_time, const double* initial_measured_control_input, 
    const double* initial_measured_output, 
    const double* initial_estimated_state_vec, double* initial_solution_vec) {
  for (int i=0; i<dim_solution_; ++i) {
    initial_solution_vec[i] = initial_guess_solution_vec_[i];
  }
  int num_itr = 0;
  double optimality_error = newton_.errorNorm(initial_time, 
                                              initial_measured_control_input,
                                              initial_measured_output,
                                              initial_estimated_state_vec,
                                              initial_solution_vec);
  while (optimality_error > newton_residual_tolerance_ 
         && num_itr < max_newton_iteration_) {
    mfgmres_.solveLinearProblem(newton_, initial_time, 
                                initial_measured_control_input,
                                initial_measured_output,
                                initial_estimated_state_vec,
                                initial_solution_vec, solution_update_vec_);
    for (int i=0; i<dim_solution_; ++i) {
      initial_solution_vec[i] += solution_update_vec_[i];
    }
    optimality_error = newton_.errorNorm(initial_time, , 
                                         initial_measured_control_input,
                                         initial_measured_output,
                                         initial_estimated_state_vec,
                                         initial_solution_vec);
    ++num_itr;
  }
}

int MHEInitializer::dim_solution() const {
  return dim_solution_;
}