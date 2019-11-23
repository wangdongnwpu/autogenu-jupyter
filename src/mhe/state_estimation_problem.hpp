// Abstruct class for state estimation problems for moving horizon estimation. 
// Loads model of MHE and define dimensions.

#ifndef STATE_ESTIMATION_PROBLEM_H
#define STATE_ESTIMATION_PROBLEM_H

#include "mhe_model.hpp"

// Abstruct class for state estimation problems. This class loads model of MHE  
// and define dimensions.
class StateEstimationProblem {
public:
  // Loads model of MHE and define dimensions.
  StateEstimationProblem();
  virtual ~StateEstimationProblem() = default;

  // Prohibits copy.
  StateEstimationProblem(const StateEstimationProblem&) = delete;
  StateEstimationProblem& operator=(const StateEstimationProblem&) = delete;

  // Returns dimension of the state.
  int dim_state() const;

  // Returns dimension of the control input.
  int dim_control_input() const;

  // Returns dimension of the constraints.
  int dim_constraints() const;

  // Returns dimension of the output.
  int dim_output() const;

  // Returns dimension of the disturbance.
  int dim_disturbance() const;

  // Returns dimension of the solution of the state estimation problem.
  virtual int dim_solution() const = 0;

protected:
  MHEModel model_;
  int dim_state_, dim_control_input_, dim_output_, dim_constraints_, 
      dim_disturbance_, dim_solution_;
};

#endif // STATE_ESTIMATION_PROBLEM_H