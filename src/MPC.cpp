#include "MPC.h"
#include <math.h>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// Setting constants
const size_t N = 10;
const double dt = 0.1;
const double REFERENCE_V = 80;  // MPH

// Penalty weight constants for mpc tuning.
// The most important thing is to stay on track -- thus keep cte and epsi minimized.
const double REFERENCE_CTE_WEIGHT = 7000.0;
const double REFERENCE_EPSI_WEIGHT = 7000.0;
// Okay, these terms are not that important, actually
const double REFERENCE_V_WEIGHT = 5.0;
const double ACTUATOR_DELTA_WEIGHT = 1.0;
const double ACTUATOR_A_WEIGHT = 1.0;
// Penalizing delta change in order to avoid sudden movements on the steering wheel
const double CHANGE_DELTA_WEIGHT = 100.0;
// Okay, but we might want to break suddenly, so the penalty here is not as high as for the delta change
const double CHANGE_A_WEIGHT = 50.0;


// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
const size_t x_start = 0;
const size_t y_start = x_start + N;
const size_t psi_start = y_start + N;
const size_t v_start = psi_start + N;
const size_t cte_start = v_start + N;
const size_t epsi_start = cte_start + N;
const size_t delta_start = epsi_start + N;
const size_t a_start = delta_start + N - 1;


// The constant 'Lf' is defined in main.cpp
extern const double Lf;


class FG_eval {

 public:

  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {

    // defining 't' only once
    size_t t;

    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // The cost is stored is the first element of `fg`. Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // The part of the cost based on the reference state. Reference values:
    //  - CTE: 0
    //  - EPSI: 0
    //  - V: as defined by REFERENCE_V
    for ( t = 0; t < N; t++ ) {
      fg[0] += REFERENCE_CTE_WEIGHT * CppAD::pow(vars[cte_start + t], 2);
      fg[0] += REFERENCE_EPSI_WEIGHT * CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += REFERENCE_V_WEIGHT * CppAD::pow(vars[v_start + t] - REFERENCE_V, 2);

      // Don't drive too fast while being off track
      // fg[0] += 500 * CppAD::pow(vars[v_start + t] * vars[cte_start + t], 2);
    }

    // Minimize the use of actuators.
    for ( t = 0; t < N - 1; t++ ) {
      fg[0] += ACTUATOR_DELTA_WEIGHT * CppAD::pow(vars[delta_start + t], 2);
      fg[0] += ACTUATOR_A_WEIGHT * CppAD::pow(vars[a_start + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for ( t = 0; t < N - 2; t++ ) {
      fg[0] += CHANGE_DELTA_WEIGHT * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += CHANGE_A_WEIGHT * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    // Setup Constraints

    // Initial constraints.
    // We add 1 to each of the starting indices due to cost being located at index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // The rest of the constraints
    for ( t = 1; t < N; t++ ) {

      // The state at time t+1
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];

      // The state at time t
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      // Only consider the actuation at time t
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
      AD<double> f0_prime = coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0, 2);
      AD<double> psides0 = CppAD::atan(f0_prime);

      // Equations for the model:
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] - v[t] / Lf * delta[t] * dt               !!! mind the operator: '-'
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] - v[t] * delta[t] / Lf * dt   !!! mind the operator: '-'
      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 - v0 * delta0 / Lf * dt);
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
      fg[1 + cte_start + t] = cte1 - (f0 - y0 + v0 * CppAD::sin(epsi0) * dt);
      fg[1 + epsi_start + t] = epsi1 - (psi0 - psides0 - v0 * delta0 / Lf * dt);
    }
  }
};


//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {

  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  // Setting the number of model variables (includes both states and inputs).
  size_t n_vars = 6 * N + 2 * (N - 1);

  // Setting the number of constraints
  size_t n_constraints = N * 6;

  // The good old loop index variable
  size_t i;

  // Initial value of the independent variables. SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for ( i = 0; i < n_vars; i++ ) {
    vars[i] = 0.0;
  }

  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Setting upper and lowerlimits for the variables.

  // Set all non-actuators upper and lowerlimits to the max negative and positive values.
  for ( i = 0; i < delta_start; i++ ) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25 degrees (values in radians).
  for ( i = delta_start; i < a_start; i++ ) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  for ( i = a_start; i < n_vars; i++ ) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for ( i = 0; i < n_constraints; i++ ) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  if ( !ok ) {
    // ...
    std::cerr << "*** Something went wrong. :-(" << std::endl;
  }

  // Cost
  // auto cost = solution.obj_value;
  // std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with `solution.x[i]`.
  vector<double> result;

  // steering angle and throttle in the first place...
  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);

  // ...then the coordinates of the trajectory line
  for ( i = 0; i < N; i++ ) {

    result.push_back(solution.x[x_start + i]);
    result.push_back(solution.x[y_start + i]);
  }

  return result;
}
