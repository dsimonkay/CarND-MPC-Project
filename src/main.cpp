#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "MPC.h"

// for convenience
using json = nlohmann::json;

// defining the latency of the actuators (in milliseconds)
const int ACTUATOR_LATENCY = 100;

// reference speed (MPH)
const double REFERENCE_V = 90.0;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
extern const double Lf = 2.67;


// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }


// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}


// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {

  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }

  return result;
}


// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order) {

  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}


// Main program module
int main() {

  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  // Debug info
  std::cout << std::endl << "Actuator latency: " << ACTUATOR_LATENCY << " milliseconds" << std::endl << std::endl;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {

    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    // cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {

      string s = hasData(sdata);
      if (s != "") {

        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {

          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double delta = j[1]["steering_angle"];
          double a = j[1]["throttle"];

          // Debug output
          std::cout << "[x, y, psi, v, delta, a]: " << px << ", " << py <<
                       ", " << psi << ", " << v << ", " << delta << ", " << a;

          // Common index variable
          size_t i;

          // Transforming waypoint coordinates to coordinates in the vehicle's coordinate system
          Eigen::VectorXd ptsx_car(ptsx.size());
          Eigen::VectorXd ptsy_car(ptsy.size());
          double waypoint_x, waypoint_y;
          for( i = 0;  i < ptsx.size();  i++ ) {

            // The map coordinates of the ith waypoint relative to the car
            waypoint_x = ptsx[i] - px;
            waypoint_y = ptsy[i] - py;

            // Performing the homogeneous transformation from the car's point of view (=base is (0,0))
            ptsx_car[i] = 0.0 + cos(-psi) * waypoint_x - sin(-psi) * waypoint_y;
            ptsy_car[i] = 0.0 + sin(-psi) * waypoint_x + cos(-psi) * waypoint_y;
          }

          // Fitting the polynomial
          auto coeffs = polyfit(ptsx_car, ptsy_car, 3);

          // Defining the "original" initial state. We're in the car's coordinate system now:
          // the vehicle is located at (0, 0).
          double x = 0.0;
          double y = 0.0;
          psi = 0.0;
          // The initial value of 'v' is known, v = v

          // polyeval(coeffs, 0) delivers simply coeffs[0], because x = 0
          double cte = coeffs[0]; // polyeval(coeffs, 0); 

          // f'(x) = coeffs[1] + 2*coeffs[2]*x + 3*coeffs[3]*x^2  ==>  for x = 0 only coeffs[1] remains
          double epsi = -atan(coeffs[1]);

          // Debug output
          std::cout << "    [cte, epsi]: " << cte << ", " << epsi;

          // The original initial state vector would look like:
          // state << 0.0, 0.0, 0.0, v, cte, epsi;

          // Considering actuator latency: let the vehicle continue driving for the specified amount of time,
          // then update the state variables and run the solver with the new "initial" state.
          const double latency = ACTUATOR_LATENCY / 1000.0;  // in seconds

          // Determining the new initial state variables
          x += v * latency;                     // _[t+1] = x[t] + v[t] * cos(psi[t]) * dt
          // y remains unchanged, as there's no motion along the Y axis in the car's coordinate system.
          // In other words: y[t+1] = y[t] + v[t] * sin(psi[t]) ==> equals to zero because sin(psi) = 0
          psi -= v / Lf * delta * latency;      // psi_[t+1] = psi[t] - v[t] / Lf * delta[t] * dt

          cte += v * sin(epsi) * latency;       // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
          epsi -= v * delta / Lf * latency;     // epsi[t+1] = psi[t] - psides[t] - v[t] * delta[t] / Lf * dt

          v += a * latency;                     // v_[t+1] = v[t] + a[t] * dt
          // the speed must not exceed the reference value
          v = std::min(v, REFERENCE_V);

          // Defining the state vector, then 
          Eigen::VectorXd state(6);
          state << x, y, psi, v, cte, epsi;

          // Solving the curent problem
          auto vars = mpc.Solve(state, coeffs);

          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          double steer_value = vars[0] / deg2rad(25);
          double throttle_value = vars[1];

          // Due to the delayed nature of the actuator execution (and due to the car having a given
          // turning radius and an inertia) we have to take the actuator latency into consideration.
          // We practically instruct the vehicle to do a desired angle change more slowly, in smaller
          // steps in order to prevent oversteering and losing control over the car.
          // 2.0 as divider seems to work well for a latency of 100ms (ref. speed: 80 mph)
          // 4.0 as divider is still okay for a latency of 200ms (ref. speed: 80 mph)
          if ( ACTUATOR_LATENCY > 50 ) {
            steer_value /= 2.0 * (ACTUATOR_LATENCY / 100.0);
          }

          // Debug output
          std::cout << "    ==>  [delta, a]: " << steer_value << ", " << throttle_value << std::endl;

          json msgJson;
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          // Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          // Extracting trajectory (x,y) values from the result returned by the solver
          size_t const N = (vars.size() - 2) / 2;
          for( i = 0; i < N; i++ ) {

            mpc_x_vals.push_back(vars[2 + 2 * i]);
            mpc_y_vals.push_back(vars[2 + 2 * i + 1]);
          }

          // Add (x,y) points to list here, points are in reference to the vehicle's coordinate system.
          // The points in the simulator are connected by a green line
          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          // Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          // A simple solution would be to pass on the transformed waypoint coordinates
          // next_x_vals.resize(ptsx_car.size());
          // next_y_vals.resize(ptsx_car.size());
          // Eigen::VectorXd::Map(&next_x_vals[0], ptsx_car.size()) = ptsx_car;
          // Eigen::VectorXd::Map(&next_y_vals[0], ptsy_car.size()) = ptsy_car;

          // ...but a smoother line looks nicer, of course :)
          for ( x = 0.0; x < 91.0; x += 3.0 ) {
            next_x_vals.push_back(x);
            next_y_vals.push_back(polyeval(coeffs, x));
          }

          // Add (x,y) points to list here, points are in reference to the vehicle's coordinate system.
          // The points in the simulator are connected by a yellow line
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          // std::cout << msg << std::endl;

          // Latency
          // The purpose is to mimic real driving conditions where the car does actuate the commands instantly.
          // Feel free to play around with this value but should be to drive around the track with 100ms latency.
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(ACTUATOR_LATENCY));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }

      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {

    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());

    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;

  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }

  h.run();
}
