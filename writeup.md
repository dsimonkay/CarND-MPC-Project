# Model Predictive Control Project

The goal of this project was to build an MPC controller that can be used to steer and drive the the vehicle in the Term 2 simulator around the lake.

The simulator delivers information over a local websocket, providing the following values to the program:
 * global waypoint positions \[(x<sub>1</sub>, y<sub>1</sub>), (x<sub>2</sub>, y<sub>2</sub>), ... (x<sub>n</sub>, y<sub>n</sub>)\],
 * orientation of the vehicle (in radians),
 * global position of the vehicle (x, y),
 * current steering angle (in radians),
 * current throttle value (\[-1, 1\]),
 * current speed (in MPH).

The simulator expects a steering angle and a throttle value. Both values have to be normalized into the \[-1, 1\] range.

## [Rubric](https://review.udacity.com/#!/rubrics/896/view) Points:


### Compilation

#### Your code should compile
The code compiles without errors or warnings. I've applied some changes to the Eigen library source code in order to get rid of the warnings which I encountered during compilation. I made changes to the following files:
 * [src/Eigen-3.3/Eigen/src/Core/AssignEvaluator.h](https://bitbucket.org/eigen/eigen/commits/b4f969795d1b0adbb43ebdd8c6bbcb42cb559687?at=3.3)
 * [src/Eigen-3.3/Eigen/src/Core/products/GeneralMatrixVector.h](https://bitbucket.org/eigen/eigen/commits/131da2cbc6958c6576d23872946aac2da5a8467c?at=3.3)

Source: [1402 - Warning: enum constant in boolean context](http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1402#c3)


### Implementation

#### The Model
I used the simple kinematic vehicle model that was introduced in the classroom. It uses a state vector by which we can keep track of the vehicle's current state and which enables us to make predictions about the vehicle's future states using actuated inputs.

 The model consists of the following elements:
 
 * the x coordinate of the vehicle (*x*),
 * the y coordinate of the vehicle (*y*),
 * the orientation of the vehicle (&psi;),
 * the velocity of the vehicle (*v*)
 * the current cross track error (*cte*),
 * the current orientation error (*e&psi;*),
 * the control input value for the steering angle (&delta;),
 * the control input value for acceleration/deceleration (*a*).

The state vector consists of the first six elements of the model components listed above: \[ x, y, &psi;, v, cte, e&psi; \]

If the state vector and the control inputs are known, we can calculate the values of the new state variables in *dt* time using the following equations:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;x<sub>t+1</sub> = x<sub>t</sub> + v<sub>t</sub> \* cos(&psi;<sub>t</sub>)  \* dt<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;y<sub>t+1</sub> = y<sub>t</sub> + v<sub>t</sub> \* sin(&psi;<sub>t</sub>)  \* dt<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&psi;<sub>t+1</sub> = &psi;<sub>t</sub> + v<sub>t</sub> / Lf \*   &delta;<sub>t</sub> \* dt<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;v<sub>t+1</sub> = v<sub>t</sub> + a<sub>t</sub> \* dt<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cte<sub>t+1</sub> = cte<sub>t</sub> + v<sub>t</sub> \* sin(e&psi;<sub>t</sub>) \* dt<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;e&psi;<sub>t+t</sub> = e&psi;<sub>t</sub> + v<sub>t</sub> / Lf  \* &delta;<sub>t</sub> \* dt<br/>

cte<sub>t</sub> is the difference between the itinary and the current vehicle position. If we have a polynomial *f* which describes the refrence line, then the cross track error at time *t* is: cte<sub>t</sub> = f(x<sub>t</sub>) - y<sub>t</sub>

e&psi;<sub>t</sub> is the difference between the current and the desired orientation of the vehicle: e&psi;<sub>t</sub> = &psi;<sub>t</sub> - &psi;des<sub>t</sub>
&psi;des<sub>t</sub> can be calculated as the tangential angle of the above mentioned polynomial *f* at *x<sub>t</sub>*, or arctan(f'(x<sub>t</sub>)), so e&psi;<sub>t</sub> = &psi;<sub>t</sub> - arctan(f'(x<sub>t</sub>))


#### Timestep Length and Elapsed Duration (N & dt)
At first I was brave enough to set the prediction horizon to 2.0 seconds and tried these *(N, dt)* values:
 * 40 &ndash; 0.05 &mdash; did not work: the vecicle drove off-track before the first curve
 * 20 &ndash; 0.1 &mdash; worked, but was not stable enough
 * 10 &ndash; 0.2 &mdash; didn't work

Then I set the prediction horizon to 1.5 seconds and tried these *(N, dt)* values:
 * 20 &ndash; 0.075 &mdash; yielded un unstable behaviour
 * 10 &ndash; 0.15 &mdash; it did work, although the vehicle slowed down considerably in the curves
 * 5 &ndash; 0.3 &mdash; the car left the track right after start
 
 After that I kept on decreasing the prediction horizon until I got a stable behaviour around 1 second. A few tries:
 * 12 &ndash; 0.08 &mdash; erratic trajectory lines in the second and third curves (although driving in the curves faster than in the "stable" solution)
 * 11 &ndash; 0.09 &mdash; erratic trajectory lines in the second and third curves (and also slowing down in the curves)
 * 9 &ndash; 0.1 &mdash; stable solution
 * 10 &ndash; 0.1 &mdash; stable solution
 * 11 &ndash; 0.1 &mdash; stable solution
 * 12 &ndash; 0.1 &mdash; pretty stable solution
 
 At the end picked a stable solution, choosing 10 for N and 0.1 for dt.


#### Polynomial Fitting and MPC Preprocessing
Before fitting a polynomial to the waypoint coordinates, I transformed these into the vehicle's coordinate system using the homogeneous transformation. Using the vehicle's coordinate system sets *x*, *y* and *&psi;* to zero in the state vector (as the car is at the origo of its coordinate system, heading towards the **x** axis). It also greatly simplifies the further calculation of *cte* and *e&psi;*.


#### Model Predictive Control with Latency
The idea behind incorporating actuator latency into the system is simple: after having calculated the initial (=current) state of the car, we virtually let the vehicle keep on driving for the predefined amount of time (=latency) and update the state vector by the changes made to its components during the elapsed time. We then run the MPC solver with the updated state vector as the new initial state. As a result we receive exactly the control input values that have to be executed after the predefined time delay in order to keep the vehicle on the desired track.


### Simulation

#### The vehicle must successfully drive a lap around the track
The car can successfully drive around the lake with a reference speed of 80 MPH.