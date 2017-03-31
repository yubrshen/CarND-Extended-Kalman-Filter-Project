# Extended Kalman Filter Project
Self-Driving Car Engineer Nanodegree Program

---
## Introduction

This is my practice of building extended Kalman filter with the help of Udacity's starter code. 

The output of the RMSE met the project rubric. Here are the RMSE:

- ./data/sample-laser-radar-measurement-data-1.txt
Accuracy - RMSE:
0.0651651
0.0605925
 0.531474
 0.544462

- ./data/sample-laser-radar-measurement-data-2.txt
Accuracy - RMSE:
0.186482
0.191275
0.502008
0.807094

I found that the key is to have measurements good enough for Kalman filter operations. 
Converting the internal state in Cartesian coordinates to Polar coordinates together with 
the rate of change in the distance $\rho$ also contributes to improving RMSE. 

Caution is applied to save computation, and also keep algorithm consistent to avoid implementation
error. 

The code might be further refined with simpler module architecture with functional programming flavor. 
The class hierarchy (FusionEKE/KalmanFilter/Tools) may be an over-kill for the straight forward algorithm implementation. 

The following shows the building and running instructions. 

## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make` 
   * On windows, you may need to run: `cmake .. -G "Unix Makefiles" && make`
4. Run it: `./ExtendedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./ExtendedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt`

