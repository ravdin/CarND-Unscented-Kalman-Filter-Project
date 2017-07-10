# Extended Kalman Filter Project
Self-Driving Car Engineer Nanodegree Program

[image1]: ./screenshot.png "Screenshot"

This project implements an Unscented Kalman Filter (UKF) in C++ to predict and track the position and velocity of a moving object.  The goal is to make use of both lidar and radar data.  The algorithm makes a prediction based on previous measurements and updates the known state.  The predictions become more accurate over time as measured by the root mean squared error (RMSE).

The RSME for the given dataset is [0.07, 0.08, 0.34, 0.22], within the acceptable constraints.  See the screenshot for the final result.

![alt text][image1]
