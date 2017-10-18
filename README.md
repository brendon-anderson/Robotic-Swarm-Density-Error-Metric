# Robotic-Swarm-Density-Error-Metric

This repository contains the code necessary to generate two benchmarks for quantitative assessment of robotic swarm coverage over a prescribed target distribution. These benchmarks make use of an error metric defined as the L1 norm of the difference between the target distribution and the swarm blob function, where the swarm blob function is the normalized superposition of Gaussian blobs at each robot location. To use these MATLAB scripts, the user must input various quantities, such as swarm parameters, numerical computation parameters, and target distribution parameters and definitions.

The first benchmark for comparison is the probability density function of the error metric. The script provided utilizes Monte Carlo integration to generate the cumulative density function of the error metric for a given swarm size and effective robot radius. The CDF is then curve fit, then the curve fit is differentiated to generate the probability density function of interest.

The second benchmark for comparison is the optimal error metric value for a prescribed swarm size and robot radius. This is computed using MATLAB's fmincon optimization function. Since fmincon provides a local minimum upon convergence, the script is written to allow for looping of the solution to implement multistart optimization.

A third script is also provided in which the robot radius is included in the optimization. This allows for a method of optimal swarm design by assisting the designer in finding the ideal combination of swarm size and robot radius to achieve a desired coverage over the target distribution.
