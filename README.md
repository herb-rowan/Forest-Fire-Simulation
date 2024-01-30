Forest Fire Simulation
Overview
This project implements a forest fire simulation using a parallelized approach with OpenMP. The forest fire problem is an example of a celluar automata problem. 
The simulation models the spread of fire in a grid-based forest environment, considering factors such as wind direction.
DynamicSceduling.cpp is the main code, with no wind direction implemented. 
wind.cpp is the version of code with wind implemented.
I've also attached an MP4 file where you can see the results animated, I did all the animation in python.
I also attached a report where you can find a more in depth anaylsis of my results.

You can run the code using:
g++ -o forest_fire_omp forest_fire_omp.cpp -fopenmp



