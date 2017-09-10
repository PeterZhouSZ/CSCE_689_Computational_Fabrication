Name:
	- Ying Wang

Citations:
	- C++ code for Möller–Trumbore intersection algorithm.
	- https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm

Notes:
	- Added the 5th argument(Input 0: one ray; Input 1: 6 rays).
	- Added some parallel code that might need OpenMP library.(main.cpp: line 7, line 193).
		- Install OpenMP:
			brew install gcc --without-multilib
			brew install openmpi

	- Also changed the Cmakelist. (g++-7 xxx.c -o xxx -fopenmp)

