# gCLMSolverSuite
A collection of code for solving the generalized Constantin-Lax-Majda equation. Please read the introduction to my dissertation at http://www.pearcewashabaugh.com for an introduction to these equations and their importance.

Current usage:

- Open main.py and enter relevant initial data and mesh resolutions

- Run python main.py in a command line. Use Python 3. This will output the solutions to the Output folder.

- To plot the velocity field u, the vorticity w, or the Lagrangian trajectory eta, open plotter_u and adjust the input file to the corresponding variable in the output folder and make sure the correct number of time and space steps are used. Then run python plotter_u.py in a command line. 

- To plot the E-P trajectories, run python plotter_E-P.py, make sure to use the correct number of time steps.

