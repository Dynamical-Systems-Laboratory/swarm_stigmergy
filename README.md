# Swarm stigmergy repository

The code in these folders allows to coordinate a swarm of robots through stigmergy, by releasing "traces" in the environment. All codes are in Matlab.
All computations for coordination are made at a continuum scale, such that we require the density of the robots and we compute the necessary density of the traces.
There are two main folders:
1. 1D
   These are all simulations in 1D, over a circle.
   1.1. Continuum
       These are all simulations at a continuum level only. We have three subfolders:
       1.1.1. Static
              This simulation requires a static density of the robots, which does not vary over time.
              The main file to run the simulation is "main.m".
              The other files are necessary function, each well detailed in their initial comment.
       1.1.2. Dynamic_fixed
              This simulation requires a time-varying density of the robots, which simulates a wave in a non-dispersive medium.
              The main file to run the simulation is "main.m".
              The other files are necessary function, each well detailed in their initial comment.
       1.1.3. Dynamic_modulating_profile
              This simulation requires a time-varying density of the robots, which simulates a wave in a dispersive medium.
              The main file to run the simulation is "main.m".
              The other files are necessary function, each well detailed in their initial comment.
   1.2. Discrete
       These are all simulations at a discrete level, with discrete robots. We continuify the density and then discretize the density of the traces to be deployed. We          have three subfolders:
       1.2.1. Static
              This simulation requires a static density of the robots, which does not vary over time.
              The main file to run the simulation is "main.m".
              The other files are necessary function, each well detailed in their initial comment.
       1.2.2. Dynamic_fixed
              This simulation requires a time-varying density of the robots, which simulates a wave in a non-dispersive medium.
              The main file to run the simulation is "main.m".
              The other files are necessary function, each well detailed in their initial comment.
       1.2.3. Dynamic_modulating_profile
              This simulation requires a time-varying density of the robots, which simulates a wave in a dispersive medium.
              The main file to run the simulation is "main.m".
              The other files are necessary function, each well detailed in their initial comment.
2. 2D
   These are all simulations in 2D, over a periodic rectangle.
   The simulation requires the swarm to track the trajectory of Leonardo da Vinci's lion automaton.
   The main file to run the simulation is "main.m".
   The other files are necessary function, each well detailed in their initial comment.
   2.1 q_T_compute
       This folder generates the file "qTs.mat" to be used in the main simulation folder. This is a discretized version of the kernel.
       You will need to run this script if you change the features of the mesh in the main simulation.
       The main file to run the simulation is "main_q_T.m".
       The other file in the folder is a necessary function, well detailed in its initial comment.
