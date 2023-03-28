Gravitational System Simulation Operating Instructions

------------------------------------------------------------------------------------------

Run using:

full_source_code.py input planetary data file} {simulation timestep} {number of timesteps} {output energy file} {output xyz file}

------------------------------------------------------------------------------------------
Input file:

- .txt format.
- Mass in units of earth masses (M☉).
- Positions in units of astronomical units (AU)
- Velocity in units of AU / earth day. 
- Initial conditions of simulated bodies given in input file with following format:

{Planet 1} {mass} {x} {y} {z} {x-velocity} {y-velocity} {z-velocity}
{Planet 2} {mass} {x} {y} {z} {x-velocity} {y-velocity} {z-velocity}
.
.
.
------------------------------------------------------------------------------------------
Output files:

////// Output xyz file: 
- .xyz format.
- Positions given in astronomical units (AU). 
- Rounded to 7 decimal places (approx. 10km).
- Gives the position of all simulation particles at each timestep in standard xyz format:

{Number of particles}
Point = {timestep number}
{Planet 1} {x} {y} {z}
{Planet 2} {x} {y} {z}
.
.
.

////// Output energy file:
-.csv format.
- Energies given in units of the solar gravitational parameter (Earth masses * AU^2 / earth Days^2) (GM☉).
- Rounded to 9 decimal places (to show what little variation in total energy exists).
- Gives the potential, kinetic and total energy of simulated system at each timestep with the following line format:

{Time in days},{Potential Energy},{Kinetic Energy},{Total Energy}
.
.
.

-----------------------------------------------------------------------------------------
Output plot:

- Generates a matplotlib plot showing variation of total energy with time.

-----------------------------------------------------------------------------------------
Representative files:

- Representative files have been generated using solar_system.txt as the input, with a timestep of 1 day for a total of 271575 days. This is long enough to get three orbits of Pluto, the longest period in the solar system. Shorter timesteps result in an incorrect orbital period for the Moon around the Earth.