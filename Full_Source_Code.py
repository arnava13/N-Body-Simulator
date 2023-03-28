#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computer Modelling Semester 2 Project Part 2

Takes an input planetary data file an runs a verlet integration regime to model the positions of the bodies and the energy of the system
over a series of timesteps, in response to the gravitational interactions between all bodies in the system.

Gives all positions in relative to system centre of mass and does not include COM kinetic energy in energy outputs. Conducts all calculations in this frame.

Makes use of the Particle3D class created in Semester 1.

Requires a "Basic_Functions" file containing functions to calculate gravitational force between two bodies, compute the vectorial separations of
all bodies and to compute the forces on all bodies & the total potential of the system.

Writes an output trajectory file in xyz format. Writes an output energies file with the format specified in a comment line.
Returns the calculated aphelion and perhelion of each body relative to the sun to the command line, if a sun is supplied. Also returns orbital periods around the sun if supplied.
Returns the calculated perigee and apogee to the command line if the earth and moon are present, and also the orbital period of the moon around the earth.

Generates a plot of total energy in the system over time, in order to demonstrate conservation of energy.

All units are astronomical units: Earth days, AUs, Earth Masses.

Author: Arnav Agarwal
Number: S2163065.
"""

from particle3D import Particle3D
import numpy as np
import Basic_Functions as bf
import sys
import scipy as scipy
import matplotlib.pyplot as plt

def read_initial_conditions(input_filename):
    """
    Parameters
    ----------
    input_filename : string
        Name of input file in standard XYZ format to read masses, positions, velocities for simulation particles.

    Returns
    -------
    particles: 1D numpy array
        Array containing Particle3D instances for each simulation body.

    """
    #Initialise empty list of particles.
    particles = []
    
    #open input file for reading
    with open(input_filename, "r") as inputfile:
        
        #read lines in input file and iterate through each line
        particlelines = inputfile.readlines()
    
        for line in particlelines:
            #Split each line and extract particle parameters.
            parameters = line.split()
            label = parameters[0]
            mass = parameters[1]
            position = np.array([float(parameters[2]),float(parameters[3]),float(parameters[4])])
            velocity = np.array([float(parameters[5]), float(parameters[6]), float(parameters[7])])
            #Initialise Particle3D instance using parameters extracted.
            particle = Particle3D(label, mass, position, velocity)
            #Append Particle3D instances 
            particles.append(particle)
            
    #Convert particles list to numpy array.
    particles = np.array(particles)
    
    return particles
             
def compute_periods(particles, trajectories, dt):
    """

    Parameters
    ----------
    particles : 1D numpy array
        Array of particle3D instances representing bodies to be simulated.
    trajectories : 3D Numpy array
       Array of positions of each simulated body at each timestep [time, particle, pos[xyz]].
    dt : float
        Simulation timestep used.

    Returns
    -------
    perhelia_dict : dictionary
        Dictionary containing periods of each simulated body around sun (body:period). None if sun not present.

    """
    #Enumerate number of particles, number of timsteps
    numparticles = len(particles)
    numstep = np.size(trajectories, axis = 0)
    
    #Initialise variable to store index of sun in particles array.
    sunindex = None

    #find index of sun particles array:
    for p in range(numparticles):
        if particles[p].label == "Sun":
                sunindex = p
    
    #Return None if Sun not in particle list.
    if sunindex == None:
        periods_dict = None
        
    else:
        #Generate empty dictionary of periods with body names as keys.
        periods_dict = {particle.label : None for particle in particles}
        #Remove sun from periods dictionary.
        del periods_dict["Sun"]
        
        #iterate through particles which are not the sun
        for p in range(numparticles):
            if p != sunindex:
                #Find initial position vector of planet from sun to measure angle from
                init_pos = particles[p].position - particles[sunindex].position
                #Create 1D array to store particle angle at each timestep.
                particle_angles = np.zeros(numstep)
                
                for i in range(numstep):
                    #Calculate new position vector from sun at each timestep.
                    new_pos = trajectories[i, p] - trajectories[i, sunindex]
                    
                    #Use dot product formula to find angle from initial position vector.
                    dotproduct = np.dot(init_pos, new_pos)
                    magsproduct = np.linalg.norm(init_pos) * np.linalg.norm(new_pos)
                    angle = np.arccos(min(dotproduct/magsproduct,1))
                    #Append angle to particle angles array.
                    particle_angles[i] = angle
                    
                #Find timestep number of minima using Scipy find peaks on inverted particle angles array.
                angle_minima_locs = scipy.signal.find_peaks(-particle_angles)[0]
                #Initialise list to successive periods between angle minima for each particle.
                periods = []
                #iterate through entries in array of angle minima timestep numbers.
                
                for i in range(len(angle_minima_locs) - 1):
                    #Calculate period as number of timesteps between each peak multiplied by timestep dt.
                    period = (angle_minima_locs[i+1] - angle_minima_locs[i]) * dt
                    #Append entry to period list.
                    periods.append(period)
                    
                #If insufficient periods in period list, set average period of particle as None.
                if np.size(periods) < 1:
                    avg_period = None
                    
                else:
                    #Find average period for each particle
                    periods = np.array(periods)
                    avg_period = np.average(periods)
                    
                #Append average period of each particle to periods dictionary.
                periods_dict[particles[p].label] = avg_period
    
    return periods_dict

def compute_moon_period(particles, trajectories, dt):
    """

    Parameters
    ----------
    particles : 1D numpy array
        Array of particle3D instances representing bodies being simulated.
    trajectories : 3D Numpy array
       Array of positions of each simulated body at each timestep [time, particle, pos[xyz]].
    dt : float
        Simulation timestep used.

    Returns
    -------
    moon_period : float
        

    """
    #Enumerate number of particles, number of timsteps
    numparticles = len(particles)
    numstep = np.size(trajectories, axis = 0)
    
    #Initialise variable to store index of earth and moon in particles array.
    moonindex = None
    earthindex = None

    #find index of sun particles array:
    for p in range(numparticles):
        if particles[p].label == "Moon":
                moonindex = p
        if particles[p].label == "Earth":
                earthindex = p
    
    #Return None if moon or earth not in particle list.
    if moonindex == None and earthindex != None:
        moon_period = "No Moon"
    elif moonindex != None and earthindex == None:
        moon_period = "No Earth"
    elif moonindex == None and earthindex == None:
        moon_period = "No Moon Earth"
        
    else:
        #Find initial position vector of moon from earth to measure angle from
        init_pos = particles[moonindex].position - particles[earthindex].position
        
        #Create 1D array to store moon angle at each timestep.
        particle_angles = np.zeros(numstep)
        
        for i in range(numstep):
            #Calculate new position vector from sun at each timestep.
            new_pos = trajectories[i, moonindex] - trajectories[i, earthindex]
            #Use dot product formula to find angle from initial position vector.
            dotproduct = np.dot(init_pos, new_pos)
            magsproduct = np.linalg.norm(init_pos) * np.linalg.norm(new_pos)
            angle = np.arccos(min(dotproduct/magsproduct,1))
            #Append angle to particle angles array.
            particle_angles[i] = angle
            
        #Find timestep number of minima using Scipy find peaks on inverted particle angles array.
        angle_minima_locs = scipy.signal.find_peaks(-particle_angles)[0]
        
        #Initialise list to successive periods between angle minima for each particle.
        periods = []
        
        #iterate through entries in array of angle minima timestep numbers.
        for i in range(len(angle_minima_locs) - 1):
            #Calculate period as number of timesteps between each peak multiplied by timestep dt.
            period = (angle_minima_locs[i+1] - angle_minima_locs[i]) * dt
            #Append entry to period list.
            periods.append(period)
            
        #If insufficient periods in period list, set moon period as None.
        if np.size(periods) < 1:
            moon_period = "Insuff Period"
            
        else:
            #Find average period.
            periods = np.array(periods)
            moon_period = np.average(periods)
            
    return moon_period

def compute_energy_deviations(potentials, kinetic_energies):
    """
    Generates a plot of total energy of a simulated planetary system over time.
    Returns an array containing the total energy at each timestep.
    
    Parameters
    ----------
    potentials : 1D Numpy Array
        Array storing the total potential energy of the simulated system at each timestep.
    kinetic_energies : 1D Numpy Array
        Array storing the total kinetic energy of the simulated system at each timestep.

    Returns
    -------
    total_energies : 1D Numpy Array
        Array storing the total energy of the simulated system at each timestep..

    """
    #Enumerate number of timesteps.
    numstep = len(potentials)
    
    #Initialise array to store total energy at each timestep.
    total_energies = np.zeros(numstep)
    
    #Iterate through timesteps and calculate total energy at each time.
    for i in range(numstep):
        total_energies[i] = potentials[i] + kinetic_energies[i]
    
    #Generate plot of total energy over time.
    plt.figure()
    plt.title('Total Energy of Planetary System')
    plt.xlabel('Time (Earth Days)')
    plt.ylabel('Total Energy ($GM_☉$)')
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    plt.plot(range(numstep), total_energies, label = "Total Energy")
    plt.show()
    
    return(total_energies)
    
def compute_perhelia_aphelia_perigee_apogee(particles, trajectories):
    """
    Parameters
    ----------
    particles : 1D numpy array
        Array of particle3D instances representing bodies being simulated.
    trajectories : 3D Numpy array
        Array of positions of each simulated body at each timestep [time, particle, pos].

    Returns
    -------
    perhelia_dict : dictionary
        Dictionary containing perhelia of each simulated body (planet:perhelion). None if sun not present.
    aphelia_dict : dictionary
        Dictionary containing aphelia of each simulated body (planet:aphelion). None if earth or moon not present.
    perigee : float
        Perigee of earth and moon, None if one or both not present.
    apogee : float
        apogee of earth and moon, None if one or both not present.

    """
    #Enumerate number of particles and timesteps.
    numparticles = len(particles)
    numstep = np.size(trajectories, axis = 0)
    
    #initialise variables to store index of moon, sun and earth in particles array.
    moonindex = None
    sunindex = None
    earthindex = None
    
    #find index of sun and moon in particles array:
    for p in range(numparticles):
        if particles[p].label == "Sun":
                sunindex = p
        if particles[p].label == "Moon":
                moonindex = p
        if particles[p].label == "Earth":
                earthindex = p
                
    #Return None for output perhelia and aphelia dictionaries if no sun present.
    if sunindex == None:
        perhelia_dict = None
        aphelia_dict = None
        
    else:
        #Generate empty dictionaries of perhelia and aphelia with only body names as keys.
        perhelia_dict = {particle.label : None for particle in particles}
        aphelia_dict = {particle.label : None for particle in particles}
        #Delete sun as entries in both dictionaries.
        del perhelia_dict["Sun"]
        del aphelia_dict["Sun"]
        
        #Iterate through particles which are not the sun.
        for p in range(numparticles):
            if p != sunindex:
                #Set aphelion variable as zero and initiate perhelion variable to test inequalities later.
                aphelion = 0
                perhelion = float('inf')
                
                #Iterate through timesteps.
                for i in range(numstep):
                    #Calculate distance of each particle from sun at each timestep.
                    calc_dist = np.linalg.norm(trajectories[i, p] - trajectories[i, sunindex])
                    #Set initial value of perhelion variable to first timestep distance.
                    #If new distance is the larger than previous largest, or smaller than previous smallest, set as new aphelion/perhelion, respectively.
                    if calc_dist > aphelion:
                        aphelion = calc_dist
                    if calc_dist < perhelion:
                        perhelion = calc_dist
                        
                #Set perhelia and aphelia for each particle in output dictionaries.
                perhelia_dict[particles[p].label] = perhelion
                aphelia_dict[particles[p].label] = aphelion
                
    #If moon or earth not present in input particle array, return perigee and apogee as None.
    if moonindex == None and earthindex != None:
        apogee = "No Moon"
        perigee = None
    if moonindex != None and earthindex == None:
        apogee = "No Earth"
        perigee = None
    if moonindex == None and earthindex == None:
        apogee = "No Moon Earth"
        perigee = None
        
    elif moonindex != None and earthindex != None:
        #Set apogee as zero and perigee as 'inf' to test inequalities later.
        apogee = 0
        perigee = float("inf")
        
        #Iterate through timesteps
        for i in range(numstep):
            #Calculate distance between sun and moon at each timestep.
            calc_dist = np.linalg.norm(trajectories[i, moonindex] - trajectories[i, earthindex])
            #If new distance is the larger than previous largest, or smaller than previous smallest, set as new apogee/perigee, respectively.
            if calc_dist > apogee:
                apogee = calc_dist
            if calc_dist < perigee:
                perigee = calc_dist
                
    return perhelia_dict, aphelia_dict, perigee, apogee

def main():
    #Read input parameters from command line.
    # Here we expect three things:
    #    the name of this file
    #    the name of the input planetary data file
    #    the timestep parameter to be used for the simulation
    #    total number of timesteps to run simulation for
    #    filename of output planet trajectories file
    #    filename of output energies file
    # So we start by checking that all three are specified and quit if not,
    # after giving the user a helpful error message.
    if len(sys.argv) !=6 :
        print("\nYou did not provide the correct number of arguments, the correct syntax is:\n")
        print("\n{full_source_code.py} {input planetary data file} {simulation timestep} {number of timesteps} {output energy file} {output xyz file}\n")
        sys.exit(1)

    #Initialise variables for input filename, timestep (dt), number of timesteps (numstep), and filenames of output energies and trajectories files.
    input_filename = sys.argv[1]
    dt = float(sys.argv[2])
    numstep = int(sys.argv[3])
    out_energies_filename = sys.argv[4]
    out_trajectories_filename = sys.argv[5]

    #Read initial conditions from input file.
    particles = read_initial_conditions(input_filename)

    #Open output files in write mode.
    out_energies_file = open(out_energies_filename, "w")
    out_trajectories_file = open(out_trajectories_filename, "w")
            
    #Calculate and subtract COM velocity from each particle's velocity.
    total_momentum = sum([particle.velocity * particle.mass for particle in particles])
    total_mass = sum([particle.mass for particle in particles])
    com_vel = total_momentum/total_mass
    for i in range(len(particles)):
        particles[i].velocity = particles[i].velocity - com_vel
    
    #Calculate initial separations, forces and potentials
    separations = bf.compute_separations(particles)
    forces, potentials = bf.compute_forces_potential(particles, separations)
    
    #Initialise integer to quantify point numbers.
    pointnum = 0
    
    #Initialise arrays storing particle positions, total KE and system potential at each timestep.
    numparticles = int(len(particles))
    trajectories = np.zeros([numstep, numparticles, 3])
    potentials = np.zeros(numstep)
    kinetic_energies = np.zeros([numstep])
    
    #Create format strings for each output file entry.
    xyzentryheadformat = str(numparticles) + "\nPoint = {point}\n"
    energyentryformat = "{time},{pot:.9f},{kin:.9f},{tot:.9f}\n"
    
    
    #Write comment line on output energies file explaining format.
    out_energies_file.write("Time in days,Potential Energy,Kinetic Energy,Total Energy\n")
    
    #Run Verlet Integration.
    
    #Iterate through times in range of numstep.
    for i in range(numstep):
        
        #Increase point number for writing.
        pointnum += 1
        
        #Write entry headers for each output file.
        out_trajectories_file.write(xyzentryheadformat.format(point = pointnum))
        
        #Initialise variable to store total kinetic energy at each timestep.
        kinetic_energy = 0
        
        #Iterate through each particle to find new position.
        for p in range(numparticles):
            #Update position to 2nd order using current force.
            particles[p].update_position_2nd(forces[p], dt)
            #Write new positions to output trajectory file.
            out_trajectories_file.write(str(particles[p]) + "\n")
            #Write new positions to trajectories array.
            trajectories[i, p] = particles[p].position
            
        
        #Calculate new separations, forces and potential and new particle positions.
        newseparations = bf.compute_separations(particles)
        newforces, potential = bf.compute_forces_potential(particles, newseparations)
        
        #Write new potential to potentials array, noting sign convention.
        potentials[i] = -potential
        
        #Update velocity of each particle at new position.
        for p in range(len(particles)):
            particles[p].update_velocity(0.5*(forces[p] + newforces[p]), dt)
            #Find total KE at new point.
            kinetic_energy += particles[p].kinetic_energy()
        
        #Write new kinetic energy to kinetic energies array.
        kinetic_energies[i] = kinetic_energy
        
        #Write current time and energies to output energy file.
        total_energy = kinetic_energy - potential #Noting sign convention.
        time = dt * i 
        out_energies_file.write(energyentryformat.format(time = time, pot = potential, kin = kinetic_energy, tot = total_energy)) 
        
        #Update forces array for next iteration.
        forces = newforces
    
    #Close output energies and trajectories file.
    out_energies_file.close()
    out_trajectories_file.close()
    
    #Compute perhelia, aphelia, perigee and apogee from function
    perhelia_dict, aphelia_dict, perigee, apogee = compute_perhelia_aphelia_perigee_apogee(particles, trajectories)
    
    #Print results of perhelia, aphelia, perigee and apogee for each planet.
    if apogee == "No Moon":
        print("\nThere is no \"Moon\" in the collection of bodies supplied, so the perigee and apogee cannot be found.\n")
    elif apogee == "No Earth":
        print("\nThere is no \"Earth\" in the collection of bodies supplied, so the perigee and apogee cannot be found.\n")
    elif apogee == "No Moon Earth":
        print("\nThere is no \"Moon\" and no \"Earth\" in the collection of bodies supplied, so the perigee and apogee cannot be found.\n")
    else: 
        print("\nThe perigee of the Earth's moon is " + str(round(perigee, 3)) + " AU." )
        print("The apogee of the Earth's moon is " + str(round(apogee, 3)) + " AU.\n" )
        
    if aphelia_dict == None or perhelia_dict == None:
        print("There is no \"Sun\" in the collection of bodies supplied, so perhelia and aphelia cannot be found.\n")
    else: 
        for particle in particles:
            if particle.label != "Sun":
                print("The perhelion of " + particle.label + " is " + str(round(perhelia_dict[particle.label], 3)) + " AU." )
                print("The aphelion of " + particle.label + " is " + str(round(aphelia_dict[particle.label], 3)) + " AU.\n" )
    
    
    #Use compute_periods to compute dictionary storing periods of each particle.
    periods_dict = compute_periods(particles, trajectories, dt)
    
    #Print orbital periods around Sun:
    if periods_dict == None:
        print("There is no \"Sun\" in the collection of bodies supplied, so orbital periods around the Sun cannot be found.\n")
    else:
        for particle in particles:
            if particle.label != "Sun":
                if periods_dict[particle.label] == None:
                    print("The simulation time was too short to find the orbital period of " + particle.label + " around the Sun.\n")
                else:
                    print("The orbital period of " + particle.label + " around the sun is " + str(round(periods_dict[particle.label], 3)) + " days.\n")
                    
    #Use compute_moon_period to find period of moon around earth.
    moon_period = compute_moon_period(particles, trajectories, dt)
    
    #Print period of Moon around Earth:
    if moon_period == "No Moon":
        print("There is no \"Moon\" in the collection of bodies supplied, so orbital period of the Moon around the Earth cannot be found.\n")
    elif moon_period == "No Earth":
        print("There is no \"Earth\" in the collection of bodies supplied, so orbital period of the Moon around the Earth cannot be found.\n")
    elif moon_period == "No Moon Earth":
        print("There is no \"Earth\" and no \"Moon\" in the collection of bodies supplied, so orbital period of the Moon around the Earth cannot be found.\n")
    elif moon_period == "Insuff Period":
        print("The simulation time was too short to find the orbital period of the moon around the earth.\n")
    else:
        print("The orbital period of the Moon around the Earth is " + str(round(moon_period, 3)) + " days.\n")
                    
    
    #Plot total energy of system over time using compute_energy_deviations.
    total_energies = compute_energy_deviations(potentials, kinetic_energies)
        
    
if __name__ == "__main__":
    main()