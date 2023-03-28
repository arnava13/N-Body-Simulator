#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computer Modelling Semester 2 Project Part 2

Requires Particle3D class from last semester.

Three functions to be used in simulations of gravitational systems. One function 
finds the gravitational force between two Particle3D instances. Next function 
computes the separations between all bodies in a set of Particle3D instances.
Final function calculates the gravitational forces between all bodies in a set
of Particle3D instances and the total potential of the system.

All units are astronomical units: Earth days, AUs, Earth Masses

Author: Arnav Agarwal
Number: S2163065
"""

import numpy as np
from particle3D import Particle3D

def gravitational_force(separation, m1, m2):
    """
    Computes the Newtonian gravitational force between two bodies.

    Parameters
    ----------
    separation : np.ndarray
        Separation vector in astronomical units of the two bodies' centres 
        of masses.
    m1, m2: float:
        Masses of each body in earth masses.

    Returns
    -------
    force : Force vector between bodies in units of M_earth * AU / Days**2

    """
    G = 8.887692593e-10 #space units (AUs, Earth Masses, Days)
    r = np.linalg.norm(separation) #Scalar separation of bodies.
    r_hat = separation/r #Unit vector in separating bodies
    force = G*m1*m2*r_hat/(r**2) #Force vector by Newton's law of gravitation.
    return force

def compute_separations(particles):
    """
    Generates an n x n x 3 array where each element represents 
    a cartesian component of the separation vector between two 
    particles in the input particle set.

    Parameters
    ----------
    particles : list
        List of Particle3D instances representing each particle.

    Returns
    -------
    separations : np.ndarray 
        n x n x 3 array representing the a cartesian component k
        of the separation vector between the ith particle and the 
        jth particle where n is the number of particles.
    """
    n = len(particles)
    separations = np.zeros([n, n, 3]) #Initialise n x n x 3 numpy array.
    for i in range(n): #Iterate through each particle pair.
        for j in range(n): 
            if i <= j: #Calculate separations for each element below the n x n diagonal.
                sep_vector = particles[j].position - particles[i].position
                separations[i,j] = sep_vector
            if i > j: #Use Newton's 3rd law for each element above the n x n diagonal.
                separations[i,j] = -separations[j,i]
    return separations


def compute_forces_potential(particles, separations):
    """
    Generates a matrix where each element represents a cartesian component 
    of the total force vector on a particle due to all other particles in the 
    input set.

    Parameters
    ----------
    particles : list
        List of Particle3D instances representing each particle.
    separations : np.ndarray
       n x n x 3 array representing the a cartesian component k
       of the separation vector between the ith particle and the 
       jth particle where n is the number of particles, as given by
       compute_separations.

    Returns
    -------
    forces : np.ndarray
        n x 3 array where each column represents a cartesian component 
        of the force vector on the particle in the each row.
    potential : float
        Total gravitational potential energy of the system in units of 
        M_earth * AU**2 / Days**2

    """
    G = 8.887692593e-10 #space units (AUs, Earth Masses, Days)
    n = len(particles)
    #Initialise n x n x 3 numpy array containing the 
    #Cartesian component k of the force vector between the ith particle
    #and the jth particle where n is the number of particles.
    force_on_i_due_to_j = np.zeros([n, n, 3])
    forces = np.zeros([n,3]) #Initialise n x 3 numpy array.
    potential = 0
    for i in range(n): #Iterate through each particle pair
        for j in range(n):
            if i < j:
                #Add potential due to pair to total potential.
                potential += G * particles[i].mass * particles[j].mass / np.linalg.norm(separations[i,j])
                #Compute each component below the n x n diagonal of force array using gravitational_force.
                force_on_i_due_to_j[i,j] = gravitational_force(separations[i,j], particles[i].mass, particles[j].mass)
            if i > j:
                #Compute each component above the n x n diagonal using Newton's 3rd law. 
                force_on_i_due_to_j[i,j] = -force_on_i_due_to_j[j,i]
        #Add up the ith row of vectors in forces matrix to find total force on ith particle. 
    forces = np.sum(force_on_i_due_to_j, axis=1) 
    return forces, potential
