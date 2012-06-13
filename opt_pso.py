#!/usr/bin/env python
# Particle Swarm Optimization Toolkit
# Fredrik Boulund
# Originally MATLAB 2010,
# Python Object Oriented implementation 2011


from sys import exit
from random import random, uniform

__doc__ = """
    PSO module
    Fredrik Boulund, 2011
    Version 0.43
    2011-03-24
    """

class Particle:
    """
    A simple particle.

    It has randomly initialized starting position in three dimensions
    and random starting velocity according to maximumVelocity.
    """
    
    def __init__(self, variableRange, maximumVelocity):
        # Set particle random starting position within valid variable range
        # Set random starting velocity within (0,maximumVelocity)
        self.pos = []
        self.velocity = []
        if len(maximumVelocity) == 1:
            for dim in variableRange:
                self.pos.append(uniform(dim[0],dim[1]))
                self.velocity.append(uniform(0,maximumVelocity[index]))
        else:
            for dim in xrange(0,len(variableRange)):
                self.pos.append(uniform(variableRange[dim][0],variableRange[dim][1]))
                self.velocity.append(uniform(-maximumVelocity[dim],maximumVelocity[dim]))
        self.bestPos = self.pos
        self.bestFitness = None

    def __str__(self):
        return ''.join(["pos=",str(self.pos),"\nv=",str(self.velocity),"\nbestPos=",str(self.bestPos),"\nbestFitness=",str(self.bestFitness)])


    def setBestPos(self, pos):
        """
        Sets the best position for current particle.
        """
        self.bestPos = pos

    def setBestFitness(self, fitness):
        """
        Sets the best fitness value for current particle.
        """
        self.bestFitness = fitness


    def updatePosition(self,variableRange):
        """
        Updates the position of a particle.

        Sums the current position with the current velocity to compute the new position
        """

        newPos = []
        for dim in xrange(0,len(self.pos)):
            newPos.append(self.pos[dim]+self.velocity[dim])
            if newPos[dim] < variableRange[dim][0]:
                newPos[dim] = variableRange[dim][0]
            elif newPos[dim] > variableRange[dim][1]:
                newPos[dim] = variableRange[dim][1]
        self.pos = newPos
        

    def updateVelocity(self, globalBest, maximumVelocity, c1, c2):
        """
        Updates the velocity of current particle.

        Takes into account the particle best and global best positions found
        in the swarm. It also enforces a restriction on the velocity composant
        in any direction.
        """

        currentPos = self.pos
        currentVelocity = self.velocity
        bestPos = self.bestPos

        # Account for the global and particle best position, ELEMENTWISE!
        newVelocity = []
        for dim in xrange(0,len(currentVelocity)):
            newVelocity.append(currentVelocity[dim] + 
                               c1*random()*(bestPos[dim] - currentPos[dim]) + 
                               c2*random()*(globalBest[dim] - currentPos[dim]))


        # Impose restrictions on maximum velocity and cap
        if len(maximumVelocity) > 1:
            for dim in xrange(0,len(newVelocity)):
                if newVelocity[dim] > maximumVelocity[dim]:
                    newVelocity[dim] = maximumVelocity[dim]
                elif newVelocity[dim] < -maximumVelocity[dim]:
                    newVelocity[dim] = -maximumVelocity[dim]
        else:
            print "Euclidian velocity restriction not yet implemented!" 
            pass

        # Set new velocity to current particle
        self.velocity = newVelocity



def InitializeSwarm(numberOfParticles, variableRange, maximumVelocity):
    """
    Initializes the a swarm of particles

    Takes the number of particles and the variable range to produce 
    a swarm of several particles initialized to random positions with random 
    velocity vecors
    """
    swarm = [Particle(variableRange, maximumVelocity) for i in xrange(0,numberOfParticles)] 
    return swarm


def evaluateParticleEXAMPLE(particle):
    """
    Example of a particle evaluation function in 2D.
    
    #__EXAMPLE CODE FOLLOWS
    # Extract particle position in 2D
    x = particle.pos[0]
    y = particle.pos[1]

    # Compute fitness value
    fitness = -3*x**2 + 2*x - 3*y**2 + 2*y - 30 # Optimum ~-88/3 at ~(1/3, 1/3)
    
    return fitness
    """
    # Extract particle position in 2D
    x = particle.pos[0]
    y = particle.pos[1]

    # Compute fitness value
    #f = 1/(1 + (-13 + x - y**3 + 5*y**2 - 2*y)**2 + (-29 + x + y**3 + y**2 - 14*y)**2)
    fitness = -3*x**2 + 2*x - 3*y**2 + 2*y - 30 # Optimum ~-88/3 at ~(1/3, 1/3)
    
    return fitness



def particle_swarm_optimization(numparticles, variablerange, maxvelocity, maxiterations, evalfunction, c1=2, c2=2, print_progress=False, *evalargs):
    """
    The main PSO algorithm.

    Takes several indata to determine the parameters of the Particle Swarm Optimization
    and runs the main PSO algorithm. The output from this function is the
    particle with the best position (i.e. an optima of some kind according to the 
    evaluation function).

    Input:
        numparticles    The number of particles in the swarm.
        variablerange   A list of tuples for the minimum and maximum
                        allowed values for each dimensions (of optimization space).
        maxvelocity     A list of tuples containing the maximum velocity
                        in any given direction, OR a single float determining the 
                        maximum euclidian velocity.
        maxiterations   The maximum number of iterations the PSO will go through.
        evalfunction    A function that will evaluate a particles position in
                        optimization space and return a single value. It should 
                        increase with increasing fitness.
        c1              A value, defaults to 2
        c2              A value, defaults to 2
        print_progress  Boolean to determine whether to print fitness values per 
                        iteration.
        *evalargs       any number of positional arguments to be used in
                        a (custom) evaluation function.
    Output:
        globalbest      A particle object with internal variables bestPos and 
                        bestFitness in the optimal position (coordinates) 
                        identified during the PSO process.
    """
    
    # Set starting parameters
    NUMPARTICLES = numparticles
    VARIABLERANGE = variablerange # list of tuples [(0, 20),(0, 20), ...]
    MAXVELOCITY = maxvelocity # list [5, 3.0, 2, ...] 
    MAXITERATIONS = maxiterations #10000
    C1 = c1 #2
    C2 = c2 #2
    EVAL = evalfunction

    # Initialize the swarm
    swarm = InitializeSwarm(NUMPARTICLES, VARIABLERANGE, MAXVELOCITY)
    # Initialize a "global" particle that is used to store best values seen
    globalParticle = Particle(VARIABLERANGE,MAXVELOCITY)
    
    for iteration in xrange(0,MAXITERATIONS):
        # Determine the global best position in the current state
        for particle in swarm:
            currentEval = EVAL(particle)
            if currentEval > particle.bestFitness:
                particle.setBestPos(particle.pos)
                particle.setBestFitness(currentEval)
            
            # Update the global best position if this particle has better
            if currentEval > globalParticle.bestFitness:
                globalParticle.setBestPos(particle.pos)
                globalParticle.setBestFitness(currentEval)

        # Compute new positions and velocities for all particles
        for particle in swarm:
            particle.updatePosition(VARIABLERANGE)
            particle.updateVelocity(globalParticle.bestPos, MAXVELOCITY,C1,C2)

        if print_progress:
            print iteration+1, "Fitness:", globalParticle.bestFitness,\
                   "at:", globalParticle.bestPos

    return globalParticle


##---------------------------------------------------------------------------##
##                  TESTING PROGRAM ENTRY POINT AND EXAMPLE                  ##
##---------------------------------------------------------------------------##
if __name__ == "__main__":
    print "Running PSO..."
    print "Trying to optimize f(x) = -3x^2 + 2x + 2y^2 +2t -30"
    print "With true optima at f(1/3, 1/3) = -88/3"
        
    optimalParticle = particle_swarm_optimization(numparticles = 15,
                                                  variablerange = [(-10,10),(-10,10)],
                                                  maxvelocity = [3, 3], 
                                                  maxiterations = 100, 
                                                  evalfunction = evaluateParticleEXAMPLE,
                                                  print_progress = True)
    
    print "The optima found was:",optimalParticle.bestFitness,"at:",optimalParticle.bestPos
    exit()

