# Chang Yoon Park Oct. 14, 2018
# This code is by no means optimal.

from helpers import *

# Define Particle system        
pSys = ParticleSystem()

# Scalar function f(x,y)
def assignField(particle):
    particle.f = particle.x[0]**2 + particle.x[1]**2 

# init gradient
def initStuff(particle):
    particle.kernelSum = 0.0
    particle.volume = 0.0
    particle.gradF = np.array([0.0,0.0])
    particle.laplF = 0.0

def computeKernelSum(pi, pj):
    pi.kernelSum += kernel.W_ij(pi, pj)
    # print(pi.id, pi.x, pj.id, pj.x, pi.x - pj.x)
    # print(pi.kernelSum)

def computeVolume(particle):
    particle.vol = 1./particle.kernelSum

# Procedure to compute grad_f(x,y), lapl_f(x,y)
def computeGradientAndLaplacian(pi,pj):

    if pi is not pj:        

        gW_ij = kernel.gW_ij(pi,pj)
        dist  = np.linalg.norm(pi.x - pj.x)
        dirVec = (pi.x - pj.x) / dist
        pi.gradF += (pj.f - pi.f) * pi.vol * gW_ij
        pi.laplF += 2.0 * pi.vol * ((pi.f - pj.f) / dist) * np.dot(dirVec, gW_ij)


# Assign initial temperature field.
def assignInitTemps(particle):
    if particle.x[0] < 0 or particle.x[1] > 1:
        particle.f = 1
    else:
        particle.f = 0

# Update Temperature
def updateTemp(particle):
    if particle.x[0] < 0 or particle.x[1] > 1:
        return

    particle.f += particle.laplF * 0.01

def find2DExample():
    # Add particles to the particle system.
    pid = 0
    for i in range(-3,21):
        for j in range(-3,21):
            p = Particle(pid)
            p.x = np.array([-i * dx, j * dx])
            pSys.addParticle(p)
            pid += 1

    pSys.forAllParticles(assignField)
    pSys.forAllParticles(initStuff)
    pSys.forAllParticlesAndItsNeighbors(computeKernelSum)
    pSys.forAllParticles(computeVolume)
    pSys.forAllParticlesAndItsNeighbors(computeGradientAndLaplacian)
    pSys.plot()


def heat2DExample():

    # Add particles to the particle system.
    pid = 0
    for i in range(-2,11):
        for j in range(0,12):
            p = Particle(pid)
            p.x = np.array([dx / 2.0 + i * dx, dx / 2.0 + j * dx])
            pSys.addParticle(p)
            pid += 1

    pSys.forAllParticles(assignInitTemps)
    pSys.forAllParticles(initStuff)
    pSys.forAllParticlesAndItsNeighbors(computeKernelSum)
    pSys.forAllParticles(computeVolume)
    for t in range(0,1000):
        pSys.forAllParticlesAndItsNeighbors(computeGradientAndLaplacian)
        pSys.forAllParticles(updateTemp)


# dx = 0.05
# smoothingLength = 0.1
# kernel = WendlandKernel(dx = dx, smoothingLength = smoothingLength)
# heat2DExample()


dx = 0.5
smoothingLength = 2.0
kernel = WendlandKernel(dx = dx, smoothingLength = smoothingLength)
find2DExample()