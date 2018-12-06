# Chang Yoon Park Oct. 14, 2018
# This code is by no means optimal.

import numpy as np
import matplotlib.pyplot as plot

## Note: This is the 2D Wendland Kernel.
class WendlandKernel:

    def __init__ (self, dx, smoothingLength):
        self.dx = dx
        self.smoothingLength = smoothingLength

    def W_ij(self, pi, pj):
        dist = np.linalg.norm(pi.x - pj.x) + 1.0E-10
        h = self.smoothingLength * 0.5
        q = (dist / h)
        if (q > 2.0):
            return 0.0
        else:
            return (0.55704230082163367519109317180380 / (h * h)) * ((1.0 - 0.5 * q)**4.0) * (2.0 * q + 1.0);            		

    def gW_ij(self, pi, pj):
        dist = np.linalg.norm(pi.x - pj.x) + 1.0E-10
        h = self.smoothingLength * 0.5;
        q = (dist / h);
        if (q > 2.0):
            return np.array([0,0])
        else:			 
            return (((7.0/(4.0*np.pi))/(h*h*h)*(-5.0 * q)*(1.0 - 0.5 * q)**3.0)) * ((pi.x - pj.x) / dist);

class Particle:
    
    def __init__ (self, pid):
        self.x = 0
        self.vol = 0
        self.kernelSum = 0
        self.f = 0
        self.fdot = 0
        self.gradF = np.array([0.0,0.0])
        self.laplF = 0
        self.id = pid


class ParticleSystem:

    def __init__ (self):
        self.pList = []
        self.numParticles = 0


    def plot(self):
        self.fig, self.axarr = plot.subplots(2,4,figsize=(10,5))
        plot.subplots_adjust(top = 0.95, bottom=0.05, hspace=0.25, wspace=0.05)

        x = [p.x[0] for p in self.pList]
        y = [p.x[1] for p in self.pList]
        z1 = [p.f    for p in self.pList]
        z2 = [p.gradF[0] for p in self.pList]
        z3 = [p.gradF[1] for p in self.pList]
        z4 = [p.laplF    for p in self.pList]

        plotx1 = []
        ploty1 = []
        plotx2 = []
        ploty2 = []
        plotx3 = []
        ploty3 = []
        plotx4 = []
        ploty4 = []

        for p in self.pList:
            # print(abs(p.x[1] - 0.5))
            if abs(p.x[1] - 0.5) < 1.0E-5:
                plotx1.append(p.x[0])
                ploty1.append(p.f)
                plotx2.append(p.x[0])
                ploty2.append(p.gradF[0])
                plotx3.append(p.x[0])
                ploty3.append(p.laplF)

            if abs(p.x[0] - 0.5) < 1.0E-5:
                plotx4.append(p.x[1])
                ploty4.append(p.gradF[1])

        self.axarr[0, 0].scatter(x,y,c=z1,s=25)
        self.axarr[0, 0].axis('equal')
        self.axarr[0, 0].set_title('f(x,y)')

        self.axarr[1, 0].plot(plotx1,ploty1)
        self.axarr[1, 0].set_title('f(x,y=0.5)')
        # axarr[0, 1].axis('equal')


        self.axarr[0, 1].scatter(x,y,c=z2,s=25)
        self.axarr[0, 1].axis('equal')
        self.axarr[0, 1].set_title('(grad_f(x,y))_x')

        self.axarr[1, 1].plot(plotx2,ploty2)
        self.axarr[1, 1].set_title('grad_f(x,y)|y=0.5')


        self.axarr[0, 2].scatter(x,y,c=z3,s=25)
        self.axarr[0, 2].axis('equal')
        self.axarr[0, 2].set_title('(grad_f(x,y))_y')

        self.axarr[1, 2].plot(plotx4,ploty4)
        self.axarr[1, 2].set_title('grad_f(x,y)|x=0.5')


        self.axarr[0, 3].scatter(x,y,c=z4,s=25)
        self.axarr[0, 3].axis('equal')
        self.axarr[0, 3].set_title('lapl_f(x,y)')

        self.axarr[1, 3].plot(plotx3,ploty3)
        self.axarr[1, 3].set_ylim([0,5])
        self.axarr[1, 3].set_title('lapl_f(x,y)|y=0.5')


        plot.show()


    def getParticle(self,pID):
        return self.pList[pID]

    def addParticle(self, part):
        self.pList.append(part)
        self.numParticles += 1

    def forAllParticles(self, func):
        for particle in self.pList:
            func(particle)

    def forAllParticlesAndItsNeighbors(self, func):
        for particle_i in self.pList:
            for particle_j in self.pList:
                func(particle_i,particle_j)


