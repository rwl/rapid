import numpy as np
from mpi4py import MPI
from packer import packV, unpackV
from parareal.fn import fnRK4, fnTrap, fnTrap_adap, fnADM, fnHAM
##AS from trap_integrator import TrapIntegrator 

class Coarse:
    def __init__(self, t1, t2, nStates, nNetwork, nsteps=1, itst=0):
        comm = MPI.COMM_WORLD
        self.rank            = comm.Get_rank()
        self.nProcessors     = comm.Get_size()
        self.nVector         = nStates+nNetwork*4
        self.nStates         = nStates
        self.nNetwork        = nNetwork
        self.pIteration      = itst

        self.timeWindowStart = float(t1)
        self.timeWindowEnd   = float(t2)
        self.timeWindow      = float(t2-t1)
        self.nSubsteps       = nsteps

        self.timeInterval    = self.timeWindow/float(self.nProcessors)
        self.timeIncrement   = self.timeInterval/float(self.nSubsteps)
        self.timeStart       = self.timeWindowStart + float(self.rank)*self.timeInterval
        self.timeEnd         = self.timeStart+self.timeInterval

        print("Coarse Rank ", self.rank, " Size ", self.nProcessors, " timeInterval ", self.timeInterval)

    def set_time_integrator(self, integratorName):
        self.timeIntegrator = integratorName

    def propagate(self, y):
        n=y.size
        if self.nVector is None:
            self.nVector = n
        elif ( self.nVector != n ):
            print("Error in vector length", self.nVector, "!=", n)

        yo = np.array(y)        # allocate new vector to return

        Xinit, Vbusinit, Ibusinit = unpackV(y, self.nStates, self.nNetwork)

##AS     X, Vbus = TrapIntegrator.advance_state(self.timeStart, self.timeEnd, self.nSubsteps, Xinit, Vbusinit)
#        X, Vbus = fnRK4(self.timeStart, self.timeEnd, self.nSubsteps, Xinit, Vbusinit)
        if self.timeIntegrator == 'trap':
             X, Vbus, Ibus = fnTrap(self.timeStart, self.timeEnd, self.nSubsteps, Xinit, Vbusinit)
        elif (self.timeIntegrator == 'adm'):     
             X, Vbus, Ibus = fnADM(self.timeStart, self.timeEnd, self.nSubsteps, Xinit, Vbusinit)
        elif (self.timeIntegrator == 'ham'):     
             X, Vbus, Ibus = fnHAM(self.timeStart, self.timeEnd, self.nSubsteps, Xinit, Vbusinit)
        else:
            raise NotImplementedError("You should implement this!")

        yo=packV(X, Vbus, Ibus)

        self.pIteration += 1
        return yo
