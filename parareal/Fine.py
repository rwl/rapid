import numpy as np
from mpi4py import MPI
from packer import packV, unpackV
from parareal.fn import fnRK4, fnTrap, fnTrap_adap

def safe_arange(start, stop, step):
    return step*np.arange(start / step, stop / step)

class Fine:
    def __init__(self, t1, t2, nStates, nNetwork, nsteps=1, itst=0, fault=[]):
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

        self.faults = fault

    def propagate(self, y, dbg=0):
        n=y.size
        if self.nVector is None:
            self.nVector = n
        elif ( self.nVector != n ):
            print("Error in vector length", self.nVector, "!=", n)

        yo = np.array(y)        # allocate new vector to return

        if not self.faults:
            print("Inside the Fine w/o fault")
            Xinit, Vbusinit, Ibus = unpackV(y, self.nStates, self.nNetwork)
            X, Vbus, Ibus = fnRK4(self.timeStart, self.timeEnd, self.nSubsteps, Xinit, Vbusinit)
        else:
            print("Inside the Fine with fault")
            X0, Vbus0, Ibus = unpackV(y, self.nStates, self.nNetwork)
            dt = (self.timeWindowEnd-self.timeWindowStart)/float(self.nSubsteps)
            nn = np.arange(0,self.nSubsteps,1)   # use single time increments
            print("nSubsteps ", self.nSubsteps)

            for ii in nn:
                t1 = self.timeWindowStart + float(ii)*dt
                t2 = t1 + dt
                ns = 1
                
                # processing of faults in debug loop 
                for ifault in self.faults:
                    r = ifault.run(t1)

                # other things to do
                X, Vbus, Ibus = fnRK4(t1, t2, ns, X0, Vbus0)
                X0 = X
                Vbus0 = Vbus

        # end new loop structure to enable faults
        yo=packV(X, Vbus, Ibus)

        self.pIteration += 1
        return yo
