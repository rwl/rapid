import numpy as np
# fault functions
import ffault

# fault objects
class Fault:
    def __init__(self, t1, t2, typ='bus', par=[], linetrip=[]):
        self.timeStart       = t1
        self.timeEnd         = t2

        # type of fault
        self.ftype           = typ
        # fault parameters
        self.fparm           = par
        # line trip
        self.ltrip           = linetrip

        # selection of fault function
        # selection is based on typ parameter
        # typ has to exist in ffault.py, dictionary fault_dict
        # there is no error checking on selection, but we will
        # not have many fault types anyway
        # function that executes during the fault
        # if faults are discrete events, this can be just a dummy function
        self.ffun            = ffault.fault_dict(typ)

        # if there are special things to do on the fault start and fault end,
        # we can add other functions here
        # that can be used for specific tasks associated with
        # starting and ending the fault, such as
        self.ffun_start      = ffault.fault_start_dict(typ)
        self.ffun_end        = ffault.fault_end_dict(typ)

        self.started         = 0 # indicator to check if the fault has started
        self.ended           = 0 # indicator to check if the fault has ended

    def run(self, t):
        result = 0

        if self.ended == 1:  # if fault already ended, do nothing
            return result

        if t >= self.timeStart and self.started == 0: # fault starts
            print("start fault",self.ftype)
            self.started = 1
            # call function to do something special for fault start
            result = self.ffun_start(t, self.fparm, self.ltrip)
            return result
            
        if self.started == 1 and t < self.timeEnd: # during fault  -- Is it necessary ??
            #print("run fault",self.ftype)
            result = self.ffun(t, self.fparm, self.ltrip)
            return result
        
        if t >= self.timeEnd: # fault ends
            print("end fault",self.ftype)
            self.ended = 1
            # call function to do something special for fault end
            result = self.ffun_end(t, self.fparm, self.ltrip)
            return result

        return result
