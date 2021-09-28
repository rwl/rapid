import csv
import numpy as np
import cProfile
import sys
from mpi4py import MPI
from parareal.Coarse import Coarse
from parareal.Fine import Fine
# read raw data and obtain an initial condition
from initialization.PowerModel import X0, Vbus0, Ibus, idxdelta
from Fault import Fault
from packer import packV
from pdefault import pDefaults

def profile(filename=None, comm=MPI.COMM_WORLD):
    def prof_decorator(f):
        def wrap_f(*args, **kwargs):
            pr = cProfile.Profile()
            pr.enable()
            result = f(*args, **kwargs)
            pr.disable()

            if filename is None:
                pr.print_stats()
            else:
                filename_r = filename + ".{}".format(comm.rank)
                pr.dump_stats(filename_r)

            return result
        return wrap_f
    return prof_decorator


@profile(filename="profile_out")
def parareal_main(md, comm=MPI.COMM_WORLD):

    id    = comm.Get_rank()  # processor id = coarse interval id
    N     = comm.Get_size()  # number of processors
    k     = 0                # iteration

    t1=md.t1
    t2=md.t2
    tol=md.tol
    tolcheck=md.tolcheck
    nCoarse=md.nCoarse
    nFine=md.nFine

    dbg=md.dbg
    outfile=md.outfile

    tf1 = 0
    tf2 = 0
    # Create fault objects
    # fault objects defined in Fault.py
    faults = []
    if md.fault:
        for ifault in md.fault['faults']:
            print("reading faults is not empty") 
            tf1 = ifault['start']
            tf2 = ifault['end']
            ftyp = ifault['type']
            fpar = ifault['index']
            linetrip = ifault['linetrip']
            # append Faults definition to faults object (e.g., print(faults[0].timeStart))
            faults.append(Fault(tf1,tf2,ftyp,fpar,linetrip))   

#    Cosimulation = 1 
#    idxDistBus = [2]
#    distSystem = []
    if md.dist:
       for idist in md.dist['distribution_network']:
           print("reading distribution network") 
           Cosimulation = idist['active']
           DistBus = idist['index']
           idxDistBus = list(np.array(DistBus) - 1)
#           distSystem.append(distNtk())

    # Control flags
    convergeNext = 0
    converge     = 0
    complete     = 0

    # default, lets all iterations run
    maxIteration = N+1   

    # PowerGrid
    nStates  = X0.size
    nNetwork = Vbus0.size
    nTotal   = nStates + nNetwork*4
    
    # For gathering solution from all processors
    solF = None
    y1 = None

    nVector  = nTotal 
    idxchk   = idxdelta

    # create time axis
    dx = (t2-t1)/float(N)
    x=np.arange(t1, t2+dx/2, 0.2)              
    # create sliding time window
    fWindow = x[np.logical_and(x>=tf1, x<=dx+tf2)]  
    print("sliding time Window", fWindow)

    ## start the sequential approach for the during fault simulation
    if id == 0:
        y0=packV(X0, Vbus0, Ibus)
        print("rank 0 : y0 type ", y0.dtype," size ", y0.size)

        solF = np.empty([len(fWindow)-1, nVector], dtype=np.float64)

        y1 = y0
        ti1 = fWindow[0]
        for i in range(len(fWindow)-1):
            ti2 = fWindow[i+1]
            if dbg > 0: print("i ", i, "t",ti1,ti2,nFine)
            fIntegrator  =  Fine(ti1, ti2, nStates, nNetwork, nFine, fault=faults)
            # simulation between saving time
            solF[i] = fIntegrator.propagate(y1, dbg=md.dbg)        
            y1 = solF[i]
            ti1 = ti2
    else:
        y0 = np.empty(nTotal)
        y1 = y0

    ## define the fine and the coarse operator
    mCoarse=Coarse(fWindow[-1], t2, nStates, nNetwork, nCoarse)
    mCoarse.set_time_integrator('trap')
    mFine  =  Fine(fWindow[-1], t2, nStates, nNetwork, nFine, fault=None)

    solAll = None
    if id == 0:
        solAll = np.empty([N, nVector], dtype=np.float64)

    if dbg == 1 and id == 0:
        print("Proc", id, "index check", idxchk)

    ## start the parareal algorithm
    k = 0
    if(id==0):
        solU = mCoarse.propagate(y1)      # this will be U_id+1_k-1
        print(" solU type : ", solU.dtype," size ", solU.size)
        comm.Send(solU, dest=id+1, tag=k) # sends Ut_1_0 to proc 1,
                                          # blocking send, can reuse
                                          # buffer for fine propagation next

        solU = mFine.propagate(y1)        # this will be the final
                                          # solution for the first
                                          # interval
        converge=1
        comm.send(converge, dest=id+1, tag=1000)
        comm.Send(solU,     dest=id+1, tag=2000)
        complete=1
        if dbg == 1: print("Proc", id, "iteration", k)
        # print "id solU", id, solU
        comm.Gather(solU, solAll, root=0) # wait until all the intervals
                                          # send their solutions
    else:
        comm.Recv(y1,source=id-1, tag=k) # receive coarse solution from
                                          # id-1, y0 is used as a receive
                                          # buffer
        # U_i_km = y0                       # save initial value for k!=0
        U_j_km = mCoarse.propagate(y1)    # propagate y0, y0 = U_i_km for k=1
        Ut_j_km=np.copy(U_j_km)           # create view for k=1
        if(id!=(N-1)):
            comm.Send(U_j_km, dest=id+1, tag=k) # send solution to proc id+1
        if dbg == 1: print("Proc", id, "iteration", k)
        if dbg >  1: print("Proc", id, "iteration", k, "vectors",   y0, U_j_km)

    # end of iteration k=0

    for k in range(1,N+1):
        if dbg == 1 and id == (N-1): print("Iteration", k)
        if(complete == 0):
            Uh_j_km = mFine.propagate(y1)      # propagate fine solution,
                                               # y0 is kept for receiving
                                               # current initial values
            if(convergeNext==1):  # this fine propagation started from
                                  # full fine simulation from the
                                  # beginning of the window
                if dbg == 1: print("Proc", id, "iteration", k, "In convergeNext")
                converge=1
                if(id != (N-1) ):
                    comm.send(converge, dest=id+1, tag=1000)
                    comm.Send(Uh_j_km,  dest=id+1, tag=2000)
                solU=Uh_j_km      # not copy, just pointer, does not
                                  # matter since this processor is done
                complete=1
                # print "id solU", id, solU
                comm.Gather(solU, solAll, root=0)
                break

            converge=comm.recv(source=id-1, tag=1000)
            comm.Recv(y1, source=id-1, tag=2000)  
    
            Ut_j_k = mCoarse.propagate(y1) # propagate coarse solution for
                                           # k iteration

            U_j_k = Ut_j_k + Uh_j_km - Ut_j_km  # new solution for id+1,
                                                # and k iteration
                                                # perhaps add another function
            Udiff=U_j_k[idxchk]-U_j_km[idxchk]  # change from the previous iteration
                                                # for defined indices

            if tolcheck == 'L2':
                Uerr = np.sqrt(Udiff.dot(Udiff))  # L2 error magnitude,

            if tolcheck == 'maxabs':
                Uerr = np.amax(np.absolute(Udiff)) # max absolute difference

            if dbg == 1: print("Proc", id, "iteration", k, "Uerr", Uerr)
            if(k > maxIteration):  # stop the loop gracefully and save the result
                tol = Uerr*2.0
                converge=1
            
            if(converge==1 and Uerr >= tol):
                converge=0
                convergeNext=1

            if(id!=(N-1)):
                comm.send(converge,  dest=id+1, tag=1000)
                comm.Send(U_j_k,     dest=id+1, tag=2000)

            if(converge==1):
                if dbg == 1: print("Proc", id, "iteration", k, "In converge", Uerr)
                solU=U_j_k     # not copy, just pointer, does not matter
                               # since this processor is done
                complete=1
                # print "id solU", id, solU
                comm.Gather(solU, solAll, root=0)
                break

            U_j_km  = np.copy( U_j_k)     # to use for convergence check
                                          # in the next iteration
            Ut_j_km = np.copy(Ut_j_k)     # to use for calculating new
                                          # solution in the next iteration

    if( id==0 ):
        dx=(t2-fWindow[-1])/float(N)
        x=np.arange(fWindow[-1], t2, dx)              # create time axis
        x = np.insert(x, 0, 0, axis=0)
        print("length of window : ", len(x), " with intervals ", x)

        print("shaper of solAll : ", solAll.shape)
        solAll = np.insert(solAll, 0, y1, axis=0) # insert initial value
        solAll = np.insert(solAll, 0, x, axis=1)  # insert time axis

        idxiwant = np.insert(idxdelta, 0, -1)  # idxiwant  is indices that i am interested in
        solAll = solAll[:, idxiwant + 1]       # interested in generator rotor angles

        if outfile:
            fo = open(outfile, 'w', newline='')
        else:
            fo = sys.stdout
        csvwriter = csv.writer(fo, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        csvwriter.writerows(solAll)
        fo.close()

if __name__ == "__main__":
    comm=MPI.COMM_WORLD
    id    = comm.Get_rank()  # processor id = coarse interval id
    N     = comm.Get_size()  # number of processors
    md=pDefaults(id=id, N=N)
    md.parse_arguments()
    parareal_main(md)
