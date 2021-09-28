'''
MATPOWER
Copyright (c) 1996-2016 by Power System Engineering Research Center (PSERC) by Ray Zimmerman, PSERC Cornell
This code follows part of MATPOWER.
See http://www.pserc.cornell.edu/matpower/ for more info.
Modified by Oak Ridge National Laboratory (Byungkwon Park) to be used in the parareal algorithm.
'''

import numpy as np
from scipy.sparse import csr_matrix, csc_matrix, lil_matrix, identity

def bustype(bus, gen):
    nb = len(bus.toarray())
    ng = len(gen.toarray())

    row = gen[:, 0].toarray().reshape(-1) - 1
    col = np.arange(ng)
    data = (gen[:, 0] > 0).toarray().reshape(-1)
    Cg = csc_matrix((data, (row, col)), shape=(nb, ng))
    bus_gen_status = Cg*np.ones(ng)

    ## form index lists for slack, PV, and PQ buses
    busidx = (bus.tocsc()[:, 1])
    busidx = busidx.todense().reshape(-1)

    ref = np.logical_and(busidx == 3, bus_gen_status > 0 )[0]
    ref = np.where(np.transpose(ref) == True)[0]

    pv = np.logical_and(busidx == 2, bus_gen_status > 0)[0]
    pv = np.where(np.transpose(pv) == True)[0]

    pq = np.logical_or(busidx == 1, bus_gen_status == 0)[0]
    pq = np.where(np.transpose(pq) == True)[0]

    return ref, pv, pq


