'''
MATPOWER
Copyright (c) 1996-2016 by Power System Engineering Research Center (PSERC) by Ray Zimmerman, PSERC Cornell
This code follows part of MATPOWER.
See http://www.pserc.cornell.edu/matpower/ for more info.
Modified by Oak Ridge National Laboratory (Byungkwon Park) to be used in the parareal algorithm
'''
import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix, csc_matrix, lil_matrix, identity

def makeSbus(baseMVA, bus, gen):
    on = sp.find(gen[:, 7] > 0)[0]                 ## which generators are on?
    gbus = gen[on, 0].toarray().astype('int64')-1  ## what buses are they at?

    ## form net complex bus power injection vector
    nb = len(bus.toarray())
    ng = len(gen.toarray())

    row = gen[:, 0].toarray().reshape(-1)-1
    col = np.arange(ng)
    data = np.ones(ng).reshape(-1)
    Cg = csc_matrix((data, (row, col)), shape=(nb, ng))
    Sbus = (Cg * (gen[on, 1] + 1j*gen[on, 2]) - (bus[:, 2] + 1j*bus[:, 3]) ) / baseMVA
    Sbus = Sbus.toarray().reshape(-1)

    return Sbus


