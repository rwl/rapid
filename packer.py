import numpy as np

# Entries in x and vbus vectors, 1 column (N,1) matrix
# ([1],[2],...,[N]), vertical vectors
# The entries are packed into flattened 1 row matrix
def packV(x, vbus, ibus): 
    v = np.concatenate((np.ravel(x),np.ravel(vbus.real),np.ravel(vbus.imag),np.ravel(ibus.real),np.ravel(ibus.imag))) 
    return v

# Unpacking from flat 1 row matrix into
# two vectors, 1 column matrices
def unpackV(v, nstat, nnetw):
    iPtr1 = nstat
    iPtr2 = iPtr1+nnetw
    iPtr3 = iPtr2+nnetw
    iPtr4 = iPtr3+nnetw 
    iPtr5 = iPtr4+nnetw 
    x = v[0:iPtr1]
    vbus = np.vectorize(complex)(v[iPtr1:iPtr2],v[iPtr2:iPtr3])
    ibus = np.vectorize(complex)(v[iPtr3:iPtr4],v[iPtr4:iPtr5]) 
    x = x.reshape((-1,1))
    vbus = vbus.reshape((-1,1))
    ibus = ibus.reshape((-1,1)) 
    return x, vbus, ibus 

