import numpy as np
def id2ctype(id):
        if id == -1:
            ctype = 'Low Seis.'
        elif id == 0:
            ctype = 'Undefined'
        elif id == 1:
            ctype = 'Mainshock*'
        elif id == 2:
            ctype = 'Doublet'
        elif id == 3:
            ctype = 'Triplet'
        elif id == 4:
            ctype = 'Swarms I'
        elif id == 5:
            ctype = 'Swarms II'
        return ctype

def sequence_type(m):
    m = np.array(m)
    m.sort()
    m = m[::-1]
    
    if len(m) >=4:
        m1 = np.round(m[0],1)
        m2 = np.round(m[1],1)
        m3 = np.round(m[2],1)
        m4 = np.round(m[3],1)

        if np.round(m1-m2,1) >= 0.5:
            type = 1 # Mainshock
        elif np.logical_and(np.round(m1-m2,1) <= 0.4, np.round(m2-m3,1)>=0.5):
            type = 2 # Doublets
        elif np.logical_and(np.round(m1-m3,1)<=0.4, np.round(m3-m4,1)>=0.5):
            type = 3 # Triplets
        elif np.round(m1-m4,1)<=0.5:
            type = 4 # Swarm -I
        elif np.round(m1-m4,1)<=1.0:
            type = 5 # Swarm -II
        else:
            type = 0 # un-defined
            
    else:
        type = -1 # low seismicity
            

    return type


