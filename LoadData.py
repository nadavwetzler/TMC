import numpy as np
import pandas as pd
from obspy import UTCDateTime

class CCats():
    pass

def datenum2time(datenum):
    t = datenum.apply(lambda x: (x-719529.0)*24*3600)
    A = t.apply(lambda x: UTCDateTime(x))
    A = np.asarray(A)
    return A

def trimcat(cat0, max_depth, lats , lons, years, Mc):

        Idepth = cat0.Depth < max_depth
        Ilat1 = cat0.Lat >= min(lats)
        Ilat2 = cat0.Lat <= max(lats)
        Ilon1 = cat0.Long >= min(lons)
        Ilon2 = cat0.Long <= max(lons)
        Ilat = np.logical_and(Ilat1, Ilat2)
        Ilon = np.logical_and(Ilon1, Ilon2)
        Ilatlon = np.logical_and(Ilat, Ilon)
        It1 = cat0.ot >= UTCDateTime(min(years), 1, 1)
        It2 = cat0.ot <= UTCDateTime(max(years), 1, 1)
        It = np.logical_and(It1, It2)
        I = np.logical_and(Ilatlon, It)
        I = np.logical_and(I, Idepth)
        Im = cat0.M >= Mc
        I = np.logical_and(I, Im)

        cat0.Lat = cat0.Lat[I]
        cat0.Long = cat0.Long[I]
        cat0.ot = cat0.ot[I]
        cat0.M = cat0.M[I]
        cat0.Depth = cat0.Depth[I]
        cat0.N = cat0.N[I]
        cat0.datenum = cat0.datenum[I]
        return cat0

def read_SCEDC(file, lats, lons, Mc, years, max_depth):
    CAT = CCats()
    # Datetime,ET,GT,MAG,M,LAT,LON,DEPTH,Q,EVID,NPH,NGRM,datenum
    A = pd.read_csv('data/%s' % file)
    CAT.Depth = np.asarray(A.DEPTH)
    CAT.Long = np.asarray(A.LON)
    CAT.Lat = np.asarray(A.LAT)
    CAT.M = np.asarray(np.round(A.MAG, 1))
    CAT.ot = datenum2time(A.datenum)
    CAT.N = np.arange(0,len(CAT.M), 1)
    CAT.datenum = np.asarray(A.datenum)
    print('file %s was loaded with %d events' % (file, len(CAT.Long)))
    CAT = trimcat(CAT, max_depth, lats , lons, years, Mc)
    CAT.name = file.split('.')[0]
    print('Cat %s was trimmed to %d events' % (CAT.name, len(CAT.Long)))
    return CAT
