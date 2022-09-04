import numpy as np


def WnCfaultL(m, ifault=0):
# Following Wells and Coppersmith 1994
    if ifault == 1:  # strike-slip
        al = 4.33
        bl = 1.49
        aw = 3.80
        bw = 2.59
    elif ifault == 3:  # reverse
        al = 4.49
        bl = 1.49
        aw = 4.37
        bw = 1.95
    elif ifault == 2:  # normal
        al = 4.34
        bl = 1.54
        aw = 4.04
        bw = 2.11
    else:
        al = 4.38
        bl = 1.49
        aw = 4.06
        bw = 2.25

    l = 10**((m - al)/bl) * 1000
    w = 10**((m - aw)/bw) * 1000
    return l, w


def AdjustEffectiveRadius(M):
    if np.size(M) > 1:
        M = np.array(M)
        fact_r = np.ones(len(M)) * 2
        fact_r[M >= 9.0] = 1
        fact_r[M >= 8.0] = 1.5
    else:
        if M >= 9.0:
            fact_r = 1
        elif M > 8.0:
            fact_r = 1.5
        else:
            fact_r = 2.0
    return fact_r


def windowsTMC(M):
    # t = 0.002*M**7.5
    # t = 0.75*np.exp(M)
    t = 0.1 * np.exp(M)
    # t[t > 365] = 365
    [r_predicted, _] = WnCfaultL(M, ifault=0)
    d = (AdjustEffectiveRadius(M) * r_predicted) / 1000
    t = np.array(np.ceil(t))
    d = np.array(np.ceil(d))
    t = t.astype(int)
    d = d.astype(int)
    return t, d


def MakeCircleUTM(rm, lon0, lat0):
    e2u_zone = int(divmod(lon0, 6)[0])+31
    e2u_conv = Proj(proj='utm', zone=e2u_zone, ellps='WGS84')
    #Apply the converter
    x0, y0 = e2u_conv(lon0, lat0)
    a = np.arange(0, 360, 1)
    a = np.deg2rad(a)
    Yo1 = rm*np.cos(a) + y0
    Xo1 = rm*np.sin(a) + x0
    lon, lat = e2u_conv(Xo1, Yo1, inverse=True)
    xy = np.zeros((len(Yo1), 2))
    xy[:, 0] = lon
    xy[:, 1] = lat
    return lon, lat, xy

def CheckInListC(list1, ClusterList):
    llist = len(list1)
    cols = []
    # run over all events in new list
    for jj in range(llist):
        # check in all clusters in ClusterList
        for mm in range(len(ClusterList)):
            # locate position of event in Cluster
            cols1 = np.argwhere(np.array(ClusterList[mm]) == list1[jj])
            # if event is found in a Cluster
            if len(cols1) > 0:
                # save the Cluster id
                for nn in range(len(cols1)):
                    cols.append(mm)
    cols = np.unique(cols)
    return cols

def TMC(cat0):
    sortM = np.argsort(cat0.M)
    sortM = sortM[::-1]
    cat0.M = cat0.M[sortM]
    cat0.Depth = cat0.Depth[sortM]
    cat0.Lat = cat0.Lat[sortM]
    cat0.Long = cat0.Long[sortM]
    cat0.N = cat0.N[sortM]
    cat0.datenum = cat0.datenum[sortM]
    cat0.ot = cat0.ot[sortM]
    len1 = len(cat0.M)

    [tw_AS, rpt] = windowsTMC(cat0.M)
    Mat_dt = np.zeros((len1, len1))

    for ii in range(len1):
        Mat_dt[ii, :] = cat0.datenum - cat0.datenum[ii]
    Mat_dt = np.ceil(Mat_dt)
    Mat_dt.astype(int)

    print('Calculating distances...')
    MatDist = np.zeros((len1, len1)).astype(int)
    Mat_AS = np.zeros((len1, len1)).astype(int)

    lonlat = np.zeros((len(cat0.Lat), 2))
    lonlat[:, 0] = cat0.Long
    lonlat[:, 1] = cat0.Lat
    for ii in range(len1):
        [_, _, xyr] = MakeCircleUTM(rpt[ii]*1000, cat0.Long[ii], cat0.Lat[ii])
        p = path.Path(xyr)
        MatDist[ii, :] = p.contains_points(lonlat)
        Mat_AS[ii, :] = np.logical_and(Mat_dt[ii, :] > 0, Mat_dt[ii, :] < tw_AS[ii])
    I = np.logical_and(MatDist, Mat_AS)
    
    print('Clustering...')
    IDS = np.zeros((sum(sum(I)), 5))
    IDS[:, 4] = Mat_dt[I]
    del Mat_dt
    IDS[:, 2] = MatDist[I]
    del MatDist

    [iq, jq] = np.meshgrid(np.arange(0, len1), np.arange(0,len1))

    IDS[:, 0] = iq[I]
    del iq
    IDS[:, 1] = jq[I]
    del jq
    IDV = np.unique(IDS[:, 0:2])
    ClusterList = []
    print('Associating events to clusters')
    for ii in range(len(IDV)):
        ID = IDV[ii]
        # Make list for all events in row2 associated with ID
        list1 = IDS[IDS[:, 0] == ID, 1]

        # add ID to list
        list1 = np.append(list1, ID)

        # Check for ID in Clusterlist
        cols = CheckInListC(list1, ClusterList)

        if len(cols) == 0: # NOT exists list
            list0 = np.unique(list1)

        else: # Exists in list
            list0 = []
            for mm in range(len(cols)):
                list2 = ClusterList[cols[mm]]
                for vv in range(len(list2)):
                    list0.append(list2[vv])
            for vv in range(len(list1)):
                list0.append(list1[vv])
            list0 = np.array(list0)
            list0 = np.unique(list0)

            # remove colums
            ClusterList1 = []
            for nn in range(len(ClusterList)):
                if nn not in cols:
                    ClusterList1.append(ClusterList[nn])
            ClusterList = ClusterList1

        ClusterList.append(np.array(list0))
        print('%d / %d %d' % (ii, len(IDV), len(ClusterList)))

    print('Add cluster info to EQ table')

    # all events
    c = np.zeros(len(cat0.Lat))
    n = np.zeros(len(cat0.Lat))
    posM = []
    clist = []
    for ii in range(len(ClusterList)):
        posC = ClusterList[ii].astype(int)
        c[posC] = ii+1
        n[posC] = len(ClusterList[ii])
        magsC = cat0.M[ClusterList[ii].astype(int)]
        posM.append(posC[np.argmax(magsC)])
        clist.append(ii+1)
    cat0.n = n
    cat0.c = c
    
    return cat0

