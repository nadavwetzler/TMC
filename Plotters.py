from osm import OSM as OSM
import matplotlib.pyplot as plb
import numpy as np
import glob
from obspy import UTCDateTime
import os
from TMC import WnCfaultL, AdjustEffectiveRadius, DistLatLonUTM
from obspy.geodetics import kilometers2degrees
from mpl_toolkits.basemap import Basemap
from Sequence_type import id2ctype
import time


def set_box_aspect(ax, ratio):
    #ratio = 0.3
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    # the abs method is used to make sure that all numbers are positive
    # because x and y axis of an axes maybe inversed.
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)


def plotMapC(ax, clust, Lats, Lons, Lon0, Lat0, ot, m0, cat, years):
    from TMC import WnCfaultL, AdjustEffectiveRadius, MakeCircleUTM
    # years = 1
    ylong = years*365  # years
    t0 = UTCDateTime(ot)
    T = (cat.ot - t0)/(3600*24)
    # catalog time
    Ilon = np.logical_and(cat.Long > min(Lons), cat.Long < max(Lons))
    Ilat = np.logical_and(cat.Lat > min(Lats), cat.Lat < max(Lats))
    Ireg = np.logical_and(Ilat, Ilon)
    ItlongAS = np.logical_and(cat.ot > t0, cat.ot < t0 + ylong*3600*24)
    ItlongFS = np.logical_and(cat.ot < t0, cat.ot > t0 - ylong*3600*24)

    dlat = np.abs(np.diff(Lats))[0]
    dlon = np.abs(np.diff(Lons))[0]
    m = Basemap(ax=ax,  projection='merc', llcrnrlat=min(Lats), urcrnrlat=max(Lats), llcrnrlon=min(Lons), urcrnrlon=max(Lons), lat_ts=1, resolution='i')
    # m = Basemap(ax=ax,  projection='cyl', llcrnrlat=min(Lats), urcrnrlat=max(Lats), llcrnrlon=min(Lons), urcrnrlon=max(Lons), resolution='i')
    # m = Basemap(ax=ax,  projection='lcc', llcrnrlat=min(Lats), urcrnrlat=max(Lats), llcrnrlon=min(Lons), urcrnrlon=max(Lons), resolution='i')
    # m = Basemap(ax=ax,  projection='tmerc', llcrnrlat=min(Lats), urcrnrlat=max(Lats), llcrnrlon=min(Lons), urcrnrlon=max(Lons), resolution='i', lon_0=Lon0, lat_0=Lat0)
    try:
        m.drawcoastlines(color='k', linewidth=0.2)
    except:
        pass
    try:
        m.drawmapboundary(zorder=0)
    except:
        pass
    try:
        m.fillcontinents(color='#cc9955', alpha=0.4, zorder=1)
    except:
        pass
    try:
        m.drawmapscale(lon=Lons[0]+0.2*dlat, lat=Lats[0]+0.1*dlat, lon0=x0, lat0=y0, length=int((2*dlat*0.1)*111), units='km')
    except:
        pass

    m.drawparallels(np.arange(Lats[0], Lats[1], 0.5*dlat), labels=[1,0,0,0], zorder=1)
    m.drawmeridians(np.arange(Lons[0], Lons[1], 0.5*dlon), labels=[0,0,0,1], zorder=1)
    # m.fillcontinents(color='coral',lake_color='aqua')

    # plot mainshock
    x0, y0 = m(Lon0, Lat0)
    ax.scatter(x0, y0, 200, marker='*', color='m', alpha=1.0, label='Mainshock',zorder=2)

    # plot foreshocks/aftershocks from the entire catalog
    xbg, ybg = m(cat.Long, cat.Lat)
    ax.scatter(xbg[np.logical_and(ItlongAS, Ireg)], ybg[np.logical_and(ItlongAS, Ireg)], 20, marker='o', color='grey', label='AS.Cat %d yr' % years, zorder=2)
    ax.scatter(xbg[np.logical_and(ItlongFS, Ireg)], ybg[np.logical_and(ItlongFS, Ireg)], 20, marker='^', color='grey', label='FS.Cat %d yr' % years, zorder=2)

    # cluster time
    Time = []
    for ii in range(len(clust['ot'])):
        Time.append(UTCDateTime(clust['ot'].values[ii]))
    Time = np.array(Time)
    Time = (Time - t0) / (3600*24)

    # plot foreshocks/aftershocks from the cluster
    selFS = Time < 0
    selAS = Time > 0
    x, y = m(clust['Long'].values, clust['Lat'].values)
    ax.scatter(x[selAS], y[selAS], 25, marker='o', facecolors='None', edgecolor='r', alpha=0.8, label='AS.Cl %d' % sum(selAS), zorder=3)
    ax.scatter(x[selFS], y[selFS], 25, marker='^', facecolors='None',edgecolor='b',alpha=0.8, label='FS.Cl %d' % sum(selFS), zorder=4)


    if m0 > 0:
        [r_predicted, _] = WnCfaultL(m0, ifault=0)
        fact_r = AdjustEffectiveRadius(m0)
        r_predicted = fact_r * r_predicted
        Ro = kilometers2degrees(r_predicted / 1000)
        [xr, yr, xyr] = MakeCircleUTM(r_predicted, Lon0, Lat0)
        xr0, yr0 = m(xr, yr)
        ax.plot(xr0, yr0, c='m', label='E.R.: %d km' % int(r_predicted / 1000))
        # ax.add_patch(mpatches.Circle(xy=[Lon0, Lat0], radius=Ro, color='red', alpha=0.3, zorder=30))

    # ax.legend(loc='upper left')

def scaleT(T, len, t0, t1):
        T = T[np.logical_and(T>t0, T<t1)] - t0
        T = T * (len / t1) + t0
        return T

def plotLatTimeClust(ax, clust, ot,Lat0, M0, cat, Lats,Lons):
    years = 1
    ylong = years*365 # years
    t0 = UTCDateTime(ot)
    T = (cat.ot - t0)/(3600*24)
    # catalog time
    Ilon = np.logical_and(cat.Long>min(Lons), cat.Long<max(Lons))
    Ilat = np.logical_and(cat.Lat>min(Lats), cat.Lat<max(Lats))
    Ireg = np.logical_and(Ilat, Ilon)
    It60 = np.logical_and(cat.ot > t0- 60*3600*24, cat.ot < t0 + 60*3600*24)
    ItlongAS = np.logical_and(cat.ot > t0 + 60*3600*24, cat.ot < t0 + (ylong)*3600*24)
    ItlongFS = np.logical_and(cat.ot < t0 - 60*3600*24, cat.ot > t0 - (ylong)*3600*24)

    # cluster time
    Time = []
    for ii in range(len(clust['Time'])):
        Time.append(UTCDateTime(clust['Time'].values[ii]))
    Time = np.array(Time)
    Time = (Time - t0) / (3600*24)

    selFS60 = np.logical_and(Time < 0, Time > -60)
    selFSlong = np.logical_and(Time <= -60, Time > -(ylong))
    selAS60 = np.logical_and(Time > 0, Time < 60)
    selASlong = np.logical_and(Time >= 60, Time < (ylong))


    ax.scatter(T[np.logical_and(It60,Ireg)], cat.Lat[np.logical_and(It60,Ireg)],15,marker='o', color='gray',alpha=0.6)
    ax.scatter(scaleT(T[np.logical_and(ItlongAS,Ireg)],30,60,ylong), cat.Lat[np.logical_and(ItlongAS,Ireg)],15,marker='o', color='gray',alpha=0.6)
    ax.scatter(-scaleT(-T[np.logical_and(ItlongFS,Ireg)],30,60,ylong), cat.Lat[np.logical_and(ItlongFS,Ireg)],15,marker='o', color='gray',alpha=0.6)

    ax.scatter(Time[selFS60], clust['Lat'][selFS60].values,s=25, facecolors='none', edgecolors='b',label='FS: %d' % sum(selFS60))
    ax.scatter(Time[selAS60], clust['Lat'][selAS60].values,s=25, facecolors='none', edgecolors='r',label='AS: %d' % sum(selAS60))
    ax.scatter(-scaleT(-Time[selFSlong],30,60,ylong), clust['Lat'][selFSlong].values,s=25, edgecolors='b', facecolors='none')
    ax.scatter(scaleT(Time[selASlong],30,60,ylong), clust['Lat'][selASlong].values,s=25, edgecolors='r', facecolors='none')
    # ax.scatter(np.log10(Time[selFS]), clust['Lat'][selFS].values,s=25, facecolors='none', edgecolors='b',label='FS: %d' % sum(np.logical_and(Time < 0, Time > -60)))
    # ax.scatter(-np.log10(Time[selAS]), clust['Lat'][selAS].values,s=25, facecolors='none', edgecolors='r',label='AS: %d' % sum(np.logical_and(Time > 0, Time < 60)))
    ax.scatter(0, Lat0, 100, marker='*', color='m',alpha=1.0)
    ax.plot([30,30],Lats,'--k')
    ax.plot([-30,-30],Lats,'--k')
    ax.plot([60,60],Lats,'-k')
    ax.plot([-60,-60],Lats,'-k')
    ax.set_xlim([-90, 90])
    ax.set_xticks([-90,-60,-30,0,30,60,90])
    ax.set_xticklabels(['-%d years' % years,'-60d','-30d','ot','30d','60d','%d years' % years])
    ax.set_ylim(Lats)
    ax.legend(loc='upper left')
    ax.set_xlabel('Time from t0 (days)')
    ax.set_ylabel('Latitude')



def plotR_TimeClust(ax, clust, ot,Lat0, Lon0, M0, cat, Lats, Lons, years):
    # years = 1
    ylong = years*365 # years
    t0 = UTCDateTime(ot)
    T = (cat.ot - t0)/(3600*24)
    # catalog time
    Ilon = np.logical_and(cat.Long > min(Lons), cat.Long < max(Lons))
    Ilat = np.logical_and(cat.Lat > min(Lats), cat.Lat < max(Lats))
    Ireg = np.logical_and(Ilat, Ilon)
    It60 = np.logical_and(T > -60, T < 60)
    # It60 = np.logical_and(cat.ot > t0 - 60*3600*24, cat.ot < t0 + 60*3600*24)
    # ItlongAS = np.logical_and(cat.ot > t0 + 60*3600*24, cat.ot < t0 + ylong*3600*24)
    # ItlongFS = np.logical_and(cat.ot < t0 - 60*3600*24, cat.ot > t0 - ylong*3600*24)
    ItlongAS = np.logical_and(T > 60, T < ylong)
    ItlongFS = np.logical_and(T < -60, T > -ylong)

    # R = DistLatLon(Lat0, Lon0, cat.Lat, cat.Long) * 1000
    R = DistLatLonUTM(Lat0, Lon0, cat.Lat, cat.Long)
    # Rc = DistLatLon(Lat0, Lon0, clust.Lat.values, clust.Lon.values) * 1000
    Rc = DistLatLonUTM(Lat0, Lon0, clust.Lat.values, clust.Long.values)
    [r_predicted, _] = WnCfaultL(M0,ifault=0)
    fact_r = AdjustEffectiveRadius(M0)
    r_predicted = np.round(fact_r * r_predicted/1000) * 1000

    # cluster time
    Time = []
    for ii in range(len(clust['ot'])):
        Time.append(UTCDateTime(clust['ot'].values[ii]))
    Time = np.array(Time)
    Time = (Time - t0) / (3600*24)

    selFS60 = np.logical_and(Time < 0, Time > -60)
    selFSlong = np.logical_and(Time <= -60, Time > -ylong)
    selAS60 = np.logical_and(Time > 0, Time < 60)
    selASlong = np.logical_and(Time >= 60, Time < ylong)

    ax.scatter(T[np.logical_and(It60, Ireg)], R[np.logical_and(It60, Ireg)], 15, marker='o', color='gray', alpha=0.6)
    ax.scatter(scaleT(T[np.logical_and(ItlongAS, Ireg)], 30, 60, ylong), R[np.logical_and(ItlongAS,Ireg)], 15, marker='o', color='gray', alpha=0.6)
    ax.scatter(-scaleT(-T[np.logical_and(ItlongFS, Ireg)], 30, 60, ylong), R[np.logical_and(ItlongFS,Ireg)], 15, marker='o', color='gray', alpha=0.6)

    ax.scatter(Time[selFS60], Rc[selFS60],s=25, facecolors='none', edgecolors='b',label='FS: %d' % (sum(selFS60)+sum(selFSlong)))
    ax.scatter(Time[selAS60], Rc[selAS60],s=25, facecolors='none', edgecolors='r',label='AS: %d' % (sum(selAS60)+sum(selASlong)))
    ax.scatter(-scaleT(-Time[selFSlong],30,60,ylong), Rc[selFSlong],s=25, edgecolors='b', facecolors='none')
    ax.scatter(scaleT(Time[selASlong],30,60,ylong), Rc[selASlong],s=25, edgecolors='r', facecolors='none')

    ax.plot([0,0], [min(R), max(R)], '--m')
    ax.plot([30,30], [min(R), max(R)], '--k')
    ax.plot([-30,-30], [min(R), max(R)], '--k')
    ax.plot([60,60], [min(R), max(R)], '-k')
    ax.plot([-60,-60], [min(R), max(R)], '-k')
    ax.plot([-90,90], [r_predicted, r_predicted], '-m', linewidth=1)
    ax.set_xlim([-90, 90])
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_xticklabels(['-%d years' % years, '-60d', '-30d', 'ot', '30d', '60d', '%d years' % years])

    ax.legend(loc='upper left')
    ax.set_xlabel('Time from t0 (days)')
    ax.set_ylabel('Distance from mainshock* (m)')
    ax.set_yscale('log')
    ax.set_ylim([min(R[R>0]), max(R)])
    ax.set_yscale('log')


def stemW(ax, x, y, c, s, lbl):
    for ii in range(len(x)):
        ax.plot([x[ii], x[ii]], [0,y[ii]], c='grey', alpha=0.8)
    ax.scatter(x, y, edgecolors=c, facecolors='w', marker=s, label=lbl)


def plotTimesMclust(ax, clust, ot, M0, cat,Lats, Lons, years):
    # years = 1
    ylong = years*365 # years
    t0 = UTCDateTime(ot)
    T = (cat.ot - t0)/(3600*24)
    # # catalog time
    Ilon = np.logical_and(cat.Long>min(Lons), cat.Long<max(Lons))
    Ilat = np.logical_and(cat.Lat>min(Lats), cat.Lat<max(Lats))
    Ireg = np.logical_and(Ilat, Ilon)
    It60 = np.logical_and(cat.ot > t0- 60*3600*24, cat.ot < t0 + 60*3600*24)
    ItlongAS = np.logical_and(cat.ot > t0 + 60*3600*24, cat.ot < t0 + ylong*3600*24)
    ItlongFS = np.logical_and(cat.ot < t0 - 60*3600*24, cat.ot > t0 - ylong*3600*24)

    # cluster time
    Time = []
    for ii in range(len(clust['ot'])):
        Time.append(UTCDateTime(clust['ot'].values[ii]))
    Time = np.array(Time)
    Time = (Time - t0) / (3600*24)

    selFS60 = np.logical_and(Time < 0, Time > -60)
    selFSlong = np.logical_and(Time <= -60, Time > -ylong)
    selAS60 = np.logical_and(Time > 0, Time < 60)
    selASlong = np.logical_and(Time >= 60, Time < ylong)

    stemW(ax, Time[selFS60], clust['M'][selFS60].values, 'b', 'o', 'FS: %d' % sum(selFS60))
    stemW(ax, Time[selAS60], clust['M'][selAS60].values, 'r', 'o', 'AS: %d' % sum(selAS60))
    stemW(ax,-scaleT(-Time[selFSlong],30,60,ylong), clust['M'][selFSlong].values,'b','o','')
    stemW(ax, scaleT(Time[selASlong],30,60,ylong), clust['M'][selASlong].values,'r','o','')

    ax.scatter(T[np.logical_and(It60, Ireg)], cat.M[np.logical_and(It60,Ireg)], facecolors='grey', edgecolors='none', s=25, alpha=0.5)
    ax.scatter(scaleT(T[np.logical_and(ItlongAS, Ireg)],30,60,ylong), cat.M[np.logical_and(ItlongAS,Ireg)], facecolors='grey', edgecolors='none', s=25, alpha=0.5)
    ax.scatter(-scaleT(-T[np.logical_and(ItlongFS, Ireg)],30,60,ylong), cat.M[np.logical_and(ItlongFS,Ireg)], facecolors='grey', edgecolors='none', s=25, alpha=0.5)

    stemW(ax,[0], [M0],'m','*','')
    ax.plot([30,30],[0, M0+1],'--k')
    ax.plot([-30,-30],[0, M0+1],'--k')
    ax.plot([60,60],[0, M0+1],'-k')
    ax.plot([-60,-60],[0, M0+1],'-k')
    ax.set_xlim([-90, 90])
    # ax.set_xlim([-10, 10])
    ax.set_xticks([-90,-60,-30,0,30,60,90])
    ax.set_xticklabels(['-%d years' % years,'-60d','-30d','ot','30d','60d','%d years' % years])
    # ax.legend(loc='upper left')
    ax.set_xlabel('Time from t0')
    ax.set_ylabel('Magnitude')
    ax.set_ylim([min(cat.M)-1.0, M0+1])


def plotworldloc2(ax, Lon0, Lat0,Lats,Longs):
    m = Basemap(projection='merc', llcrnrlat=min(Lats), urcrnrlat=max(Lats), llcrnrlon=min(Longs), urcrnrlon=max(Longs), lat_ts=1, resolution='i')
    m.drawcoastlines(color='k', linewidth=0.2)
    m.fillcontinents(color='#cc9955', alpha=0.4, zorder=1)
    x0, y0 = m(Lon0, Lat0)
    ax.scatter(x0, y0, 50, marker='*', color='m')


def plot_sequence(cat0):
    CLUST = cat0.CLUST
    deCLUST = cat0.deCLUST
    yearsLong = 2
    files = glob.glob('figs/*')
    for f in files:
        os.remove(f)
    for ii in range(len(CLUST['m0'])):

        fig = plb.figure(10000 + ii, figsize=(8, 8), dpi=150)
        clust = deCLUST[deCLUST['Cid'] == CLUST['cid'][ii]]
        print('Mw%s %s' % (CLUST['m0'][ii], UTCDateTime(CLUST['ot0'][ii]).date))

        [r_predicted,_] = WnCfaultL(CLUST['m0'][ii], ifault=0)
        fact_r = AdjustEffectiveRadius(CLUST['m0'][ii])
        r_predicted = fact_r * r_predicted
        Ro = kilometers2degrees(r_predicted / 1000)
        difx = Ro*2.1

        Lats0 = CLUST['lat0'][ii]
        Lons0 = CLUST['lon0'][ii]

        Lats = [Lats0-difx, Lats0+difx]
        Lons = [Lons0-difx, Lons0+difx]

        ax1a = fig.add_subplot(1, 2, 1)
        ax1b = fig.add_subplot(2, 2, 2)
        ax1c = fig.add_subplot(2, 2, 4)
        plotMapC(ax1a, clust, Lats, Lons, CLUST['lon0'][ii], CLUST['lat0'][ii], CLUST['ot0'][ii], CLUST['m0'][ii], deCLUST, yearsLong)
        plotTimesMclust(ax1b, clust, CLUST['ot0'][ii], CLUST['m0'][ii], deCLUST, Lats, Lons, yearsLong)
        plotR_TimeClust(ax1c, clust, CLUST['ot0'][ii], CLUST['lat0'][ii], CLUST['lon0'][ii], CLUST['m0'][ii], deCLUST, Lats, Lons, yearsLong)
        set_box_aspect(ax1b, 0.5)
        # set_box_aspect(ax1c, 0.5)
        ax1b.set_title(id2ctype(CLUST['c_type'][ii]))
        plb.suptitle('%s Mw%s' % (UTCDateTime(CLUST['ot0'][ii]).date, CLUST['m0'][ii]))

        fname = '%sOT_%s_Mw%s.pdf' % ('figs/', UTCDateTime(CLUST['ot0'][ii]).date, CLUST['m0'][ii])
        fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.2) # pad=0.4, w_pad=0.5, h_pad=0.9
        ax10 = fig.add_subplot(3, 6, 2)
        Lons0 = [min(deCLUST.Long), max(deCLUST.Long)]
        Lats0 = [min(deCLUST.Lat), max(deCLUST.Lat)]
        plotworldloc2(ax10, CLUST['lon0'][ii], CLUST['lat0'][ii], Lats0, Lons0)
        fig.savefig(fname)
        plb.close(fig)


def map_seis(cat0):
    fig = plb.figure(500)
    ax = fig.add_subplot(1, 2, 1)
    Msel = np.max(cat0.M)-1.0
    sel7 = cat0.M >=Msel

    ax.plot(cat0.Long, cat0.Lat, 'ko', ms=1, alpha=0.2)
    ax.plot( cat0.Long[sel7], cat0.Lat[sel7], 'ro', ms=8, mew=1.5, mfc='none', label='M>%2.1f' % Msel)
    ax.set_ylim([cat0.Lat.min(), cat0.Lat.max()])
    ax.set_xlim([cat0.Long.min(), cat0.Long.max()])
    ax.osm = OSM(ax)

    legend1 = ax.legend(loc='upper left')
    legend1._legend_title_box._children[0]._fontproperties._family = 'Impact'
    legend1._legend_title_box._children[0]._fontproperties.set_size(14)
    legend1.set_title(cat0.name)

    ax2 = fig.add_subplot(1, 2, 2)
    I_bg = cat0.c == 0
    ax2.plot(cat0.Long[I_bg], cat0.Lat[I_bg], 'ko', ms=1, alpha=0.5, label='Background')
    ax2.plot(cat0.Long[~I_bg], cat0.Lat[~I_bg], 'ro', ms=1, alpha=0.5, label='Clustered')
    ax2.set_ylim([cat0.Lat.min(), cat0.Lat.max()])
    ax2.set_xlim([cat0.Long.min(), cat0.Long.max()])
    ax2.osm = OSM(ax2)
    legend1 = ax2.legend(loc='upper left')
    legend1._legend_title_box._children[0]._fontproperties._family = 'Impact'
    legend1._legend_title_box._children[0]._fontproperties.set_size(14)
    legend1.set_title(cat0.name)
    # time.sleep(5)
    fig.tight_layout()
    time.sleep(5)
    fig.savefig('figs/map.pdf')



