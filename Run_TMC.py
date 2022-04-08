from LoadData import read_SCEDC as read_SCEDC
from Plotters import map_seis, plot_sequence
import matplotlib.pyplot as plb
from TMC import TMC, MakaClusterData

years = [1985, 2020]  # i.e. 1985_01-01 to 2020-01-01
lats = [32.0, 37.0]
lons = [-120, -114]
max_depth = 70
Mc = 4.0  # catalog completness level
dM = 2.0  # mim mainshock mag = Mc + dM
file = 'SCEDC.csv'


def main():
    cat0 = read_SCEDC(file, lats, lons, Mc, years, max_depth)

    cat0 = TMC(cat0)

    cat0 = MakaClusterData(cat0, Mc, dM)

    plot_sequence(cat0)

    map_seis(cat0)

    plb.show()

if __name__ == "__main__":
    main()
