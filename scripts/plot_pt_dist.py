import numpy as np
import matplotlib.pyplot as plt


def get_par_index(par):
    indexes = {"X":0, "Y":1, "Z":2, "Vx":3, "Vy":4, "Vz":5, "VC":6, "SC":7}
    return indexes[par]


def plot_1D_distribution(data, par, gbins, file_name):
    idx = get_par_index(par)
    plt.figure()
    plt.grid()
    plt.title("FILE-> %s  ;  PAR-> %s" % (file_name, par))
    plt.xlabel("%s values" % par)
    plt.ylabel("probability density")
    plt.hist(data[:,idx], bins=gbins, density=True, histtype='step', lw=1.5)
    plt.tight_layout()
    return 0


def plot_2D_distribution(data, pars, gbins, file_name):
    x_idx = get_par_index(pars[0])
    y_idx = get_par_index(pars[1])
    plt.figure()
    plt.grid()
    plt.title("FILE-> %s  ;  PAR-> %s  %s" % (file_name, pars[0], pars[1]))
    plt.xlabel("%s values" % pars[0])
    plt.ylabel("%s values" % pars[1])
    plt.hist2d(data[:,x_idx], data[:,y_idx], bins=gbins)
    plt.colorbar()
    plt.tight_layout()
    return 0


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", action="store", type=str, default="", help="file with particle statistics [def = '']")
    parser.add_argument("-p", "--parameters", action="store", type=str, default="", 
                        help="parameter names which distributions will be plotted [def='']\nPOSSIBLE NAMES: X Y Z Vx Vy Vz VC(volume_count) SC(surface_count)")
    parser.add_argument("-b", "--bins", action="store", type=int, default=50, help="number of bins along the axis [def=50]")
    parser.add_argument("-p2", "--dist_over_2_pars", action="store", default="", help="two parameters according to which you want to plot 2D hist, [def='']")
    args = parser.parse_args()
    
    if args.file=='':
        raise Warning("File with statistics was not set!")
        
    data = np.loadtxt(args.file)
    fname = args.file.split('/')[-1]
    if args.parameters!="":
        par_list = args.parameters.split(' ')
        for par in par_list:
            plot_1D_distribution(data, par, args.bins, fname)
    
    if args.dist_over_2_pars!='':
        pars = args.dist_over_2_pars.split(' ')
        plot_2D_distribution(data, pars, args.bins, fname)
    
    plt.show()
    
    