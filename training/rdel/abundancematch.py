from __future__ import print_function, division
from AbundanceMatching import AbundanceFunction, calc_number_densities, add_scatter, rematch, LF_SCATTER_MULT
import matplotlib as plt


def abundanceMatchSnapshot(proxy, scatter, lf, boxsize,
                            maxmag=-27., minmag=-5., debug=False,
                            figname=None):

    af = AbundanceFunction(lf[:,0], lf[:,1], (maxmag, minmag))

    #check the abundance function
    if debug==True:
        plt.semilogy(lf[:,0], lf[:,1])
        x = np.linspace(-27, -5, 101)
        plt.semilogy(x, af(x))
        plt.savefig(figname)

    #deconvolution and check results (it's a good idea to always check this)
    remainder = af.deconvolute(scatter*LF_SCATTER_MULT, 20)
    x, nd = af.get_number_density_table()
    plot(x, remainder/nd);

    #get number densities of the halo catalog
    nd_halos = calc_number_densities(proxy, box_size)

    #do abundance matching with some scatter
    catalog_sc = af.match(nd_halos, scatter*LF_SCATTER_MULT)

    return catalog_sc
