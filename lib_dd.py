import os, sys, itertools
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import fits as pyfits
from astropy import wcs as pywcs
import astropy.units as u
import pyregion
from pyregion.parser_helper import Shape
from matplotlib.path import Path
from scipy.ndimage import binary_dilation, generate_binary_structure
from scipy.ndimage.measurements import label, center_of_mass
try:
    from scipy.spatial import Voronoi, voronoi_plot_2d
except:
    logger.error("Load latest scipy with 'use Pythonlibs'")
    sys.exit(1)

from LiLF.lib_log import logger
from LiLF import lib_img

class Direction(object):

    def __init__(self, name):
        self.name = name 
        self.position = None # [deg, deg]
        self.size = None # [deg, deg]
        self.fluxes = None # Jy
        self.spidx_coeffs = None
        self.ref_freq = None

        # lib_img.Image objects:
        self.image = None
        self.image_res = None
        self.image_low = None
        self.image_high = None
        self.skymodel = None
        self.skydb = None

        # associated h5parm
        self.h5parms = {}


    def is_in_beam(self):
        """
        Return true if the direction is in the beam or an outsider
        """
        pass

    def set_position(self, position):
        """
        """
        self.position = position

    def set_flux(self, fluxes, spidx_coeffs=None, ref_freq=None, freq='mid'):
        """
        fluxes: list of flues of various components
        spidx_coeffs: spectral index coefficients to extract the flux for each component
        ref_freq: reference frequency (https://sourceforge.net/p/wsclean/wiki/ComponentList/)
        """
        self.fluxes = fluxes
        self.spidx_coeffs = spidx_coeffs
        self.ref_freq = ref_freq

    def get_flux(self, freq):
        """
        freq: frequency to evaluate the flux
        """
        fluxes = np.copy(self.fluxes)
        for i, term in enumerate(self.spidx_coeffs[0]):
            fluxes += self.spidx_coeffs[:,i] * ( freq/self.ref_freq - 1 )**(i+1)
        return np.sum(fluxes)

    def set_size(self, size):
        """
        size: [deg,deg]
        """
        self.size = size


class Grouper( object ):
    """
    Based on: http://www.chioka.in/meanshift-algorithm-for-the-rest-of-us-python/
    """

    def __init__(self, coords, fluxes, kernel_size=0.2, look_distance=0.3, grouping_distance=0.03):
        """
        coords: x,y coordinates for source positions
        fluxes: total flux for each source
        kernel_size: attenuate attraction, it this the flux times a gaussian of the distance with this as sigma [deg]
        look_distance: max distance to look for nearby sources [deg]
        grouping_distance: [deg]
        """
        self.coords = np.array(coords)
        self.fluxes = fluxes
        self.kernel_size = kernel_size # deg
        self.look_distance = look_distance # deg
        self.grouping_distance = grouping_distance # deg orig: 0.01
        self.past_coords = [np.copy(self.coords)]
        self.n_iterations = 100
        self.clusters = []
        logger.debug("Grouper: kernel_size=%.1f; look_distance=%.1f; grouping_distance=%.2f" % (kernel_size,look_distance,grouping_distance) )

    def euclid_distance(self, coord, coords):
        """
        Simple ditance from coord to all coords
        """
        return np.sqrt(np.sum((coord - coords)**2, axis=1))
    
    def neighbourhood_points(self, centroid, coords, max_distance):
        """
        Find close points, this reduces the load
        """
        distances = self.euclid_distance(centroid, coords)
        #print('Evaluating: [%s vs %s] yield dist=%.2f' % (x, x_centroid, distance_between))
        return np.where(distances < max_distance)
    
    def gaussian_kernel(self, distance):
        """
        """
        return (1/(self.kernel_size*np.sqrt(2*np.pi))) * np.exp(-0.5*((distance / self.kernel_size))**2)
    
    def run(self):
        """
        Run the algorithm
        """

        for it in range(self.n_iterations):
            logger.info("Grouper: Starting iteration %i" % it)
            for i, x in enumerate(self.coords):
                ### Step 1. For each datapoint x in X, find the neighbouring points N(x) of x.
                idx_neighbours = self.neighbourhood_points(x, self.coords, max_distance = self.look_distance)
                
                ### Step 2. For each datapoint x in X, calculate the mean shift m(x).
                distances = self.euclid_distance(self.coords[idx_neighbours], x)
                weights = self.gaussian_kernel(distances)
                weights *= self.fluxes[idx_neighbours]**2 # multiply by flux**1.5 to make bright sources more important
                numerator = np.sum(weights[:,np.newaxis] * self.coords[idx_neighbours], axis=0)
                denominator = np.sum(weights)
                new_x = numerator / denominator
                
                ### Step 3. For each datapoint x in X, update x <- m(x).
                self.coords[i] = new_x

            self.past_coords.append(np.copy(self.coords))

            #if it>1: 
            #    print (np.max(self.euclid_distance(self.coords,self.past_coords[-2])))

            # if things changes little, brak
            if it>1 and np.max(self.euclid_distance(self.coords, self.past_coords[-2])) < self.grouping_distance/2.: 
                break
            

    def grouping(self):
        """
        Take the last coords set and group sources nearby, then return a list of lists. 
        Each list has the index of one cluster.
        """
        coords_to_check = np.copy(self.coords)
        while len(coords_to_check) > 0:
            idx_cluster = self.neighbourhood_points(coords_to_check[0], self.coords, max_distance = self.grouping_distance)
            idx_cluster_to_remove = self.neighbourhood_points(coords_to_check[0], coords_to_check, max_distance = self.grouping_distance)

            # remove all coords of this clusters from the global list
            mask = np.ones(coords_to_check.shape[0], dtype=bool)
            mask[idx_cluster_to_remove] = False
            coords_to_check = coords_to_check[mask]

            # save this cluster indexes
            self.clusters.append(idx_cluster)

        logger.info('Grouper: Creating %i groups.' % len(self.clusters))
        return self.clusters


    def plot(self):
        """
        Plot the status of the distribution
        """
        import matplotlib as mpl
        mpl.use("Agg")
        import matplotlib.pyplot as plt
       
        # decent colors
        import cycler, random
        color_idx = np.linspace(0, 1, len(self.clusters))
        random.shuffle(color_idx)
        color = plt.cm.rainbow(color_idx)
        mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)

        logger.info('Plotting grouped sources: grouping_xxx.png')
        for i, X in enumerate(self.past_coords):
            fig = plt.figure(figsize=(8, 8))
            fig.subplots_adjust(wspace=0)
            ax = fig.add_subplot(111)

            initial_x = self.past_coords[0][:,0]
            initial_y = self.past_coords[0][:,1]

            ax.plot(initial_x,initial_y,'k.')
            ax.plot(X[:,0],X[:,1],'ro')

            ax.set_xlim( np.min(initial_x), np.max(initial_x) )
            ax.set_ylim( np.min(initial_y), np.max(initial_y) )

            #print ('Saving plot_%i.png' % i)
            ax.set_xlim(ax.get_xlim()[::-1]) # reverse RA
            fig.savefig('grouping_%00i.png' % i, bbox_inches='tight')

        # plot clustering
        fig = plt.figure(figsize=(8, 8))
        fig.subplots_adjust(wspace=0)
        ax = fig.add_subplot(111)
        for cluster in self.clusters:
            ax.plot(initial_x[cluster],initial_y[cluster], marker='.', linestyle='')

        ax.set_xlim( np.min(initial_x), np.max(initial_x) )
        ax.set_ylim( np.min(initial_y), np.max(initial_y) )
        ax.set_xlim(ax.get_xlim()[::-1]) # reverse RA

        logger.info('Plotting: grouping_clusters.png')
        fig.savefig('grouping_clusters.png', bbox_inches='tight')

def cut_skymodel(skymodel_in, skymodel_out, d, do_skydb=True, do_regions=False):
    """
    Load full skymodel and extract sources in the square around the calibrator of the given size
    """
    lsm = lsmtool.load(skymodel_in)
    # select all sources within a sqare of patch size
    lsm.select('Ra > %f' % (d.position[0]-(d.size[0]/2)/np.cos(d.position[1]* np.pi / 180.)))
    lsm.select('Ra < %f' % (d.position[0]+(d.size[0]/2)/np.cos(d.position[1]* np.pi / 180.)))
    lsm.select('Dec > %f' % (d.position[1]-d.size[1]/2))
    lsm.select('Dec < %f' % (d.position[1]+d.size[1]/2))
    if do_regions: lsm.write('ddcal/masks/regions-c%02i/%s.reg' % (cmaj,d.name), format='ds9', clobber=True)
    lsm.write(dir_skymodel, format='makesourcedb', clobber=True)
    lib_util.check_rm(dir_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (dir_skymodel, dir_skydb), log='makesourcedb_cl.log', commandType='general' )
    s.run(check=True)
