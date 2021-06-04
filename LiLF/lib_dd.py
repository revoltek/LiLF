import os, sys, glob
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import pyregion
from pyregion.parser_helper import Shape
import lsmtool

from LiLF.lib_log import logger
from LiLF import lib_img, lib_util


class Direction(object):

    def __init__(self, name):
        self.name = name 
        self.position = None  # [deg, deg]
        self.size = None  # deg (1 value)
        self.localrms = None  # Jy/b (1 value)
        self.fluxes = None  # Jy - for each component
        self.spidx_coeffs = None  # 1st order - for each component
        self.ref_freq = None  # for each component
        self.converged = None  # bool
        self.peel_off = None
        self.region_file = None

        self.model = {}
        self.h5parms = {'ph':[],'fr':[],'amp1':[],'amp2':[]}

        # for debug
        self.avg_t = 0
        self.avg_f = 0
        self.rms_noise = []
        self.mm_ratio = []

    def clean(self):
        # TODO: remove files?

        # associated h5parms
        self.h5parms = {'ph':[],'fr':[],'amp1':[],'amp2':[]}

    def set_region(self, loc):
        """
        Creates a ds9 regionfile that covers the DD-cal model
        """
        assert self.size is not None and self.position is not None  # we need this to be already set

        self.region_file = loc+'/'+self.name+'.reg'
        s = Shape('circle', None)
        s.coord_format = 'fk5'
        s.coord_list = [ self.position[0], self.position[1], self.size/2. ]  # ra, dec, radius
        s.coord_format = 'fk5'
        s.attr = ([], {'width': '2', 'point': 'cross',
                       'font': '"helvetica 16 normal roman"'})
        if not self.peel_off: s.comment = 'color=red text="%s"' % self.name
        else: s.comment = 'color=blue text="%s"' % self.name

        regions = pyregion.ShapeList([s])
        lib_util.check_rm(self.region_file)
        regions.write(self.region_file)

    def get_region(self):
        if self.region_file is None:
            raise "Missing region file."

        return self.region_file

    def set_model(self, root, typ, apply_region=True):
        """
        apply_region: Isolate the clean components of a model fits file to those under self.region
        typ = init / best / regrid / ...
        """

        if apply_region:
            region_file = self.get_region()
            for model_file in glob.glob(root+'*[0-9]-model.fits'):
                lib_img.blank_image_reg(model_file, region_file, model_file, inverse=True, blankval=0.)

        self.model[typ] = root

    def get_model(self, typ):
        """
        Return the root name for the modelfile
        """
        try:
            return self.model[typ]
        except:
            raise "Model '%s' not set." % typ

    def add_h5parm(self, typ, h5parmFile):
        """
        typ can be 'ph', 'fr', 'amp1', or 'amp2'
        h5parmFile: filename
        """
        assert (typ == 'ph' or typ == 'fr' or typ == 'amp1' or typ == 'amp2')
        self.h5parms[typ].append(h5parmFile)

    def get_h5parm(self, typ, pos=-1):
        """
        typ can be 'ph', 'fr', 'amp1', or 'amp2'
        pos: the cycle from 0 to last added, negative numbers to search backwards, if non exists returns None
        """
        assert (typ == 'ph' or typ == 'fr' or typ == 'amp1' or typ == 'amp2')
        l = self.h5parms[typ]
        try:
            return l[pos]
        except:
            return None

    def add_rms_mm(self, rms_noise, mm_ratio):
        """
        track rms noise and mm ratio
        """
        self.rms_noise.append(rms_noise)
        self.mm_ratio.append(mm_ratio)

    def set_position(self, position, distance_peeloff, phase_center):
        """
        distance_peeloff: in deg, used to decide if the source is to peel_off
        phase_center: in deg
        """
        self.position = [round(position[0], 5), round(position[1], 5)]
        c1 = SkyCoord(position[0]*u.deg, position[1]*u.deg, frame='fk5')
        c2 = SkyCoord(phase_center[0]*u.deg, phase_center[1]*u.deg, frame='fk5')
        if c1.separation(c2).deg > distance_peeloff:
            self.peel_off = True
        else:
            self.peel_off = False

    def get_flux(self, freq):
        """
        freq: frequency to evaluate the flux
        """
        return np.sum(np.array(self.fluxes) * (freq/np.array(self.ref_freq))**(np.array(self.spidx_coeffs)))

    def set_size(self, ras, decs, majs, img_beam):
        """
        Calculate the size (diameter) of this calibrator measuring the distance of each component from the mean
        :param ras:  list of ras
        :param decs: list of decs
        :param majs: major axis sizes for sources [deg] - radius
        """
        ncomp = len(ras)
        if ncomp > 1:
            maxdist = 0
            center = SkyCoord(self.position[0]*u.deg, self.position[1]*u.deg, frame='fk5')
            for ra, dec, maj in zip(ras, decs, majs):
                comp = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
                dist = center.separation(comp).deg + maj
                if dist > maxdist:
                    maxdist = dist
            size = maxdist * 2
        else:
            size = majs[0]*2

        self.size = size * 1.2  # increase 20%

        if size < 3*img_beam:
            self.size = 3*img_beam
        #elif ncomp > 1 and size < 10*img_beam:
        #    # for complex sources force a larger region
        #    self.size = 8*img_beam


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
        return np.flatnonzero(distances < max_distance)
    
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
            if it > 1 and np.max(self.euclid_distance(self.coords, self.past_coords[-2])) < self.grouping_distance/2.:
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


    def merge_ids(self, ids):
        """ Merge groups containing ids"""
        if len(ids) < 2:
            return None
        clusters = self.clusters
        contains_id = []  # indices of clusters which contain one or more of the ids
        for id in ids:
            isin = [id in cluster for cluster in clusters]  # cluster_ids for clusters containing source ids
            contains_id.append(np.nonzero(isin))
        contains_id = np.unique(contains_id)
        if len(contains_id) == 1:  # all sources are already in the same cluster!
            return None
        else:
            merged = np.concatenate([clusters[id] for id in contains_id]) # this will be the merged cluster
            clusters = list(np.delete(clusters, contains_id)) # delete clusters that are merged so they don't apper twice
            logger.info('Merge groups in same mask island: {}'.format(merged))
            clusters.append(merged.astype(int))
            self.clusters = clusters


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
    lsm.select('Ra > %f' % (d.position[0]-(d.size/2)/np.cos(d.position[1]* np.pi / 180.)))
    lsm.select('Ra < %f' % (d.position[0]+(d.size/2)/np.cos(d.position[1]* np.pi / 180.)))
    lsm.select('Dec > %f' % (d.position[1]-d.size/2))
    lsm.select('Dec < %f' % (d.position[1]+d.size/2))
    if do_regions: lsm.write('ddcal/masks/regions-c%02i/%s.reg' % (cmaj,d.name), format='ds9', clobber=True)
    lsm.write(dir_skymodel, format='makesourcedb', clobber=True)
    lib_util.check_rm(dir_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (dir_skymodel, dir_skydb), log='makesourcedb_cl.log', commandType='general' )
    s.run(check=True)
