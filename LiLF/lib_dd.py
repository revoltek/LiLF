import os, sys, glob
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.convolution import convolve_fft
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import pyregion
from pyregion.parser_helper import Shape
import lsmtool

from LiLF.lib_log import logger
from LiLF import lib_img, lib_util


class Direction(object):

    def __init__(self, name, soltypes=['ph','ph1','ph2','fr','amp1','amp2']):
        self.name = name 
        self.soltypes = soltypes
        self.position = None  # [deg, deg]
        self.size = None  # deg (1 value)
        self.localrms = None  # Jy/b (1 value)
        self.fluxes = None  # Jy - for each component
        self.spidx_coeffs = None  # 1st order - for each component
        self.ref_freq = None  # for each component
        self.converged = None  # bool
        self.dist_from_centre = None # distance from the phase centre in deg
        self.peel_off = None
        self.region_file = None

        self.model = {}
        self.h5parms = {soltype:[] for soltype in self.soltypes}

        # for debug
        self.avg_t = 0
        self.avg_f = 0
        self.rms_noise = []
        self.mm_ratio = []

    def clean(self):
        # TODO: remove files?

        # associated h5parms
        self.h5parms = {soltype:[] for soltype in self.soltypes}

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
        else: s.comment = 'color=yellow text="%s"' % self.name

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
            for model_file in glob.glob(root+'*[0-9]-model.fits') + glob.glob(root+'*[0-9]-model-pb.fits'):
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
        assert typ in self.soltypes
        self.h5parms[typ].append(h5parmFile)

    def get_h5parm(self, typ, pos=-1):
        """
        pos: the cycle from 0 to last added, negative numbers to search backwards, if non exists returns None
        """
        assert typ in self.soltypes
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
        self.dist_from_centre = c1.separation(c2).deg
        if self.dist_from_centre > distance_peeloff:
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
        self.coords = np.array([SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='fk5') for ra,dec in coords])
        self.fluxes = fluxes
        self.kernel_size = kernel_size # deg
        self.look_distance = look_distance # deg
        self.grouping_distance = grouping_distance # deg orig: 0.01
        self.past_coords = [np.copy(self.coords)]
        self.n_iterations = 100
        self.clusters = []
        logger.debug("Grouper: kernel_size=%.1f; look_distance=%.1f; grouping_distance=%.2f" % (kernel_size,look_distance,grouping_distance) )

    def sky_distance(self, coord, coords):
        """
        Simple ditance from coord to all coords
        """
        #dist1 = np.sqrt(np.sum((coord - coords)**2, axis=1))
        catalogue = SkyCoord(ra = [c.ra for c in coords], dec = [c.dec for c in coords])
        dist = np.array([d.deg for d in coord.separation(catalogue)])
        return dist
    
    def neighbourhood_points(self, centroid, coords, max_distance):
        """
        Find close points, this reduces the load
        """
        distances = self.sky_distance(centroid, coords)
        #print('Evaluating: [%s vs %s] yield dist=%.2f' % (x, x_centroid, distance_between))
        return np.flatnonzero(distances < max_distance)
    
    def gaussian_kernel(self, distance):
        """
        """
        return (1/(self.kernel_size*np.sqrt(2*np.pi))) * np.exp(-0.5*((distance / self.kernel_size))**2)
    
    def weighted_spherical_mean(self, coords_list, weights):
        """
        Calculate the weighted spherical mean of a list of SkyCoord objects.

        Parameters:
        - coords_list: List of SkyCoord objects.
        - weights: List of weights corresponding to each coordinate.

        Returns:
        - Weighted spherical mean SkyCoord object.
        """
        # Convert coordinates to Cartesian representation
        cartesian_coords = np.array([coord.represent_as('cartesian').xyz.value for coord in coords_list])

        # Apply weights to Cartesian coordinates
        weighted_cartesian_coords = cartesian_coords * np.array(weights)[:, np.newaxis]

        # Calculate the weighted mean Cartesian coordinates
        weighted_mean_cartesian = np.sum(weighted_cartesian_coords, axis=0) / np.sum(weights)

        # Convert back to spherical coordinates
        mean_coord = SkyCoord(x=weighted_mean_cartesian[0], y=weighted_mean_cartesian[1], z=weighted_mean_cartesian[2], representation_type='cartesian').represent_as('spherical')

        return SkyCoord(ra=mean_coord.lon.deg*u.deg, dec=mean_coord.lat.deg*u.deg, frame='fk5')

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
                distances = self.sky_distance(x, self.coords[idx_neighbours])
                weights = self.gaussian_kernel(distances)
                weights *= self.fluxes[idx_neighbours]**2 # multiply by flux**1.5 to make bright sources more important
                new_x = self.weighted_spherical_mean(self.coords[idx_neighbours], weights)
                #print("xnewx",x,new_x)

                ### Step 3. For each datapoint x in X, update x <- m(x).
                self.coords[i] = new_x

            self.past_coords.append(np.copy(self.coords))

            if it > 1:
                dists = [self.coords[i].separation(self.past_coords[-2][i]).deg for i in range(len(self.coords))]
                if np.max(dists) < self.grouping_distance/2.:
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
            ax = fig.add_subplot(111, polar=True)

            initial_x = np.array([s.ra.wrap_at(180 * u.deg).rad for s in self.past_coords[0]])
            initial_y = np.array([np.pi/2 - s.dec.rad for s in self.past_coords[0]])
            ax.plot(initial_x,initial_y,'k.')
            ax.plot([s.ra.wrap_at(180 * u.deg).rad for s in X], [np.pi/2-s.dec.rad for s in X],'ro')

            ax.set_xlim( np.min(initial_x), np.max(initial_x) )
            ax.set_ylim( np.min(initial_y), np.max(initial_y) )
            ax.set_rorigin(-1*np.min(initial_y))

            fig.savefig('grouping_%00i.png' % i, bbox_inches='tight')

        # plot clustering
        fig = plt.figure(figsize=(8, 8))
        fig.subplots_adjust(wspace=0)
        ax = fig.add_subplot(111, polar=True)
        for cluster in self.clusters:
            ax.plot(initial_x[cluster],initial_y[cluster], marker='.', linestyle='')

        ax.set_xlim( np.min(initial_x), np.max(initial_x) )
        ax.set_ylim( np.min(initial_y), np.max(initial_y) )
        ax.set_rorigin(-1*np.min(initial_y))

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

def make_subfield_region(name, MS, sm, min_flux, debug_dir=None):
    """
    Identify the smallest region of sky model sources that contains a certain flux.

    Parameters
    ----------
    name: name of the ds9 region file to be created
    MS: lib_ms.MS object, MS to find center and size of field of view
    sm: lsmtool.skymodel, apparent skymodel of the field
    min_flux: float, minimum flux in calibration subfield in Jy
    debug_dir: str, default=None. Save debug fits files in this dir.

    Returns
    -------
    bestbox_coord: nparray (float), ra and dec of bestbox center in deg
    bestbox_size: float, size of bestbox in deg
    """
    c_ra, c_dec = MS.getPhaseCentre()
    fwhm = 2*MS.getFWHM(freq='mid')
    cellsize  = 1/60 # 1 arcmin
    # TODO padding?/offset?
    size_pix = int(np.round(1.1*fwhm/cellsize)) # size in pix, 10% pad
    freq = np.mean(MS.getFreqs())
    steps = 40 # how many scales to consider
    boxsizes = np.linspace(fwhm/(steps),fwhm,steps)

    # iterate over box sizes and convolve with boxcar kernel to find per scale the location of the box containing the
    # maximum flux
    smt = sm.table
    type, ra, dec, I, si = smt['Type'], smt['Ra'], smt['Dec'], smt['I'], smt['SpectralIndex']
    reff = sm.getColValues('ReferenceFrequency')
    sc = SkyCoord(ra, dec)
    # NOTE: ignore higher order SI terms for now...
    assert len(np.unique(reff)) == 1
    fluxes = I*(freq/reff[0])**si[:,0]

    def create_cone_kernel(size, edge_value=0.9):
        # helper function to create the kernel. The kernel is 1 at the center and radially drops to edge_value at the edge
        # Create a meshgrid to represent the indices of the kernel
        x, y = np.meshgrid(np.arange(size), np.arange(size))
        # Calculate the distance of each point from the center
        distance = np.sqrt((x - size // 2) ** 2 + (y - size // 2) ** 2)
        distance_factor = distance/np.max(distance)
        cell_value = 1 - (1-edge_value)*distance_factor

        return cell_value
    
    max_location, max_flux_in_field = np.zeros((len(boxsizes),2)), np.zeros_like(boxsizes)

    # iterate over box sizes, for each boxsize find the best region (center location = maximum of convolved map)
    for i,boxsize in enumerate(boxsizes):
        kernel_size = int(np.rint(boxsize/cellsize))
        if kernel_size % 2 == 0:
            kernel_size += 1
        # Kernel is not perfect 2D box, but slightly drops towards the edges to center bright sources.
        kernel = create_cone_kernel(kernel_size)
        # create template image object
        hdu = lib_util.get_template_image(c_ra, c_dec, size_pix, size_pix, cellsize)
        imdata = hdu[0].data
        w = WCS(hdu[0].header)
        pixcrd = np.rint(sc.to_pixel(w)).astype(int)
        if np.min(pixcrd) < 0 or np.max(pixcrd) > size_pix:
            logger.warning('There are skymodel sources outside the region considered for calibration!')
            #sys.exit()
        for flux, pix in zip(fluxes, pixcrd.T):
            imdata[pix[1],pix[0]] += flux

        hdu[0].data = convolve_fft(imdata,kernel,normalize_kernel=False)
        if debug_dir:
            hdu.writeto(f'{debug_dir}/flux_region_map_{int(np.rint(boxsize*60)):03}amin.fits', overwrite=True)
        max_location[i] = np.unravel_index(np.argmax(hdu[0].data), hdu[0].data.shape)
        max_flux_in_field[i] = np.max(hdu[0].data)
    # print(max_flux_in_field)
    mask = max_flux_in_field > min_flux # points above min flux
    id_min = np.argwhere(mask)[0,0]
    # find smallest region containing min flux
    bestbox_coord = w.wcs_pix2world([max_location[id_min]], 0)[0]
    bestbox_flux = max_flux_in_field[id_min]
    bestbox_size = boxsizes[id_min]
    logger.info(f'Flux {bestbox_flux:.1f}Jy > min_flux ({min_flux:.1f}Jy) in {bestbox_size:.2f}deg box at ra={bestbox_coord[0]:.2f}d dec={bestbox_coord[1]:.2f}d.')

    # search for points with strong jumps using gradient and boxsize^eidx
    # eidx - linear case: gradient in Jy per box diameter (finds larger boxes)
    #      - square case: gradient in Jy per box area ( finds smaller boxes)
    # eidx = 1.0
    # grad = (max_flux_in_field[1:] - max_flux_in_field[:-2]) / (boxsizes[1:]**eidx - boxsizes[:-2]**eidx)
    # print(grad)
    # maxgrad, boxsize_maxgrad, flux_maxgrad = np.max(grad), boxsizes[np.argmax(grad)+1], fluxes[np.argmax(grad)+1]
    # if fluxes[0] > min_flux:
    #     logger.info(' flux.')
    # if flux_maxgrad <= min_flux:
    #     logger.info('Box size with max gradient <= min flux. Use box containing min flux.')
    # else:

    #
    # fig, ax1 = plt.subplots()
    # ax1.plot(boxsizes, max_flux_in_field, marker='x')
    # ax1.hlines(min_flux, boxsizes[0], boxsizes[-1], label='minflux', c='k')
    # ax2 = ax1.twinx()
    # ax2.plot(boxsizes, np.gradient(max_flux_in_field, boxsizes), marker='x', c='C1')
    # ax1.set_xlabel('box size [deg]')
    # ax1.set_ylabel('flux [Jy]')
    # ax2.set_ylabel('gradient [Jy/deg]')
    # plt.legend()
    # plt.savefig('flux_size.png')
    # plt.close()
    # fig, ax1 = plt.subplots()
    # ax1.plot(boxsizes**2, max_flux_in_field, marker='x')
    # ax1.hlines(min_flux, (boxsizes**2)[0], (boxsizes**2)[-1], label='minflux', c='k')
    # ax2 = ax1.twinx()
    # ax2.plot(boxsizes**2, np.gradient(max_flux_in_field, boxsizes**2), marker='x', c='C1')
    # # ax2.plot(boxsizes**2, max_flux_in_field/boxsizes**2, marker='x')
    # ax1.set_xlabel('box area [degree^2]')
    # ax1.set_ylabel('flux [Jy]')
    # ax2.set_ylabel('gradient [Jy/deg^2]')
    # plt.legend()
    # plt.savefig('flux_area.png')
    region_string = f"""# Region file format: DS9 version 4.1
                        global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
                        fk5
                        box({bestbox_coord[0]},{bestbox_coord[1]},{1.02*bestbox_size},{1.02*bestbox_size},0.0)"""
    region = pyregion.parse(region_string)
    region.write(name)
    return bestbox_coord, bestbox_size



