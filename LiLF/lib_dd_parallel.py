import os, sys
import numpy as np
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

def merge_faintest_patch(skymodel, applyBeam):
    fluxes = skymodel.getColValues('I', aggregate='sum', applyBeam=applyBeam)
    names =skymodel.getPatchNames()
    faintest = names[np.argmin(fluxes)]
    pos = skymodel.getPatchPositions()[faintest]
    distances = skymodel.getDistance(*pos, byPatch=True)
    closest = names[np.argsort(distances)[1]] # 0 would be the faintest patch itself
    skymodel.merge([closest, faintest])
    skymodel.setPatchPositions(method='wmean', applyBeam=applyBeam)
    return skymodel


def merge_faint_facets(skymodel, min_flux, applyBeam=False):
    """
    Merge all patches of a skymodel with a flux density below a threshold with the closest patch.

    Parameters
    ----------
    skymodel
    min_flux
    applyBeam

    Returns
    -------
    merged skymodel
    """
    skymodel = skymodel.copy()
    i = 0
    while min(skymodel.getColValues('I', aggregate='sum', applyBeam=applyBeam)) < min_flux:
        print(i)
        merge_faintest_patch(skymodel, applyBeam=applyBeam)
        i += 1
    return skymodel

def closest_distance_between_patches(skymodel):
    """
    Return the name of the two patches which are closes and their distance

    Parameters
    ----------
    skymodel

    Returns
    -------
    closest_name: (str, str) - names of the two closest patches
    closest_distance: float - distance
    """

    closest_patch = np.zeros(len(skymodel.getPatchNames()))
    closest_name = []
    names = skymodel.getPatchNames()
    for i, (name, pos) in enumerate(skymodel.getPatchPositions().items()):
        distances = skymodel.getDistance(*pos, byPatch=True)
        closest_patch[i] = np.sort(distances)[1]
        nearby_name = names[distances == closest_patch[i]]
        # Case multiple patches at distance = 0
        if len(nearby_name) > 1:
            logger.warning(f'Possible issue, multiple patches with same distance {nearby_name}!')
            if nearby_name[0] == name:
                nearby_name = nearby_name[1]
            else:
                nearby_name = nearby_name[0]
        else: nearby_name = nearby_name[0]

        closest_name.append([name,nearby_name])
    return closest_name[np.argmin(closest_patch)], np.min(closest_patch)


def merge_nearby_bright_facets(skymodel, max_distance, min_flux, applyBeam=False):
    """
    Merge all bright patches of a skymodel that are within min_distance of another patch

    Parameters
    ----------
    skymodel
    max_distance: max distance to merge
    min_flux: min flux of facets to be considered bright
    applyBeam

    Returns
    -------
    merged skymodel
    """
    skymodel = skymodel.copy()
    skymodel_bright = skymodel.copy()
    skymodel_bright.select(f'I>{min_flux}', aggregate='sum', applyBeam=applyBeam)
    if len(skymodel_bright) > 1:
        while closest_distance_between_patches(skymodel_bright)[1] < max_distance:
            closest_patches, closest_distance = closest_distance_between_patches(skymodel_bright)
            logger.info(f'Merging nearby bright patches {closest_patches[0]} {closest_patches[1]} (distance={closest_distance*3600:.2f}arcsec')
            skymodel_bright.merge(closest_patches)
            skymodel.merge(closest_patches)
    else:
        logger.warning(f'Only one bright source - nothing to merge.')
    skymodel.setPatchPositions(method='wmean', applyBeam=applyBeam)
    return skymodel


def rename_skymodel_patches(skymodel, applyBeam=False):
    """
    Rename the patches in the input sky model according to flux

    Parameters
    ----------
    skymodel : LSMTool skymodel.SkyModel object
        Input sky model
    applyBeam : bool, intrinsic/apparent
    """
    if not skymodel.hasPatches:
        raise ValueError('Cannot rename patches since the input skymodel is not grouped '
                         'into patches.')
    patch_names = skymodel.getPatchNames()
    patch_fluxes = skymodel.getColValues('I', aggregate='sum', applyBeam=applyBeam)
    patch_pos = skymodel.getPatchPositions()

    old_new_dict = {}
    for i, id in enumerate(np.argsort(patch_fluxes)[::-1]):
        old_new_dict[patch_names[id]] = f'patch_{i:02.0f}'

    patch_col = skymodel.getColValues('Patch')
    for old_name, new_name in old_new_dict.items():
        patch_col[patch_col == old_name] = new_name
        patch_pos[new_name] = patch_pos.pop(old_name)

    skymodel.setColValues('Patch', patch_col)
    skymodel.setPatchPositions(patch_pos)




# class Direction(object):

#     def __init__(self, name):
#         self.name = name 
#         self.isl_num = int(name.split('_')[-1])
#         self.mask_voro = None
#         self.position_facet = None # [deg, deg]
#         self.position_cal = None # [deg, deg]
#         self.flux_cal = None # Jy
#         self.flux_facet = None # Jy
#         self.region_facet = None
#         self.size_cal = None # [deg, deg]
#         self.size_facet = None # [deg, deg]
#         self.cal_has_facet = None # Bool that tells if the cal is within the mask_voro
#         # lib_img.Image objects:
#         self.image = None
#         self.image_res = None
#         self.image_low = None
#         self.image_high = None

#         self.h5parms = {}
#         self.skymodel = None
#         self.skydb = None

#     def is_in_beam(self):
#         """
#         Return true if the direction is in the beam or an outsider
#         """
#         pass

#     def set_position(self, position, cal=True):
#         """
#         cal: if Tue is position of the central calibrator otherwise central position of the facet.
#         """
#         if cal:
#             self.position_cal = position
#         else:
#             self.position_facet = position

#     def set_flux(self, flux, cal=True, freq='mid'):
#         """
#         cal: if Tue is flux of the central calibrator otherwise of the whole facet.
#         freq: 'mid' frequency or 'min' frequency
#         """
#         if cal:
#             if freq == 'mid':
#                 self.flux_cal = flux
#             elif freq == 'min':
#                 self.flux_cal_minfreq = flux
#         else:
#             self.flux_facet = flux

#     def set_size(self, size, cal=True):
#         """
#         size: [deg,deg]
#         """
#         if cal:
#             self.size_cal = size
#         else:
#             self.size_facet = size

#     def add_mask_voro(self, mask_voro):
#         """
#         """
#         # read mask
#         fits = pyfits.open(mask_voro)
#         hdr, data = lib_img.flatten(fits)
#         w = pywcs.WCS(hdr)
#         pixsize_ra = -1*hdr['CDELT1']
#         pixsize_dec = hdr['CDELT2']

#         coord = np.where(data.T == self.isl_num)
#         if len(coord[0]) == 0:
#             self.cal_has_facet = False
#             self.size = [0.1,0.1]
#             self.position_facet = self.position_cal
#         else:
#             self.cal_has_facet = True
#             # calculate size
#             size_ra = (np.max(coord[0])-np.min(coord[0]))*pixsize_ra
#             size_dec = (np.max(coord[1])-np.min(coord[1]))*pixsize_dec
#             self.size_facet = [size_ra, size_dec]
#             # calculate position 
#             dir_x = np.mean([ np.max(coord[0]), np.min(coord[0]) ])
#             dir_y = np.mean([ np.max(coord[1]), np.min(coord[1]) ])
#             ra, dec =  w.all_pix2world(dir_x, dir_y, 0, ra_dec_order=True)
#             self.position_facet = [float(ra), float(dec)]


def make_voronoi_reg(directions, fitsfile, outdir_reg='regions', out_mask=None, png=None):
    """
    Take a list of coordinates and an image and voronoi tesselate the sky.
    It saves ds9 regions + fits mask of the facets

    directions : array of Direction objects
    firsfile : mask fits file to tassellate (used for coordinates and as template for the out_mask)
    outdir_reg : dir where to save regions
    out_mask : output mask with different numbers in each facet
    png : output png file that shows the tassellation
    """

    def closest_node(node, nodes):
        """
        Return closest values to node from nodes
        """
        nodes = np.asarray(nodes)
        dist_2 = np.sum((nodes - node)**2, axis=1)
        return np.argmin(dist_2)

    logger.debug("Image used for tasselation reference: "+fitsfile)
    fits = pyfits.open(fitsfile)
    hdr, data = lib_img.flatten(fits)
    w = pywcs.WCS(hdr)
    pixsize = np.abs(hdr['CDELT1'])

    # Get facets central pixels
    ras = np.array([d.position_cal[0] for d in directions])
    decs = np.array([d.position_cal[1] for d in directions])
    x_fs, y_fs = w.all_world2pix(ras, decs, 0, ra_dec_order=True)
    # keep trak of numbers in the direction names to name correctly patches in the fits files
    # in this way Isl_patch_12 will have "12" into the fits for that patch.
    nums = [d.isl_num for d in directions]

    x_c = data.shape[0]/2.
    y_c = data.shape[1]/2.

    # Check if dir is in img, otherwise drop
    idx_for_facet = []
    for i, direction in enumerate(directions):
        x, y = w.all_world2pix(ras[i], decs[i], 0, ra_dec_order=True)
        if x < 0 or x > data.shape[0] or y < 0 or y > data.shape[1]:
            logger.info('Direction %s is outside the fitsfile and will not have a facet.' % direction.name)
        else:
            idx_for_facet.append(i)

    # convert to pixel space (voronoi must be in eucledian space)
    x1 = 0
    y1 = 0
    x2 = data.shape[1] # note that y is before x in fits.data
    y2 = data.shape[0]

    # do tasselization
    vor = Voronoi(np.array((x_fs[idx_for_facet], y_fs[idx_for_facet])).transpose())
    box = np.array([[x1,y1],[x2,y2]])
    impoly = voronoi_finite_polygons_2d_box(vor, box)

    # create fits mask (each region one number)
    x, y = np.meshgrid(np.arange(x2), np.arange(y2)) # make a canvas with coordinates
    x, y = x.flatten(), y.flatten()
    pixels = np.vstack((x,y)).T 
    data_facet = np.zeros(shape=data.shape)
    for num, poly in zip(nums,impoly):
        p = Path(poly)
        pixels_region = p.contains_points(pixels)
        data_facet[ pixels_region.reshape(y2,x2) ] = num

    # put all values in each island equal to the closest region
    struct = generate_binary_structure(2, 2)
    data = binary_dilation(data, structure=struct, iterations=3).astype(data.dtype) # expand masks
    blobs, number_of_blobs = label(data.astype(int).squeeze(), structure=[[1,1,1],[1,1,1],[1,1,1]])
    center_of_masses = center_of_mass(data, blobs, list(range(number_of_blobs+1)))
    for blob in range(1,number_of_blobs+1):
        # get closer facet
        facet_num = closest_node(center_of_masses[blob], np.array([y_fs[idx_for_facet],x_fs[idx_for_facet]]).T)
        # put all pixel of that mask to that facet value
        data_facet[ blobs == blob ] = nums[facet_num]

    # save regions
    if not os.path.isdir(outdir_reg): os.makedirs(outdir_reg)

    all_s = []
    for i, poly in enumerate(impoly):
        ra, dec = w.all_pix2world(poly[:,0],poly[:,1], 0, ra_dec_order=True)
        coords = np.array([ra,dec]).T.flatten()

        s = Shape('Polygon', None)
        s.coord_format = 'fk5'
        s.coord_list = coords # ra, dec, radius
        s.coord_format = 'fk5'
        s.attr = ([], {'width': '2', 'point': 'cross',
                       'font': '"helvetica 16 normal roman"'})
        s.comment = 'color=red'
        all_s.append(s)

        regions = pyregion.ShapeList([s])
        regionfile = outdir_reg+'/'+directions[idx_for_facet[i]].name+'.reg'
        regions.write(regionfile)

    # add names for all.reg
    for d in directions:
        s = Shape('circle', None)
        s.coord_format = 'fk5'
        s.coord_list = [ d.position_cal[0], d.position_cal[1], 0.01 ] # ra, dec, radius
        s.coord_format = 'fk5'
        s.attr = ([], {'width': '1', 'point': 'cross',
                       'font': '"helvetica 16 normal roman"'})
        s.comment = 'color=white text="%s"' % d.name
        all_s.append(s)

    regions = pyregion.ShapeList(all_s)
    regionfile = outdir_reg+'/all.reg'
    regions.write(regionfile)
    logger.debug('There are %i regions within the PB and %i outside (no facet).' % (len(idx_for_facet), len(directions) - len(idx_for_facet)))

    # save fits mask
    if out_mask is not None:
        pyfits.writeto(out_mask, data_facet, hdr, overwrite=True)

    # plot tesselization
    if png is not None:
        import matplotlib.pyplot as pl
        pl.figure(figsize=(8,8))
        ax1 = pl.gca()
        voronoi_plot_2d(vor, ax1, show_vertices=True, line_colors='black', line_width=2, point_size=4)
        for i, d in enumerate(directions): ax1.text(x_fs[i], y_fs[i], d.name, fontsize=15)
        ax1.plot([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1])
        ax1.set_xlabel('RA (pixel)')
        ax1.set_ylabel('Dec (pixel)')
        ax1.set_xlim(x1,x2)
        ax1.set_ylim(y1,y2)
        logger.debug('Save plot: %s' % png)
        pl.savefig(png)


def voronoi_finite_polygons_2d_box(vor, box):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    box.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    box : (2,2) float array
        corners of bounding box
        numpy.array([[x1,y1],[x2,y2]])
        example: [[0,0],[2000,2000]]

    Returns
    -------
    poly : array of M (N,2) arrays
        polygon coordinates for M revised Voronoi regions.

    """
    import matplotlib.transforms as mplTrans
    import matplotlib.path as mplPath

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")
    if box.shape != (2,2):
        raise ValueError("Bounding box should be 2x2 array ((x1,y1),(x2,y2))")

    radius = np.max(box)

    # define the bounding box transform from the box extent - to be used to intersect with the regions
    bbox = mplTrans.Bbox(box)

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Matrix of distances from corners to points
    corners = [(box[0,0],box[0,1]), (box[0,0], box[1,1]), (box[1,0], box[0,1]), (box[1,0], box[1,1]) ] # BL, TL, BR, TR
    dists = []
    for p in vor.points:
        dists.append([ np.sqrt( (c[0]-p[0])**2 + (c[1]-p[1])**2 ) for c in corners ])
    dists = np.array(dists)

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge
            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # add corners of the image as vertices to prevent corners from being excluded
        # before expand corners so the ordering wil be right
        corners = [(box[0,0]-radius,box[0,1]-radius), (box[0,0]-radius, box[1,1]+radius), (box[1,0]+radius, box[0,1]-radius), (box[1,0]+radius, box[1,1]+radius) ]
        for c, dist in enumerate(dists[p1]): # cycle on corners
            if dist == np.min(dists[:,c]):
                # add that corner to the region
                #print "add corner %i to region %i" % (c,p1)
                new_region.append(len(new_vertices))
                new_vertices.append(list(corners[c]))

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    regions, imvertices = new_regions, np.asarray(new_vertices)

    ## now force them to be in the bounding box
    poly = np.asarray([imvertices[v] for v in regions])

    newpoly = []

    for p in poly:
        polyPath = mplPath.Path(p)
        newpolyPath = polyPath.clip_to_bbox(bbox)
        # makes vertices on the edge of the image 1 pixel outside the image
        # this is to be sure to include the pixels on the edge
        for i, vertice in enumerate(newpolyPath.vertices):
            if vertice[0] == box[0,0]: newpolyPath.vertices[i][0] = box[0,0]-1
            if vertice[1] == box[0,1]: newpolyPath.vertices[i][1] = box[0,1]-1
            if vertice[0] == box[1,0]: newpolyPath.vertices[i][0] = box[1,0]+1
            if vertice[1] == box[1,1]: newpolyPath.vertices[i][1] = box[1,1]+1
        coords, idx = np.unique([ x+1j*y for (x,y) in newpolyPath.vertices], return_index=True) # using complex; a way for grouping
        coords = [[x.real,x.imag] for x in coords[np.argsort(idx)]] # preserve order
        coords += [coords[0]] # close the line
        newpoly.append( np.array(coords) )

    return np.asarray(newpoly)

# class Grouper( object ):
#     """
#     Based on: http://www.chioka.in/meanshift-algorithm-for-the-rest-of-us-python/
#     """

#     def __init__(self, coords, fluxes, kernel_size=0.2, look_distance=0.3, grouping_distance=0.03):
#         """
#         coords: x,y coordinates for source positions
#         fluxes: total flux for each source
#         kernel_size: attenuate attraction, it this the flux times a gaussian of the distance with this as sigma [deg]
#         look_distance: max distance to look for nearby sources [deg]
#         grouping_distance: [deg]
#         """
#         self.coords = np.array(coords)
#         self.fluxes = fluxes
#         self.kernel_size = kernel_size # deg
#         self.look_distance = look_distance # deg
#         self.grouping_distance = grouping_distance # deg orig: 0.01
#         self.past_coords = [np.copy(self.coords)]
#         self.n_iterations = 100
#         self.clusters = []
#         logger.debug("Grouper: kernel_size=%.1f; look_distance=%.1f; grouping_distance=%.2f" % (kernel_size,look_distance,grouping_distance) )

#     def euclid_distance(self, coord, coords):
#         """
#         Simple ditance from coord to all coords
#         """
#         return np.sqrt(np.sum((coord - coords)**2, axis=1))
    
#     def neighbourhood_points(self, centroid, coords, max_distance):
#         """
#         Find close points, this reduces the load
#         """
#         distances = self.euclid_distance(centroid, coords)
#         #print('Evaluating: [%s vs %s] yield dist=%.2f' % (x, x_centroid, distance_between))
#         return np.where(distances < max_distance)
    
#     def gaussian_kernel(self, distance):
#         """
#         """
#         return (1/(self.kernel_size*np.sqrt(2*np.pi))) * np.exp(-0.5*((distance / self.kernel_size))**2)
    
#     def run(self):
#         """
#         Run the algorithm
#         """

#         for it in range(self.n_iterations):
#             logger.info("Grouper: Starting iteration %i" % it)
#             for i, x in enumerate(self.coords):
#                 ### Step 1. For each datapoint x in X, find the neighbouring points N(x) of x.
#                 idx_neighbours = self.neighbourhood_points(x, self.coords, max_distance = self.look_distance)
                
#                 ### Step 2. For each datapoint x in X, calculate the mean shift m(x).
#                 distances = self.euclid_distance(self.coords[idx_neighbours], x)
#                 weights = self.gaussian_kernel(distances)
#                 weights *= self.fluxes[idx_neighbours]**2 # multiply by flux**1.5 to make bright sources more important
#                 numerator = np.sum(weights[:,np.newaxis] * self.coords[idx_neighbours], axis=0)
#                 denominator = np.sum(weights)
#                 new_x = numerator / denominator
                
#                 ### Step 3. For each datapoint x in X, update x <- m(x).
#                 self.coords[i] = new_x

#             self.past_coords.append(np.copy(self.coords))

#             #if it>1: 
#             #    print (np.max(self.euclid_distance(self.coords,self.past_coords[-2])))

#             # if things changes little, brak
#             if it>1 and np.max(self.euclid_distance(self.coords, self.past_coords[-2])) < self.grouping_distance/2.: 
#                 break
            

#     def grouping(self):
#         """
#         Take the last coords set and group sources nearby, then return a list of lists. 
#         Each list has the index of one cluster.
#         """
#         coords_to_check = np.copy(self.coords)
#         while len(coords_to_check) > 0:
#             idx_cluster = self.neighbourhood_points(coords_to_check[0], self.coords, max_distance = self.grouping_distance)
#             idx_cluster_to_remove = self.neighbourhood_points(coords_to_check[0], coords_to_check, max_distance = self.grouping_distance)

#             # remove all coords of this clusters from the global list
#             mask = np.ones(coords_to_check.shape[0], dtype=bool)
#             mask[idx_cluster_to_remove] = False
#             coords_to_check = coords_to_check[mask]

#             # save this cluster indexes
#             self.clusters.append(idx_cluster)

#         logger.info('Grouper: Creating %i groups.' % len(self.clusters))
#         return self.clusters


#     def plot(self):
#         """
#         Plot the status of the distribution
#         """
#         import matplotlib as mpl
#         mpl.use("Agg")
#         import matplotlib.pyplot as plt
       
#         # decent colors
#         import cycler, random
#         color_idx = np.linspace(0, 1, len(self.clusters))
#         random.shuffle(color_idx)
#         color = plt.cm.rainbow(color_idx)
#         mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)

#         logger.info('Plotting grouped sources: grouping_xxx.png')
#         for i, X in enumerate(self.past_coords):
#             fig = plt.figure(figsize=(8, 8))
#             fig.subplots_adjust(wspace=0)
#             ax = fig.add_subplot(111)

#             initial_x = self.past_coords[0][:,0]
#             initial_y = self.past_coords[0][:,1]

#             ax.plot(initial_x,initial_y,'k.')
#             ax.plot(X[:,0],X[:,1],'ro')

#             ax.set_xlim( np.min(initial_x), np.max(initial_x) )
#             ax.set_ylim( np.min(initial_y), np.max(initial_y) )

#             #print ('Saving plot_%i.png' % i)
#             ax.set_xlim(ax.get_xlim()[::-1]) # reverse RA
#             fig.savefig('grouping_%00i.png' % i, bbox_inches='tight')

#         # plot clustering
#         fig = plt.figure(figsize=(8, 8))
#         fig.subplots_adjust(wspace=0)
#         ax = fig.add_subplot(111)
#         for cluster in self.clusters:
#             ax.plot(initial_x[cluster],initial_y[cluster], marker='.', linestyle='')

#         ax.set_xlim( np.min(initial_x), np.max(initial_x) )
#         ax.set_ylim( np.min(initial_y), np.max(initial_y) )
#         ax.set_xlim(ax.get_xlim()[::-1]) # reverse RA

#         logger.info('Plotting: grouping_clusters.png')
#         fig.savefig('grouping_clusters.png', bbox_inches='tight')
