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
        self.isl_num = int(name[10:])
        self.mask_voro = None
        self.position_facet = None # [deg, deg]
        self.position_cal = None # [deg, deg]
        self.flux_cal = None # Jy
        self.flux_facet = None # Jy
        self.region_facet = None
        self.size = None # [deg, deg]
        self.cal_has_facet = None # Bool that tells if the cal is within the mask_voro
        # lib_img.Image objects:
        self.image = None
        self.image_res = None
        self.image_low = None
        self.image_high = None

    def is_in_bean(self):
        """
        Return true if the direction is in the beam or an outsider
        """
        return True

    def set_position(self, position, cal=True):
        """
        cal: if Tue is position of the central calibrator otherwise central position of the facet.
        """
        if cal:
            self.position_cal = position
        else:
            self.position_facet = position

    def set_flux(self, flux, cal=True):
        """
        cal: if Tue is flux of the central calibrator otherwise of the whole facet.
        """
        if cal:
            self.flux_cal = flux
        else:
            self.flux_facet = flux

    def add_mask_voro(self, mask_voro):
        """
        """
        # read mask
        fits = pyfits.open(mask_voro)
        hdr, data = lib_img.flatten(fits)
        w = pywcs.WCS(hdr)
        pixsize_ra = -1*hdr['CDELT1']
        pixsize_dec = hdr['CDELT2']

        coord = np.where(data.T == self.isl_num)
        if len(coord[0]) == 0:
            self.cal_has_facet = False
            self.size = [0.1,0.1]
            self.position_facet = self.position_cal
        else:
            self.cal_has_facet = True
            # calculate size
            size_ra = (np.max(coord[0])-np.min(coord[0]))*pixsize_ra
            size_dec = (np.max(coord[1])-np.min(coord[1]))*pixsize_dec
            self.size = [size_ra, size_dec]
            # calculate position 
            dir_x = np.mean([ np.max(coord[0]), np.min(coord[0]) ])
            dir_y = np.mean([ np.max(coord[1]), np.min(coord[1]) ])
            ra, dec =  w.all_pix2world(dir_x, dir_y, 0, ra_dec_order=True)
            self.position_facet = [float(ra), float(dec)]


#def split_directions(directions, fitsfile):
#    """
#    Return 2 direction dicts. One contains only the directions that are within the boundaries of the image.
#    The other the rest of the directions.
#    
#    directions : dict with {'Dir_0':[ra,dec], 'Dir_1':[ra,dec]...}
#    fitsfile: used to find what is in/out the image
#    """
#    fits = pyfits.open(fitsfile)
#    hdr, data = lib_img.flatten(fits)
#    w = pywcs.WCS(hdr)
#
#    directions_in = {}
#    directions_out = {}
#
#    for direction in directions:
#        # Get facets central pixels
#        ras = directions[direction][0].degree
#        decs = directions[direction][1].degree
#        x, y = w.all_world2pix(ras, decs, 0, ra_dec_order=True)
#        if x < 0 or x > data.shape[0] or y < 0 or y > data.shape[1]:
#            logger.info('Direction %s is outside the primary beam and will not have a facet (it will still be a calibrator).' % direction)
#            directions_out[direction] = directions[direction]
#        else:
#            directions_in[direction] = directions[direction]
#
#    return directions_in, directions_out
 

def make_voronoi_reg(directions, fitsfile, outdir_reg='regions', out_mask='facet.fits', png=None):
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
            logger.info('Direction %s is outside the primary beam and will not have a facet (it will still be a calibrator).' % direction.name)
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

    # save fits mask
    pyfits.writeto(out_mask, data_facet, hdr, overwrite=True)

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

    # plot tesselization
    if png is not None:
        import matplotlib.pyplot as pl
        pl.figure(figsize=(8,8))
        ax1 = pl.gca()
        voronoi_plot_2d(vor, ax1, show_vertices=True, line_colors='black', line_width=2, point_size=4)
        for i, d in enumerate(directions): ax1.text(x_fs[i], y_fs[i], d, fontsize=15)
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

#def sizes_from_mask_voro(mask_voro):
#    """
#    Compute image sizes from the mask_voro
#    mask_voro has different numbers for each patch
#
#    Returns
#    -------
#    Dict indexed by facet names with [size_ra and size_dec] in degrees
#    """
#    sizes = {}
#    # read mask
#    fits = pyfits.open(mask_voro)
#    hdr, data = lib_img.flatten(fits)
#    pixsize_ra = -1*hdr['CDELT1']
#    pixsize_dec = hdr['CDELT2']
#
#    # calculate sizes
#    for i in np.unique(data):
#        assert i != 0 # there should not be any 0 in the mask ptherwise something is wrong with tessellation
#        #print 'Working on', i
#        coord = np.where(data.T == i)
#        size_ra = (np.max(coord[0])-np.min(coord[0]))*pixsize_ra
#        size_dec = (np.max(coord[1])-np.min(coord[1]))*pixsize_dec
#        sizes['Isl_patch_%i' % i] = [size_ra, size_dec]
#        #print sizes[-1]
#
#    # return list
#    return sizes
#
#def directions_from_mask_voro(mask_voro):
#    """
#    Compute facet direction (centre) from the mask_voro
#    mask_voro has different numbers for each patch
#
#    Returns
#    -------
#    Dict indexed by facet names with direction [ra dec] in degrees
#    """
#    directions = {}
#    # read mask
#    fits = pyfits.open(mask_voro)
#    hdr, data = lib_img.flatten(fits)
#    w = pywcs.WCS(hdr)
#
#    # calculate sizes
#    for i in np.unique(data):
#        assert i != 0 # there should not be any 0 in the mask ptherwise something is wrong with tessellation
#        #print 'Working on', i
#        coord = np.where(data.T == i)
#        dir_x = np.mean([ np.max(coord[0]), np.min(coord[0]) ])
#        dir_y = np.mean([ np.max(coord[1]), np.min(coord[1]) ])
#        ra, dec =  w.all_pix2world(dir_x, dir_y, 0, ra_dec_order=True)
#        directions['Isl_patch_%i' % i] = [float(ra), float(dec)]
#        #print directions[-1]
#
#    # return list
#    return directions
