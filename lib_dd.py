import os, sys, itertools
import numpy as np
from astropy.table import Table
from astropy.coordinates import Angle, SkyCoord, match_coordinates_sky
from astropy.io import fits as pyfits
from astropy import wcs as pywcs
import astropy.units as u
import pyregion
from pyregion.parser_helper import Shape
import bdsf
try:
    from scipy.spatial import Voronoi
except:
    logger.error("Load latest scipy with 'use Pythonlibs'")
    sys.exit(1)

from lib_log import logger

def table_to_circ_region(table, outfile, racol='RA', deccol='DEC', sizecol='size', color='red', label=True):
    """
    Get a table with ra, dec, size and generate a circular ds9 region 
    TODO: if cat is given, add a small circle around all sources on the edge
    """

    regions = []
    for i, r in enumerate(table):
        s = Shape('circle', None)
        s.coord_format = 'fk5'
        s.coord_list = [ r[racol], r[deccol], r[sizecol] ] # ra, dec, radius
        s.coord_format = 'fk5'
        s.attr = ([], {'width': '2', 'point': 'cross',
                       'font': '"helvetica 16 normal roman"'})
        if label: s.comment = 'color={} text="{}"'.format(color, str(i))
        else: s.comment = 'color={}'.format(color)
        regions.append(s)

    regions = pyregion.ShapeList(regions)
    check_rm(outfile)
    regions.write(outfile)


def make_directions_from_skymodel(filename, outdir='regions/', flux_min_Jy=1.0, size_max_arcmin=3.0,
    directions_separation_max_arcmin=5.0, directions_max_num=20, flux_min_for_merging_Jy=0.2):
    """
    Selects appropriate calibrators from srl file
    Parameters
    ----------
    filename : srl file
        Skymodel made by grouping clean components of dir-independent model
    outdir: string
        Directory where to save the ds9 regions
    flux_min_Jy : float
        Minimum flux density for a calibrator in Jy
    size_max_arcmin : float
        Maximum size for a calibrator in arcmin
    directions_separation_max_arcmin : float
        Maximum separation in arcmin between two calibrators for gouping into a
        single direction
    directions_max_num : int, optional
        Limit total number of directions to this value
    flux_min_for_merging_Jy : float, optional
        Minimum peak flux for a source to be considered for merging
    Returns
    -------
    table : astropy table
        Table with direction information
    """
    # open astropy table
    t = Table.read(filename, format='fits')['RA','DEC','Maj','Peak_flux','Total_flux'] # restrict to some cols
    t.rename_column('Maj', 'dd_size') # use Maj as proxy for region size
    logger.info('# sources initial: %i' % len(t))

    # exclude sources that are too extended
    t_large = t[ (t['dd_size'] >= size_max_arcmin*u.arcmin) ]
    t = t[ (t['dd_size'] < size_max_arcmin*u.arcmin) ]
    logger.info('# sources after cut on size: %i' % len(t))
    logger.info('# large sources: %i' % len(t_large))
    if len(t) == 0:
        logger.critical("No sources found that meet the specified max size criterion.")
        sys.exit(1)

    t['dd_size'] *= 3. # now that we cut on size, enlarge all the regions to peak up artifacts and sidelobes around dd calibrators
    # min size, set to 3 arcmin
    t[ (t['dd_size'] < 3.*u.arcmin) ]['dd_size'] = 3.*u.arcmin

    # exclude sources that are too faint
    t_large = t_large[ (t_large['Peak_flux'] > flux_min_for_merging_Jy) ]
    t = t[ (t['Peak_flux'] > flux_min_for_merging_Jy) ]
    logger.info('# sources after cut min flux for merging: %i' % len(t))
    if len(t) == 0:
        logger.critical("No sources found above %f Jy." % flux_min_for_merging_Jy )
        sys.exit(1)

    t.sort('Peak_flux')
    t.reverse()

    # combine nearby sources
    for s in t:
        # if ra/dec changes, continue finding nearby sources until no-sources are found
        updated = True
        while updated:
            dists = SkyCoord(ra=s['RA']*u.degree, dec=s['DEC']*u.degree).separation(SkyCoord(ra=t['RA'], dec=t['DEC']))
            updated = False
            for i, dist in enumerate(dists):
                if dist < directions_separation_max_arcmin*u.arcmin and dist > 0.*u.degree:
                    # if a source is dominant keep that at the center of the patch
                    if t['Peak_flux'][i] > 3*s['Peak_flux']:
                        s['RA'] = t['RA'][i]
                        s['DEC'] = t['DEC'][i]
                        updated = True
                    # other wise weighted mean
                    elif t['Peak_flux'][i] > s['Peak_flux']:
                        s['RA'] = (s['RA']*s['Peak_flux'] + t['RA'][i]*t['Peak_flux'][i])/(s['Peak_flux']+t['Peak_flux'][i])
                        s['DEC'] = (s['DEC']*s['Peak_flux'] + t['DEC'][i]*t['Peak_flux'][i])/(s['Peak_flux']+t['Peak_flux'][i])
                        updated = True

                    s['dd_size'] = max(s['dd_size'], t['dd_size'][i]) + dist.degree

                    s['Total_flux'] += t['Total_flux'][i]
                    s['Peak_flux'] = max(s['Peak_flux'], t['Peak_flux'][i])

                    t.remove_rows(i)
    logger.info('# sources after combining close-by sources: %i' % len(t))

    # Filter patches on total flux density limit
    t = t[ (t['Total_flux'] > flux_min_Jy) ]
    logger.info('# sources after cut min flux: %i' % len(t))
    if len(t) == 0:
        logger.critical("No sources or merged groups found that meet the specified "
            "min total flux density criterion.")
        sys.exit(1)

    # Trim directions list to get directions_max_num of directions
    t.sort('Peak_flux')
    t.reverse()
    if directions_max_num is not None:
        t = t[:directions_max_num]
        logger.info('# sources after cut on max directions: %i' % len(t))

    for s in t_large:
        dists = SkyCoord(ra=s['RA']*u.degree, dec=s['DEC']*u.degree).separation(SkyCoord(ra=t['RA'], dec=t['DEC']))
        for i, dist in enumerate(dists):
            if dist < directions_separation_max_arcmin*u.arcmin:
                t['Total_flux'][i] += s['Total_flux']
                t['dd_size'][i] = max(s['dd_size'], t['dd_size'][i]) + dist.degree

    # sort on a weighted mix of total and peak flux
    t['Comb_flux'] = 0.33*t['Total_flux']+0.66*t['Peak_flux']
    t.sort('Comb_flux')
    t.reverse()

    # Writedd global region files
    table_to_circ_region(t, outdir+'/all.reg')

    # save source by source (size is converted in a radius, this is conservative)
    for i, s in enumerate(t):
        table_to_circ_region(t[i:i+1], outdir+'/ddcal%02i.reg' % i )

    t['name'] = ['ddcal%02i' % i for i in xrange(len(t))]

    t.remove_column('Comb_flux')
    return t


def make_directions_from_img(imagename, outdir='regions/', target_flux_jy=10, bright_source_jy=5., size_max_arcmin=3., trials=None):
    """
    fitsfile = selfcal model, used for coordinates
    outdir = Directory where to save the ds9 regions
    target_flux_jy = target flux per facet, Jy
    bright_source_jy = these sources are centered on their facet
    size_max_arcmin = Maximum size for a source to be considered, arcmin
    trials = number of sources to use as possible region center, if None, use all
    """

    # Run pybdsf
    logger.info('Finding directions...')
    if not os.path.exists('regions/DIEcatalog.fits'):
        bdsf_img = bdsf.process_image(imagename, rms_box=(55,12), \
            thresh_pix=5, thresh_isl=3, atrous_do=False, atrous_jmax=3, \
            adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(80,20), \
            quiet=True)
        check_rm('regions')
        os.makedirs('regions')
        bdsf_img.write_catalog(outfile='regions/DIEcatalog.fits', catalog_type='srl', format='fits')

    t = Table.read('regions/DIEcatalog.fits', format='fits')['RA','DEC','Maj','Peak_flux','Total_flux'] # restrict to some cols
    t.rename_column('Maj', 'size') # use Maj as proxy for region size
    # exclude sources that are too extended
    t = t[ (t['size'] < size_max_arcmin*u.arcmin) ]
    t['size'] *= 3 # enlarge all sizes
    logger.info('# sources after cut on size: %i' % len(t))
    total_flux = np.sum(t['Total_flux'])
    logger.info('# sources initial: %i -- total flux = %.2f Jy' % (len(t),total_flux) )
    if trials is None: trials = len(t)

    t.sort('Peak_flux')
    t.reverse()

    # make ddcal table 
    ddcal                    = Table()
    ddcal['RA']              = np.zeros(trials)
    ddcal['DEC']             = np.zeros(trials)
    ddcal['dd_size']         = np.zeros(trials)
    ddcal['facet_size']      = np.zeros(trials)
    ddcal['Peak_flux']       = np.zeros(trials)
    ddcal['Total_flux']      = np.zeros(trials)
    ddcal['RA'].unit         = u.degree
    ddcal['DEC'].unit        = u.degree
    ddcal['dd_size'].unit    = u.degree
    ddcal['facet_size'].unit = u.degree
    ddcal['Total_flux'].unit = u.Jy
    ddcal['Peak_flux'].unit  = 'Jy/beam'

    # find ddcal as regions of enough flux around bright sources
    idx_sources = []
    for idx_ddcal, dd in enumerate(ddcal):
        #print "### source number %i" % idx_ddcal
        #idx_ddcal = [i for i in xrange(len(t)) if not (i in idx_good_cals or i in idx_bad_cals)][0] # first usable source
        s = t[idx_ddcal] # use this source as a center for new ddcal region
        dd['RA'] = s['RA']
        dd['DEC'] = s['DEC']
        dd['dd_size'] = s['size']
        dd['Total_flux'] = 0.
        dd['Peak_flux'] = 0.
        dists = SkyCoord(ra=dd['RA']*u.degree, dec=dd['DEC']*u.degree).separation(SkyCoord(ra=t['RA'], dec=t['DEC'])).degree
        idx_closests = [idx for (dist,idx) in sorted(zip(dists,xrange(len(dists))))] # idx list from closest to farthest
        idx_sources.append([])
        for idx_closest in idx_closests:

            # first cycle matches at distance=0 the calibrator
            s = t[idx_closest]
            dd['Total_flux'] += s['Total_flux']
            dd['Peak_flux'] = max(dd['Peak_flux'], s['Peak_flux'])

            idx_sources[idx_ddcal].append(idx_closest) # keep track of all sources in this region
            if dd['Total_flux'] > target_flux_jy: break

        if len(idx_sources[idx_ddcal]) > 1:
            dd['dd_size'] = max(SkyCoord(ra=dd['RA']*u.degree, dec=dd['DEC']*u.degree).separation(SkyCoord(ra=t[idx_sources[idx_ddcal]]['RA'], dec=t[idx_sources[idx_ddcal]]['DEC'])).degree)

    # some cleaning up:
    # remove large regions
    toremove = np.where(ddcal['dd_size'] > .8)[0]
    ddcal.remove_rows(toremove)
    for r in sorted(toremove, reverse=True):
        del idx_sources[r]
    logger.info('Number of ddcal after size cut: %i' % len(ddcal))

    # for bright sources keep only centered regions
    idx_brights = np.where(t['Total_flux'] > bright_source_jy)[0]
    toremove = []
    for i in xrange(len(idx_sources)):
        for idx_bright in idx_brights:
            if idx_bright in idx_sources[i][1:]:
                toremove.append(i)
    ddcal.remove_rows(toremove)
    for r in sorted(toremove, reverse=True):
        del idx_sources[r]
    logger.info('Number of ddcal after bright cal cut: %i' % len(ddcal))

    # sort in size
    idx_sources = [idx for (size,idx) in sorted(zip(ddcal['dd_size'],idx_sources))]
    ddcal.sort('dd_size')
        
    # TODO: wrong, regions must be independent or calibration is messed up
    # finally retain only independent regions
    toremove = []
    for i, dd in enumerate(ddcal):
        for j, dd2 in enumerate(ddcal[:i]):
            if j in toremove: continue
            n_match = len([s for s in idx_sources[i] if s in idx_sources[j]])
            if n_match > 0.25*len(idx_sources[i]):
                toremove.append(i)
                break
    ddcal.remove_rows(toremove)
    for r in sorted(toremove, reverse=True):
        del idx_sources[r]
    logger.info('Number of ddcal after cut on overlaps: %i' % len(ddcal))

    ddcal['name'] = ['ddcal%02i' % i for i in xrange(len(ddcal))]

    # save ddcal regions
    for i, s in enumerate(ddcal):
        table_to_circ_region(t[idx_sources[i]], outdir+'/ddcal%02i.reg' % i, color='green')
    # save all regions in one file for inspection
    table_to_circ_region(ddcal, outdir+'/all.reg', sizecol='dd_size')
    return ddcal


def make_voronoi_reg(directions, fitsfile, outdir='regions/', beam_reg='', png=None):
    """
    Take a list of coordinates and an image and voronoi tesselate the sky.
    It saves ds9 regions of the facets and and return facet sizes

    directions : dict with {'dir0':[ra,dec], 'dir1':[ra,dec]...}
    firsfile : model fits file to tassellate (used for coordinates)
    outdir : dir where to save regions
    beam_reg : a ds9 region showing the the primary beam, exclude directions outside it
    """

    import lib_img
    logger.debug("Image used for tasselation reference: "+fitsfile)
    fits = pyfits.open(fitsfile)
    hdr, data = lib_img.flatten(fits)
    w = pywcs.WCS(hdr)
    pixsize = np.abs(hdr['CDELT1'])

    # Add facet size column
    ras = np.array([directions[d][0].degree for d in directions])
    decs = np.array([directions[d][1].degree for d in directions])
    x, y = w.all_world2pix(ras, decs, 0, ra_dec_order=True)

    x_c = data.shape[0]/2.
    y_c = data.shape[1]/2.

    if beam_reg == '':
        # no beam, use all directions for facets
        idx_for_facet = range(len(directions))
    else:
        r = pyregion.open(beam_reg)
        beam_mask = r.get_mask(header=hdr, shape=data.shape)
        beamradius_pix = r[0].coord_list[2]/pixsize
        idx_for_facet = []
        for i, dd in enumerate(t):
            if beam_mask[t['x'][i],t['y'][i]] == True:
                idx_for_facet.append(i)

    # convert to pixel space (voronoi must be in eucledian space)
    x1 = 0
    y1 = 0
    x2 = data.shape[0]
    y2 = data.shape[1]

    # do tasselization
    vor = Voronoi(np.array((x[idx_for_facet], y[idx_for_facet])).transpose())
    box = np.array([[x1,y1],[x2,y2]])
    impoly = voronoi_finite_polygons_2d_box(vor, box)

    # save regions
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
        regionfile = outdir+directions.keys()[idx_for_facet[i]]+'.reg'
        regions.write(regionfile)

        #if beam_reg != '': npix = size_from_reg(fitsfile, [regionfile, beam_reg], [ras[idx_for_facet[i]], decs[idx_for_facet[i]]])
        #else: npix = size_from_reg(fitsfile, [regionfile], [ras[idx_for_facet[i]], decs[idx_for_facet[i]]])
    
    # add names for all.reg
    for d_name, d_coord in directions.iteritems():
        s = Shape('circle', None)
        s.coord_format = 'fk5'
        s.coord_list = [ d_coord[0].degree, d_coord[1].degree, 0.01 ] # ra, dec, radius
        s.coord_format = 'fk5'
        s.attr = ([], {'width': '1', 'point': 'cross',
                       'font': '"helvetica 16 normal roman"'})
        s.comment = 'color=white text="%s"' % d_name
        all_s.append(s)

    regions = pyregion.ShapeList(all_s)
    regionfile = outdir+'all.reg'
    regions.write(regionfile)
    logger.debug('There are %i calibrator within the PB and %i outside (no facet).' % (len(idx_for_facet), len(directions) - len(idx_for_facet)))

    # plot tesselization
    if png is not None:
        import matplotlib.pyplot as pl
        pl.figure(figsize=(8,8))
        ax1 = pl.gca()
        ax1.plot(x,y,'*',color='red')
        for i, d in enumerate(directions): ax1.text(x[i], y[i], d, fontsize=15)
        if beam_reg != '':
            c1 = pl.Circle((x_c, y_c), beamradius_pix, color='g', fill=False)
            ax1.add_artist(c1)
        ax1.plot([x1,x1,x2,x2,x1],[y1,y2,y2,y1,y1])
        for p in impoly:
            pp = p.transpose()
            ax1.plot(pp[0],pp[1])
        ax1.set_xlabel('RA (pixel)')
        ax1.set_ylabel('Dec (pixel)')
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
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

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

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    regions, imvertices = new_regions, np.asarray(new_vertices)
    #return new_regions, np.asarray(new_vertices)

    ## now force them to be in the bounding box
    poly = np.asarray([imvertices[v] for v in regions])

    newpoly = []

    for p in poly:
        polyPath = mplPath.Path(p)
        newpolyPath = polyPath.clip_to_bbox(bbox)
        pp = newpolyPath.vertices.transpose()
        newpoly.append(pp.transpose())

    return np.asarray(newpoly)
