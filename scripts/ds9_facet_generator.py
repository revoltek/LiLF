#!/usr/bin/env python3

# This script was written by Jakob Maljaars and
# Reinout van Weeren.

from scipy.spatial import Voronoi, voronoi_plot_2d
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import argparse
import sys
import casacore.tables as pt

from shapely.geometry import Polygon
import shapely.geometry
import shapely.ops
import tables


def read_dir_fromh5(h5):
    """
    Read in the direction info from an H5 file
    Parameters
    ----------
    h5 : str
        h5 filename

    Returns
    ----------
    sourcedir: np.array
        Array containing directions (ra, dec in units of radians)
    """

    H5 = tables.open_file(h5, mode="r")
    sourcedir = H5.root.sol000.source[:]["dir"]
    sourcename = H5.root.sol000.source[:]["name"].astype('str')
    # if len(sourcedir) < 2:
        # print("Error: H5 seems to contain only one direction")
        # sys.exit(1)
    H5.close()
    return sourcedir, sourcename


def makeWCS(centreX, centreY, refRA, refDec, crdelt=0.066667):
    """
    Makes simple WCS object.
    Parameters
    ----------
    centreX : int
        Centre x pixel
    centreY : int
        Centre y pixel
    refRA : float
        Reference RA in degrees
    refDec : float
        Reference Dec in degrees
    crdelt: float, optional
        Delta in degrees for sky grid. Default value is 0.066667 (=4amin)
    Returns
    -------
    w : astropy.wcs.WCS object
        A simple TAN-projection WCS object for specified reference position
    """

    w = WCS(naxis=2)
    w.wcs.crpix = [centreX, centreY]
    w.wcs.cdelt = np.array([-crdelt, crdelt])
    w.wcs.crval = [refRA, refDec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.set_pv([(2, 1, 45.0)])
    return w


def convert_to_deg(array_ra, array_dec):

    try:
        # Degree format
        new_ra = Angle(array_ra, unit="degree")
        new_dec = Angle(array_dec, unit="degree")
    except ValueError:
        # Skymodel.txt format
        new_ra = Angle(array_ra, unit="hourangle")
        new_dec = Angle(np.char.replace(array_dec, ".", ":", 2), unit="degree")

    return [new_ra.deg, new_dec.deg]


def generate_centroids_from_source_catalog(catalog_file, npoints, w):
    """
    Generate centroids from a source cataloge, such as gleam-osm.

    Parameters
    ----------
    catalog_file : str
        Source catalogue file
    npoints : int
        Number of (brightest) sources to select
    w : astropy.wcs.WCS object
        [description]

    Returns
    -------
    np.2darray
        Numpy (npoints, 2) array with pixel coordinates of . Dimension:
    """

    catalog = np.genfromtxt(catalog_file, delimiter=",", dtype=str, encoding=None)
    source_idx = np.argsort(catalog[:, 2])[: -npoints - 1 : -1]
    catalog = np.char.strip(catalog)

    # Search the keywords "ra" and "dec" in the first line of the catalog
    # to get the right indexes.
    # If not found use the indexes 0 and 1 as in the "gleam.osm" catalog
    try:
        index_ra = list(np.char.lower(catalog[0, :])).index("ra")
        index_dec = list(np.char.lower(catalog[0, :])).index("dec")
    except:
        index_ra = 0
        index_dec = 1

    # Convert Ra/Dec unit to degrees
    [ra_coords, dec_coords] = convert_to_deg(
        catalog[source_idx, index_ra], catalog[source_idx, index_dec]
    )

    x, y = w.wcs_world2pix(ra_coords, dec_coords, 1)
    return np.vstack((x.flatten(), y.flatten())).T


def tessellate(x_pix, y_pix, w, dist_pix, bbox, nouter=64, plot_tessellation=True):
    """
    Returns Voronoi tessellation vertices

    Parameters
    ----------
    x_pix : array
        Array of x pixel values for tessellation centers
    y_pix : array
        Array of y pixel values for tessellation centers
    w : WCS object
        WCS for transformation from pix to world coordinates
    dist_pix : float
        Distance in pixels from center to outer boundary of facets
    nouter : int
        Number of points to generate on the outer boundary for constraining
        the Voronoi tessellation. Defaults to 64
    plot_tessellation : bool
        Plot tessellation

    Returns
    -------
    list, np.2darray
        List of shapely Polygons, and np.2darray of corresponding (Voronoi) points (ra,dec in degrees)
    """

    # Get x, y coords for directions in pixels. We use the input calibration sky
    # model for this, as the patch positions written to the h5parm file by DP3 may
    # be different
    xy = []
    for RAvert, Decvert in zip(x_pix, y_pix):
        xy.append((RAvert, Decvert))

    # Generate array of outer points used to constrain the facets
    means = np.ones((nouter, 2)) * np.array(xy).mean(axis=0)
    offsets = []
    angles = [np.pi / (nouter / 2.0) * i for i in range(0, nouter)]
    for ang in angles:
        offsets.append([np.cos(ang), np.sin(ang)])
    scale_offsets = dist_pix * np.array(offsets)
    outer_box = means + scale_offsets

    # Tessellate and clip
    points_all = np.vstack([xy, outer_box])
    vor = Voronoi(points_all)

    # Filter out the infinite regions
    region_indices = [
        region_idx
        for region_idx in vor.point_region
        if -1 not in vor.regions[region_idx]
    ]
    polygons = []
    for idx in region_indices:
        vertex_coords = vor.vertices[vor.regions[idx]]
        polygons.append(Polygon(vertex_coords))

    clipped_polygons = []
    for polygon in polygons:
        # facet_poly = Polygon(facet)
        clipped_polygons.append(polygon_intersect(bbox, polygon))

    if plot_tessellation:
        import matplotlib.pyplot as plt

        [plt.plot(*poly.exterior.xy) for poly in clipped_polygons]
        plt.xlabel("Right Ascension [pixels]")
        plt.ylabel("Declination [pixels]")
        plt.axis("square")
        plt.tight_layout()
        plt.show()

    verts = []
    for poly in clipped_polygons:
        verts_xy = poly.exterior.xy
        verts_deg = []
        for x, y in zip(verts_xy[0], verts_xy[1]):
            x_y = np.array([[y, x, 0.0, 0.0]])
            ra_deg, dec_deg = w.wcs_pix2world(x, y, 1)
            verts_deg.append((ra_deg, dec_deg))
        verts.append(verts_deg)

    # Reorder to match the initial ordering
    ind = []
    for poly in polygons:
        for j, (xs, ys) in enumerate(zip(x_pix, y_pix)):
            if poly.contains(shapely.geometry.Point(xs, ys)):
                ind.append(j)
                break
    verts = [verts[i] for i in ind]

    ra_point, dec_point = w.wcs_pix2world(x_pix, y_pix, 1)
    return [Polygon(vert) for vert in verts], np.vstack((ra_point, dec_point)).T


def generate_centroids(
    xmin, ymin, xmax, ymax, npoints_x, npoints_y, distort_x=0.0, distort_y=0.0
):
    """
    Generate centroids for the Voronoi tessellation. These points are essentially
    generated from a distorted regular grid.

    Parameters
    ----------
    xmin : float
        Min-x pixel index, typically 0
    ymin : float
        Min-y pixel index, typically 0
    xmax : float
        Max-x pixel index, typically image width
    ymax : float
        Max-y pixel index, typically image height
    npoints_x : int
        Number of points to generate in width direction
    npoints_y : int
        Number of points to generate in height direction
    distort_x : float, optional
        "Cell width" fraction by which to distort the x points, by default 0.0
    distort_y : float, optional
        "Cell height" fraction by which to distory the y points, by default 0.0

    Returns
    -------
    X,Y : np.1darray
        Flattened arrays with X,Y coordinates
    """

    x_int = np.linspace(xmin, xmax, npoints_x)
    y_int = np.linspace(ymin, ymax, npoints_y)

    np.random.seed(0)

    # Strip the points on the boundary
    x = x_int[1:-1]
    y = y_int[1:-1]
    X, Y = np.meshgrid(x, y)

    xtol = np.diff(x)[0]
    dX = np.random.uniform(low=-distort_x * xtol, high=distort_x * xtol, size=X.shape)
    X = X + dX

    ytol = np.diff(y)[0]
    dY = np.random.uniform(low=-distort_x * ytol, high=distort_y * ytol, size=Y.shape)
    Y = Y + dY
    return X.flatten(), Y.flatten()


def polygon_intersect(poly1, poly2):
    """
    Returns the intersection of polygon2 with polygon1
    """
    clip = poly1.intersection(poly2)
    return clip


def write_ds9(fname, polygons, points=None, names=None):
    """
    Write ds9 regions file, given a list of polygons
    and (optionally) a set of points attached to

    Parameters
    ----------
    fname : str
        Filename for output file
    polygons : list
        List of shapely.Polygons
    points : np.2darray, optional
        Array of point coordinates (ra, dec in degrees) that should be
        attached to a facet, by default None
    """

    if points is not None:
        assert (
            len(polygons) == points.shape[0]
        ), "Number of polygons and number of points should match"

    # Write header
    header = [
        "# Region file format: DS9 version 4.1",
        'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1',
        "fk5",
        "\n",
    ]
    with open(fname, "w") as f:
        f.writelines("\n".join(header))
        polygon_strings = []
        for i, polygon in enumerate(polygons):
            poly_string = "polygon("
            xv, yv = polygon.exterior.xy
            for (x, y) in zip(xv[:-1], yv[:-1]):
                poly_string = f"{poly_string}{x:.5f},{y:.5f},"
            # Strip trailing comma
            poly_string = poly_string[:-1] + ")"
            if points is not None:
                poly_string += f"\npoint({points[i, 0]:.5f}, {points[i, 1]:.5f})"
                if names is not None:
                    poly_string += f" # text={names[i]}"
            polygon_strings.append(poly_string)
        f.write("\n".join(polygon_strings))


def main(args):
    # get phase centre from the ms in units of degrees
    t = pt.table(args.ms + "::FIELD", ack=False)
    phasedir = t.getcol("PHASE_DIR").squeeze()
    cphasedir = SkyCoord(
        ra=phasedir[0] * u.rad, dec=phasedir[1] * u.rad
    )  # astropy coordinate
    phaseCentreRa = cphasedir.ra.degree
    phaseCentreDec = cphasedir.dec.degree

    # Pixel "resolution" (in degrees!)
    dl_dm = args.pixelscale / 60.0 / 60.0  # in units of degree

    # Image size (in pixels)
    xmin = 0
    xmax = args.imsize
    ymin = 0
    ymax = args.imsize
    centreX = (xmax - xmin) // 2 + 1
    centreY = (ymax - ymin) // 2 + 1

    # To cut the Voronoi tessellation on the bounding box, we need
    # a "circumscribing circle"
    dist_pix = np.sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2)

    # Make World Coord Stystem transform object
    w = makeWCS(centreX, centreY, phaseCentreRa, phaseCentreDec, dl_dm)

    single_dir = False # If there is only one direction in the h5 file, we should return rectangle over the full image

    if args.h5:
        # load in the directions from the H5
        sourcedir, sourcename = read_dir_fromh5(args.h5)
        if len(sourcedir) < 2:
            single_dir = True

        # make ra and dec arrays and coordinates c
        ralist = sourcedir[:, 0]
        declist = sourcedir[:, 1]
        c = SkyCoord(ra=ralist * u.rad, dec=declist * u.rad)

        # convert from ra,dec to x,y pixel
        x, y = w.wcs_world2pix(c.ra.degree, c.dec.degree, 1)
    elif args.sourcecatalog:
        # Add two points to account for outer points that will be stripped away
        x_background, y_background = generate_centroids(
            xmin, ymin, xmax, ymax, args.backgroundfacets + 2, args.backgroundfacets + 2
        )
        xy = generate_centroids_from_source_catalog(
            args.sourcecatalog[0], int(args.sourcecatalog[1]), w
        )
        x = np.hstack((x_background, xy[:, 0]))
        y = np.hstack((y_background, xy[:, 1]))
    else:
        raise Exception(
            "Use either --h5 or --sourcecatalog, to generate facets from an h5 file or a source catalogue, respectively."
        )

    bbox = Polygon([(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)])
    if single_dir:
        pointsarr = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
        skyarr = []
        for x,y in pointsarr:
            xsky,ysky = w.wcs_pix2world(x, y, 1)
            skyarr.append((xsky,ysky))
        skybox = Polygon(skyarr)

        write_ds9(
            args.outputfile, np.array([skybox]), points=np.array([[phaseCentreRa,phaseCentreDec]]) if args.writevoronoipoints else None, names=sourcename if args.h5 else None
        )

    else:
        facets, points = tessellate(
            x, y, w, dist_pix, bbox, plot_tessellation=args.plottessellation
        )

        write_ds9(
            args.outputfile, facets, points=points if args.writevoronoipoints else None, names=sourcename if args.h5 else None
        )



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Make DS9 Voronoi tessellation region file for WSClean"
    )
    parser.add_argument("--ms", help="boxfile", type=str, required=True)
    parser.add_argument(
        "--h5",
        help="H5Parm file to be used for generating the facets. Excludes the use of --source-catalog",
        type=str,
    )
    parser.add_argument(
        "--sourcecatalog",
        help="Path to source catalogue, followed by the number of brightest points that should be used.",
        nargs=2,
    )
    parser.add_argument(
        "--backgroundfacets",
        help="Set the resolution for the number of background facets, in case facets are generated from a source catalog.",
        type=int,
        default=4,
    )
    parser.add_argument(
        "--imsize",
        help="image size, required if boxfile is not used",
        type=int,
        default=8192,
    )
    parser.add_argument(
        "--pixelscale",
        help="pixels size in arcsec, default=1.5",
        type=float,
        default=1.5,
    )
    parser.add_argument(
        "--plottessellation", help="Plot tessellation", action="store_true"
    )
    parser.add_argument(
        "--writevoronoipoints",
        help="Write the Voronoi input points to the output regions file",
        action="store_true",
    )
    parser.add_argument(
        "--outputfile",
        help="Name of output file, defaults to facets.reg",
        type=str,
        default="facets.reg",
    )
    args = parser.parse_args()
    main(args)

