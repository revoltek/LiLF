#!/usr/bin/env python

import optparse
import numpy as np
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from astropy.io import fits
from astropy.wcs import WCS

import casacore.tables as ct
import losoto.h5parm as h5parm
from losoto.lib_operations import reorderAxes

def plot_screen(img, prefix='', title='', suffix='', wcs=None):
    fig = figure(figsize=(12, 12))
    if wcs is not None:
        ax = fig.add_subplot(111, projection=wcs)
    else:
        ax = fig.add_subplot(111, projection=wcs)

    ax.set_title(title)
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    i = ax.imshow(img, origin='lower', vmin=-0.1, vmax=0.1)
    cbar = fig.colorbar(i, fraction=0.046, pad=0.04)
    cbar.set_label('dTEC')
    fig.savefig(prefix + 'tecscreen_' + suffix + '.svg')
    plt.close(fig)
    del fig


# Interpolate the grid using a nearest neighbour approach.
# https://stackoverflow.com/questions/5551286/filling-gaps-in-a-numpy-array
def interpolate_station(antname, ifstep, h5_stations, size, data, wcs, RA, DEC, gains, names):
    ''' Interpolate the solutions of a given antenna into a smooth screen using radial basis function interpolation.

    Args:
        antname (str): name of the antenna to interpolate.
        ifstep (int): index of the frequency slot to be interpolated.

    Returns:
        antindex (int): index of the antenna being processed.
        screen (ndarray): numpy ndarray containing the interpolated screen.
    '''
    # print('Processing antenna {:s}.'.format(antname))
    interpidx = h5_stations.index(antname)
    # These will be the new coordinates to evaluate the screen over.
    x = np.arange(size)
    y = np.arange(size)
    xx, yy = np.meshgrid(x, y)
    # Do radial basis function interpolation for each time step.
    # First axis of data is time.
    screen = np.empty((data.shape[0], 1, 1, 4, size, size))
    X, Y = np.around(wcs.wcs_world2pix(RA, DEC, 0)).astype(int)
    # data has shape (time, freq, ant, matrix, y, x)
    # gains has shape (time, freq, ant, dir, pol)
    # Iterate over all timeslots.
    for itstep in range(data.shape[0]):
        # Interpolate the gains, not the Re/Im or Amp/Phase separately.
        # gains: ['time', 'freq', 'ant', 'dir', 'pol']
        gXX = gains[itstep, ifstep, interpidx, :, 0]
        gYY = gains[itstep, ifstep, interpidx, :, -1]
        rbfXX = Rbf(X, Y, gXX, smooth=0.0)
        rbfYY = Rbf(X, Y, gYY, smooth=0.0)
        tinterpXX = rbfXX(xx, yy)
        tinterpYY = rbfYY(xx, yy)
        if not np.allclose(tinterpXX[Y, X], gXX, rtol=1e-5):
            raise ValueError('Interpolated screen for polarization XX does not go through nodal points.')
        if not np.allclose(tinterpYY[Y, X], gYY, rtol=1e-5):
            raise ValueError('Interpolated screen for polarization YY does not go through nodal points.')
        del gXX, gYY
        matrix = np.asarray([np.real(tinterpXX), np.imag(tinterpXX), np.real(tinterpYY), np.imag(tinterpYY)])
        screen[itstep, 0, 0, :, :, :] = matrix
    '''
    didx = directions.index('P470')
    fig = figure()
    fig.suptitle(antname)
    ax = fig.add_subplot(111)
    ax.plot(np.real(gains[:, ifstep, interpidx, 0, didx]), 'C0h--', label='ReXX H5parm', markersize=12)
    ax.plot(screen[:, 0, 0, 0, 136, 129], 'C1h--', label='ReXX Rbf Interp. Y,X')
    ax.plot(screen[:, 0, 0, 0, 129, 136], 'C2h--', label='ReXX Rbf Interp. X,Y')
    ax.legend()
    ax.set_xlabel('Time'); ax.set_ylabel('Re G')
    show()
    del fig
    if ('CS' not in antname) and ('RS' not in antname):
        print('== bla ==')
        for i, (rr,dd) in enumerate(zip(RA,DEC)):
            print(rr)
            print(dd)
            orig = tstep[interpidx][i]
            interp = rbf(rr, dd)
            print('Difference TEC and Rbf: {:e}'.format(orig - interp))
            print('== end bla ==')
    fig = figure()
    fig.suptitle(antname)
    ax = fig.add_subplot(111)
    im = ax.imshow(screen[..., -1], origin='lower', extent=[ra_max, ra_min, dec_min, dec_max])
    ax.scatter(RA, DEC, marker='x', color='r')
    ax.set_xlabel('Right ascension'); ax.set_ylabel('Declination')
    fig.colorbar(im)
    fig.savefig(antname + '_{:03d}.png'.format(itstep))
    del fig
    show()
    '''
    return names.index(antname), screen

def main(options):

    size = options.size
    ms = options.ms
    h5p = options.h5parm
    output = options.output
    
    with ct.taql('SELECT NAME FROM {ms:s}::ANTENNA'.format(ms=ms)) as t:
        names = t.getcol('NAME')
    with ct.taql('SELECT REFERENCE_DIR FROM {ms:s}::FIELD'.format(ms=ms)) as t:
        phasecenter = t.getcol('REFERENCE_DIR').squeeze()
    with ct.taql('SELECT TIME FROM {ms:s}'.format(ms=ms)) as t:
        min_time = np.min(t.getcol('TIME'))
        max_time = np.max(t.getcol('TIME'))
    
    #print('Phasecentre', phasecenter)
    # Avoid negative CRVAL1
    if np.rad2deg(phasecenter[0]) < 0:
        phasecenter[0] = phasecenter[0] + (2. * np.pi)
    
    h5 = h5parm.h5parm(h5p)
    ss = h5.getSolset('sol000')
    st = ss.getSoltab('amplitude000')
    # add cut in time to restrict to this ms
    st.setSelection(time={'min':min_time,'max':max_time})
    
    time = st.getAxisValues('time')
    stime = time[0]
    dtime = time[1] - time[0]
    freq = st.getAxisValues('freq')
    sfreq = freq[0]
    dfreq = freq[1] - freq[0]
    
    Nantenna = len(names)
    Ntimes = len(time)
    Nfreqs = len(freq)
    
    # Set the frequency axis (that TEC doens't use) to 150 MHz as central frequency with 50 MHz of bandwidth.
    header = '''SIMPLE  = T / file does conform to FITS standard
BITPIX  = -32           / number of bits per data pixel
NAXIS   = 6             / number of data axes
NAXIS1  = {size:d}      / length of RA axis
NAXIS2  = {size:d}      / length of DEC axis
NAXIS3  = 4
NAXIS4  = {nant:d}      / length of ANTENNA axis
NAXIS5  = 1             / length of FREQ axis
NAXIS6  = {ntimes:d}    / length of TIME axis
EXTEND  = T             / FITS dataset may contain extensions
CTYPE1  = 'RA---SIN'    / Right ascension angle cosine
CRPIX1  = {centerpix:d}
CRVAL1  = {ra:f}
CDELT1  = -{pixelsize:.4f}
CUNIT1  = 'deg     '
CTYPE2  = 'DEC--SIN'     / Declination angle cosine
CRPIX2  = {centerpix:d}
CRVAL2  = {dec:f}
CDELT2  = {pixelsize:.4f}
CUNIT2  = 'deg     '
CTYPE3  = 'MATRIX  '
CRPIX3  = 1.
CDELT3  = 1.
CTYPE4  = 'ANTENNA '
CRPIX4  = 1.
CRVAL4  = 0.
CTYPE5  = 'FREQ    '           / Central frequency
CRPIX5  = 1.
CRVAL5  = {cfreq:f}
CDELT5  = {dfreq:f}
CUNIT5  = 'Hz      '
CTYPE6  = 'TIME    '
CRPIX6  = 1.
CRVAL6  = {stime:f} / Should be an AIPS time
CDELT6  = {dtime:f}'''.format(size=size, centerpix=int(round(size/2)), pixelsize=5/size, nant=Nantenna, ntimes=Ntimes, ra=np.rad2deg(phasecenter[0]), dec=np.rad2deg(phasecenter[1]), cfreq=sfreq, dfreq=dfreq, stime=stime, dtime=dtime)

    # FITS is transposed compared to Numpy.
    # 4 is the "matrix" entry containing Re(XX), Im(XX), Re(YY), Im(YY)
    print (Ntimes, Nfreqs, Nantenna, 4, size, size)
    data = np.empty((Ntimes, Nfreqs, Nantenna, 4, size, size), dtype=np.float32)
    
    # Read in h5parm.
    h5 = h5parm.h5parm(h5p)
    ss = h5.getSolset('sol000')
    sta = ss.getSoltab('amplitude000')
    stp = ss.getSoltab('phase000')
    # add cut in time
    sta.setSelection(time={'min':min_time,'max':max_time})
    stp.setSelection(time={'min':min_time,'max':max_time})
    
    h5_stations = list(st.getAxisValues('ant'))
    # Find nearest pixel for a given direction.
    H = fits.Header.fromstring(header, sep='\n')
    wcs = WCS(H).celestial
    directions = list(st.getAxisValues('dir'))
    sources = ss.getSou()
    RA = []
    DEC = []
    dirs = []
    amps = reorderAxes(sta.getValues()[0], sta.getAxesNames(), ['time', 'freq', 'ant', 'dir', 'pol'])
    phases = reorderAxes(stp.getValues()[0], stp.getAxesNames(), ['time', 'freq', 'ant', 'dir', 'pol'])
    gains = amps * np.exp(1j * phases)
    for d in directions:
        c = sources[d]
        diridx = directions.index(d)
        dirs.append(d)
        RAd, DECd = np.rad2deg(c)
        RA.append(RAd)
        DEC.append(DECd)
        c = wcs.wcs_world2pix(RAd, DECd, 0)
        print('Direction {:s} is at pixel coordinates {:f},{:f}'.format(d, float(c[0]), float(c[1])))
    
        print('Adding direction {:s} at {:f},{:f}'.format(d, RAd, DECd))
        pixel_coords = wcs.wcs_world2pix(RAd, DECd, 0)
        X = int(np.around(pixel_coords[0]))
        Y = int(np.around(pixel_coords[1]))
        print('== Direction corresponds to (x, y) = ({:d}, {:d})'.format(X, Y))
        # Fill it with the solutions for each antenna.
        for istation, station in enumerate(names):
            idx = h5_stations.index(station)
            refidx = h5_stations.index('CS001LBA')
            valsa_tmp = sta.getValues()[0]
            valsa_tmp_r = reorderAxes(valsa_tmp, sta.getAxesNames(), ['time', 'freq', 'ant', 'dir', 'pol'])
            del valsa_tmp
    
            valsp_tmp = stp.getValues()[0]
            valsp_tmp_r = reorderAxes(valsp_tmp, stp.getAxesNames(), ['time', 'freq', 'ant', 'dir', 'pol'])
            del valsp_tmp
            valsp_tmp = valsp_tmp_r - valsp_tmp_r[:, :, np.newaxis, refidx, :, :]
            real = valsa_tmp_r * np.cos(valsp_tmp_r)
            imag = valsa_tmp_r * np.sin(valsp_tmp_r)
    
            reXX = real[..., 0]
            imXX = imag[..., 0]
            reYY = real[..., -1]
            imYY = imag[..., -1]
            matrix = np.asarray([reXX, imXX, reYY, imYY], dtype=np.float32)
            matrix = np.transpose(matrix, (1, 2, 3, 0, 4))
            data[:, :, istation, :, Y, X] = matrix[:, :, idx, :, diridx]
            del matrix
    
    hdu = fits.PrimaryHDU(header=H)
    hdu.data = data
    hdu.writeto(output+'_raw.fits', overwrite=True)
    RA = np.asarray(RA)
    DEC = np.asarray(DEC)
    
    # Inspired by https://stackoverflow.com/questions/38309535/populate-numpy-array-through-concurrent-futures-multiprocessing
    data_int = np.empty(data.shape)
    print('Making Gain screen.')
    for ifreq in range(data.shape[1]):
        print('Processing frequency slot {:d}'.format(ifreq))
        for station in names:
            antenna, screen = interpolate_station(station, ifreq, h5_stations, size, data, wcs, RA, DEC, gains, names)
            data_int[:, ifreq, antenna, :, :, :] = screen[:, 0, 0, ...]
    
    hdu = fits.PrimaryHDU(header=H)
    hdu.data = data_int
    hdu.writeto(output+'_rbf.fits')
    print('Finished interpolating.')

opt = optparse.OptionParser()
opt.add_option('-m','--ms',help='Input MS [no default].',default='')
opt.add_option('-p','--h5parm',help='Input h5parm [no default].',default='')
opt.add_option('-s','--size', type=int, help='Size of the output screen in pixels [default: 256].',default=256)
opt.add_option('-o','--output',help='Output root name [default: output].',default='output')
options, arguments = opt.parse_args()
main(options)
