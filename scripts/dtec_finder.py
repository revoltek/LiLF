#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "J. Boxelaar"

import logging
import numpy as np
import matplotlib.pyplot as plt
import losoto.h5parm as losoto
from scipy.signal import lombscargle, find_peaks, medfilt
from scipy.interpolate import make_smoothing_spline
from scipy.optimize import curve_fit

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')
logging.info('dTEC fitter - Jort Boxelaar')

def dTEC(freqs, dtec):
    model_tec = -8067*(freqs/60)**-1 * dtec
    return model_tec - model_tec[len(model_tec)//2]

def fit_lombscargle(phases, data, tint, ant="CS002LBA", freq_split=0, plot = True, dtec_sign=1):
    phases = phases[:,:,data['ant']==ant,0].squeeze().T
    stacked_phase = np.concatenate([phases[freq_split:,time] for time in range(tint[0],tint[1])])
    stacked_freqs = np.concatenate([data['freq'][freq_split:]/1.e6 for __ in range(tint[0],tint[1])])

    index_order = np.argsort(stacked_freqs)
    freqs = stacked_freqs[index_order]
    angles = stacked_phase[index_order]

    freqs = freqs[np.isnan(angles) == False]
    angles = angles[np.isnan(angles) == False]

    if len(freqs) == 0:
        logging.warning(f"time interval {tint[0]}-{tint[1]} fully flagged for ant: {ant}")
        return np.nan, np.nan, np.array([]), np.array([])
    
    dtec_freqs = (-8067*(freqs/60)**-1) * np.pi / 180
    bandwidth = dtec_freqs.ptp()
    n = len(dtec_freqs)
    lombfreqs = np.linspace(1/bandwidth, n/bandwidth, 5*n)
    
    periodogram = lombscargle(dtec_freqs, np.sin(angles), freqs=lombfreqs)
    periodogram = periodogram/periodogram.max() 
    
    # find peaks in periodogram with largest dtec (select last five peaks with relative height > 0.1)
    # after that select the two highest peaks
    peaks, __ = find_peaks(periodogram, height=0.1)
    dtecs = peaks[np.argsort(lombfreqs[peaks])][::-1][:5]
    dtecs = dtecs[np.argsort(periodogram[dtecs])][::-1][:2]
    if plot:
        plt.cla(); plt.clf()
        print(lombfreqs[dtecs])
        for dtec in dtecs:
            plt.axvline(lombfreqs[dtec], color='red')
            plt.plot(lombfreqs, periodogram)
        plt.show()
    
    if len(dtecs) == 0:
        # if highest dtec peaks are not found, default to old method
        print("no peaks found at  ", ant, tint)
        peaks, __ = find_peaks(periodogram)
        dtecs = peaks[np.argsort(periodogram[peaks])][::-1][:2]

    if len(dtecs) == 1:
        best_dtec = lombfreqs[dtecs[0]]
    elif periodogram[dtecs[0]]/periodogram[dtecs[1]] > 2:
        best_dtec = lombfreqs[dtecs[0]]
    else:
        best_dtec = np.sum(periodogram[dtecs] * lombfreqs[dtecs])/np.sum(periodogram[dtecs])
    
    unwrapped = np.unwrap(angles) - np.unwrap(angles)[len(angles)//2]
    if best_dtec > 2 and dtec_sign != 0:
        # sign is determined by the sign of previous timesteps
        sign = dtec_sign
    else:
        minus = np.sum((unwrapped - dTEC(freqs, -best_dtec)/180*np.pi)**2)
        plus = np.sum((unwrapped - dTEC(freqs, best_dtec)/180*np.pi)**2)
        if minus < plus:
            sign = -1
        else:
            sign = 1
        
    if plot:
        fig, axs = plt.subplots(1,2, figsize=(10,5))
        axs[0].scatter(freqs, np.unwrap(angles) - np.unwrap(angles)[len(angles)//2], s=1)
        axs[0].plot(freqs, dTEC(freqs, sign * best_dtec)/180*np.pi, color='black')
        #axs[1].plot(stacked_freqs, stacked_phase, color='black',alpha=0.5, linewidth=1)
        axs[1].scatter(stacked_freqs, stacked_phase, s=1)
        plt.show()#; plt.cla(); plt.clf()

    return sign * best_dtec, 0.001, freqs, angles

def open_sols(path, soltab="phase000", constrain=False, apply_flags=False):
    sols = losoto.h5parm(path, readonly = True)
    solset = sols.getSolset(solset = "sol000")
    
    soltab = solset.getSoltab(soltab = soltab)  
    phases, data = soltab.getValues()
    if constrain:
        phase_shift = phases.copy()
        phase_shift[phases > np.pi] -= 2*np.pi
        phase_shift[phases < -np.pi] += 2*np.pi
        phases = phase_shift
    
    if apply_flags:
        flags, __ = soltab.getValues(weight=True)
        phases[flags < 0.5] = np.nan
        logging.info("flagged data percentage:", np.sum(flags < 0.5)/len(flags))
    sols.close()
    return phases, data

def quality_check(full_dtec, full_dtec_err, full_timesteps, full_ants):
    fine = 0
    for i, ant in enumerate(full_ants):
        if "DE" not in ant:
            continue
        
        bins = np.linspace(full_timesteps[i].min(), full_timesteps[i].max(), 3)
        digitized = np.digitize(full_timesteps[i], bins)
        noise_l = len(full_dtec[i][(digitized == 1) & (50*full_dtec_err[i] > 0.1)]) + 1
        noise_h = len(full_dtec[i][(digitized == 2) & (50*full_dtec_err[i] > 0.1)]) + 1
        
        if noise_h/noise_l > 2:
            fine += 1
        elif noise_h/noise_l < 0.5:
            fine -= 1

    if fine >= 2:
        logging.info("second half too noisy")
        return 'first'
    elif fine <= -2:
        logging.info("first half too noisy")
        return 'second'
    else:
        logging.info("noise levels are fine")
        return 'all'
    


def get_ref_angle(angs):
    all_nan = True
    width = 5
    while all_nan:
        angle_range = angs[len(angs)//2-width:len(angs)//2+width]
        if np.all(np.isnan(angle_range)):
            width += 5
        else:
            all_nan = False
    return np.nanmean(angle_range)

def fit_unwrap(phases, data, tint, ant="CS002LBA", freq_split=0, plot = True):
    phases = phases[:,:,data['ant']==ant,0].squeeze().T
    stacked_phase = np.concatenate([phases[freq_split:,time] for time in range(tint[0],tint[1])])
    stacked_freqs = np.concatenate([data['freq'][freq_split:]/1.e6 for __ in range(tint[0],tint[1])])

    index_order = np.argsort(stacked_freqs)
    freqs = stacked_freqs[index_order]
    angles = stacked_phase[index_order]
    angles_notnan = angles[~np.isnan(angles)]
    angles_notnan = np.unwrap(angles_notnan) * 180/np.pi
    
    new_angles = np.zeros(len(freqs))
    new_angles[np.isnan(angles)] = np.nan
    new_angles[~np.isnan(angles)] = angles_notnan
    angles = new_angles
    if angles[-1] - angles[0] < 0:
        p0 = -1
        bounds = ([-np.inf],[0])
    else:
        p0 = 1
        bounds = ([0],[np.inf])

    if "CS" in ant:
        p0 *= 0.01
    elif "St" in ant:
        p0 = 0
    elif "RS" in ant:
        p0 *= 0.1
    else:
        p0 *= 2
    
    if np.all(np.isnan(angles)):
        logging.warning(f"time interval {tint[0]}-{tint[1]} fully flagged for ant: {ant}")
        return np.nan, np.nan, np.array([]), np.array([])
    
    angles = angles - get_ref_angle(angles)
    freqs = freqs[np.isnan(angles) == False]
    angles = angles[np.isnan(angles) == False]

    if len(freqs) == 0:
        logging.warning(f"time interval {tint[0]}-{tint[1]} fully flagged for ant: {ant}")
        return np.nan, np.nan, np.array([]), np.array([])
    
    popt, pcov = curve_fit(dTEC, freqs, angles, p0=[p0], bounds=bounds)
    perr = np.sqrt(np.diag(pcov))
    if perr[0] == np.inf:
        logging.warning(f"fit failed for interval {tint[0]}-{tint[1]}, ant: {ant}")
        return np.nan, np.nan, np.array([]), np.array([])
    
    if plot:
        fig, axs = plt.subplots(1,2, figsize=(10,5))
        axs[0].scatter(freqs, angles, s=1)
        axs[0].plot(freqs, dTEC(freqs, popt), color='black')
        #axs[1].plot(stacked_freqs, stacked_phase, color='black',alpha=0.5, linewidth=1)
        axs[1].scatter(stacked_freqs, stacked_phase, s=1)
        plt.show()#; plt.cla(); plt.clf()
        
    return popt[0], perr[0], freqs, angles

def calculate_dtec(solspath, soltab="phase000", antenna="all", solution_interval=10, nstacks=None, split_band=True, constrain_phases=True, plot=False, t_start=None, select_time='all', mode='lombscargle'):
    phases, data = open_sols(solspath, soltab=soltab, constrain=constrain_phases)
    freq_split = int(len(data["freq"])/2) if split_band else 0 # exclude bottom half of band from fit (noisy)

    full_dtec = list()
    full_dtec_err = list()
    full_timesteps = list()
    full_ants = list()
    for i, ant in enumerate(data['ant']):
        if antenna != "all":
            if antenna == "dutch" and not ("CS" in ant or "RS" in ant or "St" in ant):
                continue
            elif antenna == "international" and ("CS" in ant or "RS" in ant or "St" in ant):
                continue
            elif antenna == "international_no_german" and ("CS" in ant or "RS" in ant or "DE" in ant or "St" in ant):
                continue
            elif antenna == "german" and "DE" not in ant:
                continue
        if ant != antenna and antenna not in ["all", "dutch", "international", "german", "international_no_german"]:
            continue
        
        logging.debug(f"working on ant: {ant}...")
        
        if nstacks is not None:
            timestack = nstacks
        elif 'RS' in ant or 'CS' in ant:
            timestack = 7
        else:
            timestack = 1
        
        if select_time == 'first':
            timesplit = np.arange(0, len(data['time'])//2, timestack)[:-1]
        elif select_time == 'second':
            timesplit = np.arange(len(data['time'])//2, len(data['time']), timestack)[:-1]
        else:
            timesplit = np.arange(0, len(data['time']), timestack)[:-1]
        
        dtec_arr = list()
        dtec_err_arr = list()
        time_steps = list()
        skipped = 0
        guess = 0
        for t, bin in enumerate(timesplit):
            if skipped < solution_interval:
                skipped += 1
                continue
            else:
                skipped = 0
            if t_start is not None:
                tint = [t_start, t_start+timestack]
            else:
                tint = [bin, bin+timestack]
            
            if mode == 'lombscargle':
                dtec, dtec_err, freqs, angles = fit_lombscargle(phases, data, tint, ant=ant, freq_split=freq_split, plot=plot, dtec_sign=guess)
                if abs(dtec) > 2:
                    guess = dtec/np.abs(dtec)
                else:
                    guess = 0
            elif mode == 'curvefit':
                dtec, dtec_err, freqs, angles = fit_unwrap(phases, data, tint, ant=ant, freq_split=freq_split, plot=plot)
            else:
                print("mode not recognized")
                sys.exit()
                
            dtec_arr.append(dtec)
            dtec_err_arr.append(dtec_err)
            time_steps.append(bin)

            if t_start is not None:
                break
        
        full_dtec.append(np.asarray(dtec_arr))
        full_dtec_err.append(np.asarray(dtec_err_arr))
        full_timesteps.append(np.asarray(time_steps))
        full_ants.append(ant)
    return full_dtec, full_dtec_err, full_timesteps, full_ants

def combine_dtec_tuples(tecs: list):
    full_dtec = list()
    full_dtec_err = list()
    full_timesteps = list()
    full_ants = list()
    for tec in tecs:
        full_dtec += tec[0]
        full_dtec_err += tec[1]
        full_timesteps += tec[2]
        full_ants += tec[3]
    return full_dtec, full_dtec_err, full_timesteps, full_ants

def get_smoothed_tec(dtec, timesteps, ant, new_time_axis):
    nonan_dtec = dtec.copy()[np.isnan(dtec) == False]
    t = timesteps[np.isnan(dtec) == False]
    lam = len(nonan_dtec)//20 
    if len(nonan_dtec) < 15:
        print("too few data points")
    elif lam%2 == 0:
        lam += 1

    filtered = medfilt(nonan_dtec, lam)

    if "PL" in ant or "LV" in ant or "IE" in ant:
        spl = make_smoothing_spline(t, filtered, lam=1e6)
    else:
        spl = make_smoothing_spline(t, filtered, lam=1e5)
    return spl(np.arange(len(new_time_axis)))


def write_tec_solutions(solspath, dtec, new_solspath: str = ""):
    import os
    if new_solspath == "":
        if "/" not in solspath:
            new_solspath = "cal-dtec.h5"
        else:
            new_solspath = os.path.dirname(solspath) + "/cal-dtec.h5"

    logging.debug(solspath)
    logging.debug(new_solspath)
    
    import os
    try: os.remove(new_solspath)
    except: pass
    os.system(f"cp {solspath} {new_solspath}")
    
    dtec_arr, __, __, ants = dtec
    dtec_vals = np.asarray([np.nanmedian(dtec_arr[i]) for i in range(len(ants))])
    dtec_dict = dict(zip(ants, dtec_vals))
    
    sols = losoto.h5parm(new_solspath, readonly = False)
    solset = sols.getSolset(solset = "sol000")
    soltab = solset.getSoltab(soltab = "phase000")
    
    __, data = soltab.getValues()
    
    aranged_dtec = np.zeros(len(data['ant']))
    for i, ant in enumerate(data['ant']):
        aranged_dtec[i] = dtec_dict[ant]
    
    tec = np.tile(aranged_dtec, (len(data['time']),1))

    weights = np.ones_like(tec)
    solset.makeSoltab(soltype='tec', soltabName='tec000', axesNames=['time', 'ant'], axesVals=[data['time'], data["ant"]], vals=tec, weights=weights)
    sols.close()
    

def write_smooth_tec_solutions(solspath, dtec, new_solspath: str = ""):
    import os
    if new_solspath == "":
        if "/" not in solspath:
            new_solspath = "cal-dtec.h5"
        else:
            new_solspath = os.path.dirname(solspath) + "/cal-dtec.h5"

    logging.debug(solspath)
    logging.debug(new_solspath)
    
    import os
    try: os.remove(new_solspath)
    except: pass
    os.system(f"cp {solspath} {new_solspath}")
    
    full_dtec, __, full_timesteps, full_ants = dtec
    
    sols = losoto.h5parm(new_solspath, readonly = False)
    solset = sols.getSolset(solset = "sol000")
    soltab = solset.getSoltab(soltab = "phase000")
    __, data = soltab.getValues()
    
    smoothed_dtec = np.zeros((len(data['time']), len(data['ant'])))
    print(full_ants)
    for i, ant in enumerate(data['ant']):
        idx = np.where(np.asarray(full_ants) == ant)[0][0]
        smoothed_dtec[:,i] = get_smoothed_tec(full_dtec[idx], full_timesteps[idx], full_ants[idx], data['time'])
    
    weights = np.ones_like(smoothed_dtec)
    solset.makeSoltab(soltype='tec', soltabName='tec000', axesNames=['time', 'ant'], axesVals=[data['time'], data["ant"]], vals=smoothed_dtec, weights=weights)
    sols.close()
    
  
    
def plot_dtec(dtec: tuple) -> None:
    full_dtec, full_dtec_err, full_timesteps, full_ants = dtec
    
    def weighted_mean(x, xerr):
        xerr[xerr == 0] = 0.2*x[xerr == 0]
        weight = 1./xerr**2
        weighted_mean = np.nansum(x*weight)/np.nansum(weight)
        weighted_error = np.sqrt(1./np.nansum(weight))
        return weighted_mean, weighted_error

    for i, ant in enumerate(full_ants):
        weighted_mean_dtec, weighted_error_dtec = weighted_mean(full_dtec[i], 50*full_dtec_err[i])
        plt.title(f"{ant}\n Weighted mean: {weighted_mean_dtec:.2f} +/- {weighted_error_dtec:.2f} \n Median: {np.nanmedian(full_dtec[i]):.2f}")
        plt.axhline(weighted_mean_dtec, color='black', linestyle='--')
        plt.axhline(np.nanmedian(full_dtec[i]), color='red', linestyle='--')
        plt.errorbar(full_timesteps[i], full_dtec[i], yerr=50*full_dtec_err[i], fmt="o", markersize=2, color='black', alpha=1, capsize=1, elinewidth=0.5)
        plt.show()
        print(np.nanmedian(full_dtec[i]))


def dTEC_fitter(solspath: str, mode: str = "tid") -> None:
    # get dtec of dutch stations at small intervals (every 5 datapoints) average over 7 datapoints
    dutch_dtec = calculate_dtec(solspath, antenna='dutch', solution_interval=5, nstacks=7, mode='curvefit', split_band=False)
    
    if mode == "tdt":
        # get dtec of german stations at large intervals (every 20 datapoints) average over 1 datapoint
        #  this run is to assess the quality of the data (what time part is better)
        german_dtec = calculate_dtec(solspath, antenna='german', solution_interval=20, nstacks=1, mode='curvefit')
        select_time = quality_check(*german_dtec)
    elif mode == "tid":
        select_time = 'all'

    # get dtec of international stations at large intervals (every 40 datapoints) average over 1 datapoint
    # use Lomb-Scargle to fit the data
    german_dtec = calculate_dtec(solspath, antenna='german', solution_interval=40, nstacks=1, mode='lombscargle')
    int_dtec = calculate_dtec(solspath, antenna='international_no_german', solution_interval=1, nstacks=1, select_time=select_time, mode='lombscargle') #international_no_german
    full_dtec = combine_dtec_tuples([dutch_dtec, german_dtec, int_dtec])
    if mode == "tdt":
        write_tec_solutions(solspath, full_dtec)
    elif mode == "tid":
        write_smooth_tec_solutions(solspath, full_dtec)
    



if __name__ == "__main__":
    import sys
    logging.info('Starting dTEC fitter')
    dTEC_fitter(sys.argv[1])
    logging.info('Done')
    