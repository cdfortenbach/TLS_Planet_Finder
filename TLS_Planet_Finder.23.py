#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Sat Sep 20 09:05:08 2020

@author: cfvm: Charles Fortenbach

Title: TLS_Planet_Finder

Version:  1.1.0

This project is licensed under the GNU GPLv3 License.

"""

import batman
import datetime
import lightkurve as lk
import matplotlib.pyplot as plt
import math
import numpy as np
import os
import pandas as pd
import scipy.optimize as op
import scipy.signal
import statistics
import textwrap
import time
from transitleastsquares import catalog_info
from transitleastsquares import cleaned_array
from transitleastsquares import transitleastsquares
from transitleastsquares import transit_mask


# Write summary table header to txt file in working directory
now = datetime.datetime.now()
timestr = time.strftime("%Y%m%d-%H%M%S")

fout = open(os.getcwd() + '/TLS_Smry_tables/' +'TLS_Smry_'
            + timestr + '.txt', "w", newline='')

fout.write('TLS_Planet_Finder run ID:  ' + str(now) + ' PT' + '\n' + '\n')

# Read in TLS config parameters:
print('\n')
print('TLS Planet Finder:')
print('\n')
print('Reading in TLS config parameters . . . ')
TLS_config = np.genfromtxt('TLS_config.txt', dtype=str, delimiter=', ', skip_header=6)  # all str
window_length     = int(TLS_config [0, 1])    # size of window sample pts for medfilt routine
ndist_factor      = int(TLS_config [1, 1])    # increasing this factor will reduce the distance between outlier peaks (troughs) identified
outlier_threshold = float(TLS_config [2, 1])  # remove transient if abs value of outlier depth is > outlier_threshold * the median of all other major troughs 
stdev_threshold   = float(TLS_config [3, 1])  # remove transient if the value of this std dev is > stdev_threshold * the overall std dev
sigma_upper       = int(TLS_config [4, 1])    # *** Note that the upper and lower bounds are reversed in Lightkurve (bug?)  
sigma_lower       = int(TLS_config [5, 1])    # *** Note that the upper and lower bounds are reversed in Lightkurve (bug?)
npps              = int(TLS_config [6, 1])    # planets to consider per system
osf               = int(TLS_config [7, 1])    # TLS oversampling factor
dgs               = float(TLS_config [8, 1])  # TLS transit grid step
tdm               = float(TLS_config [9, 1])  # TLS transit depth min

# Write config parameters to summary table
fout.write(' Configuration parameters for this run: ' + '\n')

fout.write('   window_length: ' 
           + str(window_length) + '\n')
fout.write('   ndist_factor: ' 
           + str(ndist_factor) + '\n')
fout.write('   outlier_threshold: ' 
           + str(outlier_threshold) + '\n')
fout.write('   stdev_threshold: ' 
           + str(stdev_threshold) + '\n')
fout.write('   sigma_upper: ' 
           + str(sigma_upper) + '\n')
fout.write('   sigma_lower: ' 
           + str(sigma_lower) + '\n')                  
fout.write('   planets per system limit: ' 
           + str(npps) + '\n') 
fout.write('   TLS oversampling factor: ' 
           + str(osf) + '\n') 
fout.write('   TLS transit grid step: ' 
           + str(dgs) + '\n') 
fout.write('   TLS transit depth min: ' 
           + str(tdm) + '\n' + '\n') 

# Define headers and formats for summary table output
head1 = ("|----------------------target info--------|-------------------------------------------------results----------------------------------------------------|\n")
head2 = (" Target          TOI    Msn   pc  Bitmask      R*      M*      Per   tdepth    tdur   edepth     inc    Rp/R*     Rp    a/R*    ecc   SDE    SNR  od/ev\n")
head3 = ("                                           (Rsun)  (Msun)      (d)    (ppm)    (hr)    (ppm)   (deg)            (Rj)                                (\u03C3)\n")
head4 = (" ------          ---    ---   --  -------    ----    ----    -----     ----    ----     ----    ----    -----   ----     ---    ---   ---    ---    ---\n")
fmt1 =  (" {a:15s} {b:5s} {c:>4s} {d:4d} {e:>8s} {f:7.2f} {g:7.2f} {h:8.4f} {i:8.0f} {j:7.3f} {k:8.0f} {l:7.1f} {m:8.3f} {n:6.2f}  {o:6.2f} {p:6.2f} {q:5.1f} {r:6.1f} {s:6.2f}\n")

fout.write(head1)
fout.write(head2)
fout.write(head3)
fout.write(head4)

# Define candidate planet names
names = ['nan', '1st', '2nd', '3rd', '4th']

# Define batman limb darkening model options
ld_opt = ['uniform', 'linear', 'quadratic', 'nonlinear', 'exponential', \
          'logarithmic', 'squareroot', 'power2', 'custom']
    
# Read in TOI list:
print('')
print('Reading in TOI list . . . ') 
TOI_data = np.loadtxt('TOI_list.txt', dtype=np.ndarray, skiprows=5)  # all str

# Solve the 1D TOI_data array problem
try:
    nshape = (TOI_data.shape[0], TOI_data.shape[1])
except IndexError:
    nshape = (1, TOI_data.shape[0]) 

# Start the TOI loop 
for n in range(0, nshape[0]): 
    if nshape[0] == 1:
        TOI = TOI_data[0]
    else:
        TOI = TOI_data[n, 0]
   
    # Now resolve the associated TIC
    url="https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=pipe"
    TOI_df = pd.read_csv(url, delimiter='|', index_col=1)
    TIC = TOI_df['TIC ID'][float(TOI)+.01]
    Target = 'TIC ' + str(TIC)
    
    if nshape[0] == 1:
        Mission = TOI_data[1]
        Bitmask = TOI_data[2]
    else:
        Mission = TOI_data[n, 1]
        Bitmask = TOI_data[n, 2]
   
    # Download the target light curve
    print('')
    print('')
    print('New Target:')
    print('Downloading light curve for the TOI '+ TOI + ' system . . .')
    print('This may take a while . . .')
    print('')
    search_result = lk.search_targetpixelfile(Target, mission=Mission)
    
    # Need to continue loop in the case of no data on MAST
    if len(search_result)==0:
        print('Data for this system not available on MAST')
        continue
    
    # Stitch available lc files
    tpfs = search_result.download_all(quality_bitmask=Bitmask)
    lcs = lk.LightCurveCollection([tpf.to_lightcurve() for tpf in tpfs])
    lc = lcs.stitch()
    print('')
    print(search_result)
    
    # Rename lc variables
    raw_time = lc.time
    raw_flux = lc.flux
 
    # Extract BTJD base from lc file
    BTJD_base = int(lc.astropy_time.jd[0] - lc.time[0])
  
    # Gather priors for stellar properties
    ab, M_star, M_star_min_err, M_star_max_err, R_star, R_star_min_err, R_star_max_err = catalog_info(TIC_ID=int(TIC))
    print('')
    print('Host star properties: ')
    print('  R_star: ' + str(format(R_star, '.3f')) + '  +' + str(format(R_star_max_err, '.3f')) + ' -' + str(format(R_star_min_err, '.3f')) + ' solar radii')
    print('  M_star: ' + str(format(M_star, '.3f')) + '  +' + str(format(M_star_max_err, '.3f')) + ' -' + str(format(M_star_min_err, '.3f')) + ' solar masses')
    print('  Limb darkening model: quadratic;  with coeffs:', ab)
    
    # Need to continue loop in the case of nan values for host star properties
    if math.isnan(R_star) == True:
        print('')
        print('Host star properties not available on MAST')
        continue
    
    # Start the planet loop
    for j in range(1, npps+1):
        
        pc = j  # jth planet candidate 
        
        # Plot raw (normalized) light curve
        fignum = npps*n +int(pc)
        plt.figure(fignum, figsize = (8, 10.5))
        ax1 = plt.subplot2grid((6, 2), (0, 0), rowspan=1, colspan=1) 
        ax1.xaxis.set_visible(False)
        ax1.set_xlim(lc.time.min(), lc.time.max())
        ax1.set_ylabel("Raw flux (rel)")
        ax1.plot(raw_time, raw_flux, "k.", markersize=1)
        
        # Clean up the light curve (C. Fortenbach's revised version of A. Mayo's cleanup functions)
        # Only need to run this once for each TIC:
        if j == 1:
            def iter_clip(lc, sigma_lower, sigma_upper):  
                new_lc, msk = lc.remove_outliers(sigma_lower, sigma_upper, return_mask=True)
                if sum(msk) == 0:
                    return new_lc
                else:
                    return iter_clip(new_lc, sigma_lower, sigma_upper)
            
            def my_custom_corrector_func(lc):
                corrected_lc = lc.normalize()
                clean_lc = iter_clip(corrected_lc, sigma_lower, sigma_upper)  
                return clean_lc
            
            lc = lcs.stitch(corrector_func=my_custom_corrector_func)

            t = lc.time  # New time array
            y = lc.flux  # New flux array
            
            # Determine typical cadence and other time stats
            delta_t = np.zeros(30)
            for ii in range(0, 30):
                delta_t[ii] = t[ii+1]-t[ii]
                delta_t[ii] = np.trunc(delta_t[ii]*10**7)/(10**7)
                
            delta_t_std = statistics.mode(delta_t) 
            total_t = (t[len(t)-1] - t[0])
            fraction_t = min(1, delta_t_std*len(t)/total_t)
            print('')
            print('Standard cadence: ' + str(format(delta_t_std*24*60, '.2f')) + ' (min)')
            print('Total time of dataset incl. gaps: ' + str(format(total_t, '.2f')) + ' (d)')
            print('Fraction of time covered by data: ' + str(format(fraction_t*100, '.2f')) + ' (%)')
            
            # Find time gaps/discontinuities in light curve
            n = 0
            tidx = np.zeros(100)
            for jj in range (1, len(t)-1):
                delta_t = t[jj+1] - t[jj]
                if (np.abs(delta_t - delta_t_std)) > 1.0:  # One day gap
                    n = n + 1
                    tidx[n] = jj
            tidx = tidx[tidx != 0]
            
            # Generate median trend line        
            trend = scipy.signal.medfilt(y, window_length)
            
            # Adjust medfilt trend line for discontinuities in time
            endspan = int(window_length/2)
            medspan = int(window_length/10)
            addspan = int(window_length/5)
            
            for n1 in range(0, len(tidx)):
                # Adjust left side of time gap
                t1 = int(tidx[n1]-endspan)
                t2 = int(tidx[n1]-addspan)
                y1 = trend[t1]
                y2 = np.median(y[int(t2-medspan):int(t2+medspan):]) 
                slope = (y2-y1)/(t2-t1)
                for n2 in range(t1, int(tidx[n1])):
                    trend[n2+1] = trend[n2] + slope

                # Adjust right side of time gap    
                t1 = int(tidx[n1]+1+addspan)
                t2 = t1 + endspan
                y1 = np.median(y[int(t1-medspan):int(t1+medspan)]) 
                y2 = trend[t1 + endspan] 
                slope = (y2-y1)/(t2-t1)
                for n3 in range(0, endspan-1):
                    trend[int(tidx[n1]+endspan-n3-1)] = trend[int(tidx[n1]+endspan-n3)]-slope
            
            # Clean up lc endpts disturbed by medfilt
            # For the left end
            t1 = medspan
            t2 = endspan
            y1 = np.median(y[t1-medspan:t1+medspan]) 
            y2 = trend[endspan] 
            slope = (y2-y1)/(t2-t1)  
            for n4 in range(0, endspan):
                trend[endspan-n4-1] = trend[endspan-n4] - slope 
                
            # For the right end
            t1 = len(t)-endspan
            t2 = len(t)-medspan
            y1 = trend[t1]
            y2 = np.median(y[t2-medspan:t2+medspan]) 
            slope = (y2-y1)/(t2-t1)
            for n5 in range(len(t)-endspan+1, len(t)-1):
                trend[n5+1] = trend[n5] + slope
                
            # Recast flux with adjusted median filter trend
            y = y /trend
            
            # Find transit minimums ('peaks')
            ndist = int(len(t)/ndist_factor)
            peaks, properties = scipy.signal.find_peaks(-y, height=None, distance=ndist)
            yorder = np.sort(y[peaks])
              
            # Remove anomalous transient outliers
            for kk in range(0, len(yorder)-1):
                if (1-yorder[kk]) > outlier_threshold*(1-np.median(yorder)):
                    ymin_idx = np.where(y == yorder[kk])
                    nymin_idx = int(ymin_idx[0][0])
                    # now deal in time; transient -1/+3 days
                    t_lbound = float(t[ymin_idx] - 1.5)  # days
                    lbound_idx = (np.abs(t - t_lbound)).argmin()
                    t_rbound = float(t[ymin_idx] + 1.5)  # days
                    rbound_idx = (np.abs(t - t_rbound)).argmin()
                    nsize = rbound_idx - lbound_idx
                    scale = (np.std(y[lbound_idx-60:lbound_idx]) + np.std(y[rbound_idx:rbound_idx+60]))/2 
                    y[lbound_idx:rbound_idx] =  np.random.normal(loc=1, scale=scale, size=nsize)  # fill with white noise
                    print('\n')
                    print('** removed anomalous transient event from stitched/filtered light curve **')
                    
            # Find std dev outliers
            nspan = int(window_length/2)
            y_stdev = np.zeros(len(t))
            for nz in range(0, len(t)):
                y_stdev[nz] = np.std(y[nz:nz+nspan])
            
            ndist = int(len(t)/ndist_factor)
            stdev_peaks, properties = scipy.signal.find_peaks(y_stdev, height=None, distance=ndist)
            y_stdev_order = -np.sort(-y_stdev[stdev_peaks])
            stdev = np.std(y)
            
            # Remove anomalous std dev outliers
            for kk in range(0, len(y_stdev_order)-1):
                if (y_stdev_order[kk]) > stdev_threshold*(np.median(y_stdev_order)):
                    ymin_idx = np.where(y_stdev == y_stdev_order[kk])
                    nymin_idx = int(ymin_idx[0][0])
                    # now deal in time; transient -1/+3 days
                    t_lbound = float(t[ymin_idx] - 1.5)  # days
                    lbound_idx = (np.abs(t - t_lbound)).argmin()
                    t_rbound = float(t[ymin_idx] + 1.5)  # days
                    rbound_idx = (np.abs(t - t_rbound)).argmin()
                    nsize = rbound_idx - lbound_idx
                    scale = (np.std(y[lbound_idx-60:lbound_idx]) + np.std(y[rbound_idx:rbound_idx+60]))/2 
                    y[lbound_idx:rbound_idx] =  np.random.normal(loc=1, scale=scale, size=nsize)  # fill with white noise
                    print('\n')
                    print('** removed anomalous increased std dev event from stitched/filtered light curve **')
                        
            # Initialize certain parameters
            tls_results = 0
            t_dur_days = 0
   
        else:          
            # Now, use the transit_mask function to mask the in-transit points of each planet in turn.
            intransit = transit_mask(t, tls_results.period, 2*t_dur_days, tls_results.T0)
            y = y[~intransit]
            t = t[~intransit]
            t, y = cleaned_array(t, y)
          
        # Plot detrended and masked light curve
        ax2 = plt.subplot2grid((6, 2), (1, 0), rowspan=1, colspan=1)
        ax2.set_xlim(t.min(), t.max())
        ax2.set_xlabel('Time [BTJD - ' + str(BTJD_base) +']')
        ax2.set_ylabel('Flux (ppm)')
        if j != 1:
            ax2.text(0.75, 0.95, 'Masked', transform=ax2.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=0.8))
        ax2.text(0.75, 0.15, 'Detrended', transform=ax2.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=0.8))
        stdev = np.std((y - 1.0)*10**6)  # Convert stdev to ppm units
        ax2.text(0.03, 0.15, '$\sigma = $' +str(format(stdev, '.0f')) + '$\,$ppm', color='red', transform=ax2.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=0.8))
        ax2.plot(t, (y - 1.0)*10**6, "k.", markersize=1)  # Convert flux to ppm units
        
        # Now run TLS on the detrended and masked light curve
        print(' ') 
        print('Search for the', names[j], 'planet using TLS:') 
        time.sleep(0.5)
        model = transitleastsquares(t, y)
        tls_results = model.power(u=ab, oversampling_factor=osf,\
                              duration_grid_step=dgs, transit_depth_min=tdm)
        
        t_folded = (tls_results.folded_phase - 0.5)*tls_results.period
        flux_data = tls_results.folded_y           #normalized folded flux
        flux_err = tls_results.folded_dy           #1-sigma error on flux

        # Set best guess batman parameters for this target
        t0_guess  = 0.                             #time of inferior conjunction (days)
        per_guess = tls_results.period             #orbital period (days), from TLS
        rp_guess  = np.sqrt(1-tls_results.depth)   #planet radius (units of stellar rad), from TLS
        a_guess   = (((per_guess/365.2422)**2)*M_star)**(0.3333333)*215.032/R_star #semi-major axis (units of stellar rad)
        b_guess     = 0.5                          #guess for value of transit impact parameter: b
        inc_guess = np.arccos(b_guess/(a_guess))*180/np.pi #orbital inclination (in degrees)
        ecc_guess = 0.0                            #eccentricity
        w_guess   = 0                              #longitude of periastron (in degrees)  
        f_guess   = 0.00001                        #data variance may be underestimated by a fractional amount: f    
        ld_model  = ld_opt[2]                      #ld model options
        ld_coeffs = ab                             #ld model coefficients from stellar parameters
        bman_int_fac = 0                           #scale factor for batman integration step size (0 default/optimal)
        
        # Need to break planet candidate loop for the case that no transits are fit (per_guess will be set to nan by TLS)
        if math.isnan(per_guess) == True:
            print('')
            print('No detection of additional planets in this system beyond the', names[j-1], 'candidate.')
            break
        
        # Define log likelihood function
        def lnlike(theta, x, y, yerr):
            t0, per, rp, a, inc, ecc, w, lnf = theta
            params = batman.TransitParams()     
            params.t0 = t0                      
            params.per = per                    
            params.rp = rp                      
            params.a = a                        
            params.inc = inc                    
            params.ecc = ecc                    
            params.w = w                        
            params.limb_dark = ld_model     
            params.u = ld_coeffs  
            m = batman.TransitModel(params, x, fac=bman_int_fac) 
            model = m.light_curve(params)          
            inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
            return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))
        
        # Determine chi2
        chi2 = lambda *args: -1 * lnlike(*args)
        
        # Define bounds on free parameters
        t0_bnds  = (-0.0000001, 0.0000001)
        per_bnds = (0.99995*per_guess, 1.00005*per_guess)
        rp_bnds  = (0.950*rp_guess, 1.050*rp_guess)
        a_bnds   = (0.950*a_guess, 1.050*a_guess)
        b_max    = 1  #max value of transit impact parameter
        inc_bnds = (np.arccos(b_max/(0.950*a_guess))*180/np.pi, 90)
        ecc_bnds = (0.0, 0.2) 
        w_bnds   = (-120, 120)
        lnf_bnds = (-10, 1)
        
        # Determine max likelihood (min chi2) solution for free parameters
        result = op.minimize(chi2, [t0_guess, per_guess, rp_guess, a_guess, inc_guess, ecc_guess, w_guess, np.log(f_guess)], args=(t_folded, flux_data, flux_err), \
                             bounds=[t0_bnds, per_bnds, rp_bnds, a_bnds, inc_bnds, ecc_bnds, w_bnds, lnf_bnds])     
        t0_ml, per_ml, rp_ml, a_ml, inc_ml, ecc_ml, w_ml, lnf_ml = result["x"]
               
        print('\n')
        print("""Maximum likelihood result (primary transit):
              t0  = {0:9.5f} (guess:{1:9.5f}) days
              per = {2:9.5f} (guess:{3:9.5f}) days
              rp  = {4:9.5f} (guess:{5:9.5f}) stellar radii
              a   = {6:9.5f} (guess:{7:9.5f}) stellar radii
              inc = {8:9.5f} (guess:{9:9.5f}) deg
              ecc = {10:9.5f} (guess:{11:9.5f})
              w   = {12:9.5f} (guess:{13:9.5f}) deg
              f   = {14:9.5f} (guess:{15:9.5f})
        """.format(t0_ml, t0_guess, per_ml, per_guess, rp_ml, rp_guess, a_ml, a_guess, inc_ml, inc_guess, ecc_ml, ecc_guess, w_ml, w_guess, np.exp(lnf_ml), f_guess))
        
        # Reset batman parameters for ml case
        params_ml = batman.TransitParams()         
        params_ml.t0 = t0_ml                        
        params_ml.per = per_ml                      
        params_ml.rp = rp_ml                        
        params_ml.a = a_ml                          
        params_ml.inc = inc_ml                      
        params_ml.ecc = ecc_ml                      
        params_ml.w = w_ml                          
        params_ml.limb_dark = ld_model              
        params_ml.u = ld_coeffs               
        
        # Run batman transit model for ml params
        ml = batman.TransitModel(params_ml, t_folded, fac=bman_int_fac) 
        
        # Estimate the transit duration
        transit_start_idx = 0
        for i in range(0, len(t_folded)):
            if ml.light_curve(params_ml)[i] < 1.0:
                transit_start_idx = i
                transit_start_time = t_folded[i]
                break
        for i in range(transit_start_idx+1, len(t_folded)):
            if ml.light_curve(params_ml)[i] == 1.0:
                transit_end_time = t_folded[i]
                break
        t_dur_days = transit_end_time - transit_start_time   # est. of transit duration, T14, in days
        t_dur_days = max(0.0000001, t_dur_days)
        print('              t_dur (hrs)     =', str(format(t_dur_days*24, '.2f')))
        
        # Convert flux values to ppm
        flux_data_ppm = (flux_data-1)*10**6
        flux_ml_ppm = (ml.light_curve(params_ml)-1)*10**6
        flux_err_ppm = flux_err*10**6
        
        # Estimate the transit depth
        t_depth_ppm = -np.amin(flux_ml_ppm)   # est. of transit depth, norm
        t_depth_ppm = max(0.0000001, t_depth_ppm)
        print('              t_depth (ppm)   =', str(format(t_depth_ppm, '.0f')))
        
        # Plot phase folded transit light curve vs time-from-mid-transit
        ax3 = plt.subplot2grid((6, 2), (2, 0), rowspan=1, colspan=1)
        nplot_slice = max(1, int(len(t_folded)/20000))   # reduce density of plot data (calcs continue with full data-set)
        ax3.scatter(t_folded[::nplot_slice], flux_data_ppm[::nplot_slice], color='darkgrey', s=1, alpha=0.5, zorder=2) 
        ax3.plot(t_folded, flux_ml_ppm, "r--", lw=1.5, alpha=1.0)
        ax3.set_xlim(-t_dur_days*1.5, t_dur_days*1.5)
        ax3.set_ylim(min(-0.01, np.amin(flux_ml_ppm)*1.50), max(0.01, -np.amin(flux_ml_ppm)*0.50))
        ax3.set_xlabel("Time from mid-transit (days)")
        ax3.set_ylabel("Flux (ppm)")
        ax3.legend(['model fit'], loc='lower right', fontsize=9)
        ax3.text(0.03, 0.30, 'Primary', transform=ax3.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=1.0))
        ax3.text(0.03, 0.15, str(names[j]) + ' planet', transform=ax3.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=1.0))
        
        # Now consider the secondary eclipse
        # Need to compound the time and flux arrays to access phase 0.5
        t_cmpd_se = np.append(t_folded, t_folded + tls_results.period) - 0.5*tls_results.period
        flux_cmpd = np.append(flux_data, flux_data)
        flux_err_cmpd = np.append(flux_err, flux_err)
        
        # Convert flux values to ppm
        flux_cmpd_ppm = (flux_cmpd-1)*10**6
        flux_err_cmpd_ppm = np.append(flux_err, flux_err)*10**6
        
        # Set best guess parameters for secondary eclipse case (reflected light only)                         
        albedo = 0.45
        phi_alpha = 1.0
        fp_guess = albedo*((rp_ml/a_ml)**2)*(phi_alpha)
        t_secondary_guess = 0.
        
        # Define log likelihood function (secondary version)
        def lnlike_sec(theta, x, y, yerr):
            t0, per, rp, a, inc, ecc, w, lnf, fp, t_secondary = theta
            params = batman.TransitParams()     
            params.t0 = t0                      
            params.per = per                    
            params.rp = rp                      
            params.a = a                        
            params.inc = inc                    
            params.ecc = ecc                    
            params.w = w   
            params.fp = fp
            params.t_secondary = t_secondary            
            params.limb_dark = ld_model     
            params.u = ld_coeffs 
            m = batman.TransitModel(params, x, fac=bman_int_fac, transittype="secondary") 
            model = m.light_curve(params)          
            inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
            return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))
        
        # Determine chi2
        chi2 = lambda *args: -1 * lnlike_sec(*args)
        
        # Define bounds on free parameters
        t0_bnds  = (-0.0000001, 0.0000001)
        per_bnds = (per_ml, 1.000001*per_ml)
        rp_bnds  = (rp_ml, 1.000001*rp_ml)
        a_bnds   = (a_ml, 1.000001*a_ml)
        inc_bnds = (0.99999*inc_ml, inc_ml)
        ecc_bnds = (ecc_ml-0.0000001, ecc_ml+0.0000001) 
        w_bnds   = (w_ml-0.0000001, w_ml+0.0000001)
        lnf_bnds = (lnf_ml-0.0000001, lnf_ml+0.0000001)
        fp_bnds  = (0.45*fp_guess, 1.75*fp_guess)
        t_secondary_bnds = (-0.0000001, 0.0000001)
        
        # Determine the max likelihood (min chi2) solution for free parameters
        result = op.minimize(chi2, [t0_ml, per_ml, rp_ml, a_ml, inc_ml, ecc_ml, w_ml, lnf_ml, fp_guess, t_secondary_guess], args=(t_folded, flux_data, flux_err), \
                             bounds=[t0_bnds, per_bnds, rp_bnds, a_bnds, inc_bnds, ecc_bnds, w_bnds, lnf_bnds, fp_bnds, t_secondary_bnds])     
        t0_se, per_se, rp_se, a_se, inc_se, ecc_se, w_se, lnf_se, fp_se, t_secondary_se = result["x"]
        
        print('\n')
        print("""Maximum likelihood result (secondary eclipse):
              fp  = {0:9.5f} (guess:{1:9.5f})
              ts  = {2:9.5f} (guess:{3:9.5f}) days
        """.format(fp_se, fp_guess, t_secondary_se, t_secondary_guess))
        
        # Reset the batman parameters for se case
        params_se = batman.TransitParams()         
        params_se.t0 = t0_se                        
        params_se.per = per_se                      
        params_se.rp = rp_se                        
        params_se.a = a_se                          
        params_se.inc = inc_se                      
        params_se.ecc = ecc_se                      
        params_se.w = w_se                          
        params_se.limb_dark = ld_model              
        params_se.u = ld_coeffs
        params_se.fp = fp_se
        params_se.t_secondary = t_secondary_se
        
        # Run the batman transit model for secondary eclipse params, ms
        ms = batman.TransitModel(params_se, t_cmpd_se, transittype="secondary")
        flux_se_ppm = (ms.light_curve(params_se)-1)*10**6
        edepth = np.amax(flux_se_ppm)
        
        # Estimate the eclipse depth
        print('              ecl_depth (ppm) =', str(format(edepth, '.0f')))
        
        # Plot the phase folded transit light curve vs time-from-mid-eclipse
        ax4 = plt.subplot2grid((6, 2), (2, 1), rowspan=1, colspan=1) 
        nplot_slice = 1   # reduce density of plot data (calcs continue with full data-set)
        ax4.scatter(t_cmpd_se[::nplot_slice], flux_cmpd_ppm[::nplot_slice], color='darkgrey', s=0.5, alpha=0.8, zorder=2)
        ax4.plot(t_cmpd_se, flux_se_ppm, "b--", lw=1.5, alpha=1.0)
        ax4.set_xlim(-t_dur_days*1.5, t_dur_days*1.5)
        ax4.set_ylim(min(-200, -np.amax(flux_se_ppm)*3.25), max(200, np.amax(flux_se_ppm)*2.75))
        ax4.set_xlabel("Time from central eclipse (days)")
        ax4.legend(['model fit'], loc='lower right', fontsize=9)
        
        # Plot the phase folded transit light curve vs phase
        ax5 = plt.subplot2grid((6, 2), (3, 0), rowspan=1, colspan=2)
        ax5.plot((tls_results.folded_phase[0:int(0.9*len(tls_results.folded_phase))] - 0.5), flux_ml_ppm[0:int(0.9*len(tls_results.folded_phase))], 'r--')
        nplot_slice = max(1, int(len(t_folded)/20000))   # reduce density of plot data (calcs continue with full data-set)
        ax5.scatter((tls_results.folded_phase[::nplot_slice] - 0.5), (tls_results.folded_y[::nplot_slice] - 1.0)*10**6, color='darkgrey', s=0.5, alpha=0.5, zorder=2)
        ax5.scatter((tls_results.folded_phase[::nplot_slice] + 0.5), (tls_results.folded_y[::nplot_slice] - 1.0)*10**6, color='darkgrey', s=0.5, alpha=0.5, zorder=2)
        ax5.set_xlim(-0.2, 0.8)
        ax5.set_ylim(min(-0.01, np.amin(flux_ml_ppm)*1.50), max(0.01, -np.amin(flux_ml_ppm)*0.60))
        ax5.set_xlabel('Phase')
        ax5.set_ylabel('Flux (ppm)')
        ax5.text(0.01, 0.15, str(names[j]) + ' planet', transform=ax5.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=1.0))
           
        # Plot the TLS periodogram:
        ax6 = plt.subplot2grid((6, 2), (5, 0), rowspan=1, colspan=2)
        ax6.axvline(tls_results.period, alpha=0.4, lw=3)
        ax6.set_xlim(np.min(tls_results.periods), np.max(tls_results.periods))
        # Highlight the harmonics of the peak period
        for n in range(2, 10):
            ax6.axvline(n*tls_results.period, alpha=0.4, lw=1, linestyle="dashed")
            ax6.axvline(tls_results.period / n, alpha=0.4, lw=1, linestyle="dashed")
        ax6.set_ylabel(r'SDE')
        ax6.set_xlabel('Period (days)')
        ax6.text(0.60, 0.85,'TLS Periodogram, ' + str(names[j]) + ' planet', transform=ax6.transAxes, fontsize=9)
        ax6.text(0.65, 0.70,'Period (d) = ' + str(format(tls_results.period, '.4f')) + ';  SDE = ' + str(format(tls_results.SDE, '.1f')), transform=ax6.transAxes, fontsize=9)
        ax6.plot(tls_results.periods, tls_results.power, color='black', lw=0.5)
        ax6.set_xlim(0, max(tls_results.periods))
                
        # Plot the full light curve, together with the best transit model (in red)
        ax7 = plt.subplot2grid((6, 2), (4, 0), rowspan=1, colspan=2)
        in_transit = transit_mask(
            t,
            tls_results.period,
            t_dur_days,
            tls_results.T0)
        ax7.scatter(
            t[in_transit],
            (y[in_transit] - 1.0)*10**6,
            color='red',
            s=2,
            zorder=0)
        ax7.scatter(
            t[~in_transit],
            (y[~in_transit] - 1.0)*10**6,
            color='darkgrey',
            alpha=0.5,
            s=1,
            zorder=0)
        ax7.plot(
            tls_results.model_lightcurve_time,
            (tls_results.model_lightcurve_model - 1.0)*10**6, alpha=0.5, color='red', zorder=1)
        ax7.set_xlim(min(lc.time), max(lc.time))
        ax7.set_xlabel('Time [BTJD - ' + str(BTJD_base) +']')
        ax7.set_ylabel('Flux (ppm)')
        ax7.text(0.01, 0.15, ' Full light curve', transform=ax7.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=0.8))
        ax7.text(0.81, 0.15, str(names[j]) + ' planet transits', transform=ax7.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=0.8))
          
        # Print some statistics to the console:
        print(' ')     
        print('Signal detection eff.(SDE):', format(tls_results.SDE,'.1f'))
        print('SNR:', format(tls_results.snr,'.1f'))
        print(' ') 
          
        # Create subplot with stats only, no axes
        ax8 = plt.subplot2grid((6, 2), (0, 1), rowspan=2, colspan=1)
        ax8.axis('off')
        ax8.text(0.10,0.96, Target +'; TOI-' + TOI + '.xx', transform=ax8.transAxes)
        ax8.text(0.13,0.88, 'Bitmask = ' + str(Bitmask), transform=ax8.transAxes, fontsize=9)
        ax8.text(0.13,0.80, 'Planet candidate = ' + str(pc), transform=ax8.transAxes, fontsize=9)
        ax8.text(0.13,0.72, 'Period (days) = ' + str(format(per_ml, '.5f')), transform=ax8.transAxes, fontsize=9)
        ax8.text(0.13,0.64, 'Transit depth (ppm) = ' + str(format(t_depth_ppm, '.0f')), transform=ax8.transAxes, fontsize=9)
        ax8.text(0.13,0.56, 'Transit duration (hr) = ' + str(format(t_dur_days*24, '.4f')), transform=ax8.transAxes, fontsize=9)
        ax8.text(0.13,0.48, 'SDE = ' + str(format(tls_results.SDE,'.1f')) + ';  SNR = ' + str(format(tls_results.snr,'.1f')), transform=ax8.transAxes, fontsize=9)
        ax8.text(0.13,0.40, 'Eclipse depth (ppm) = ' + str(format(edepth, '.0f')), transform=ax8.transAxes, fontsize=9)
        ax8.text(0.13,0.32, '$R_{*}$ = ' + str(format(R_star, '.2f')) + ' (+' + str(format(R_star_max_err, '.2f')) + ', -' + str(format(R_star_min_err, '.2f')) +') $R_{\odot}$', transform=ax8.transAxes, fontsize=9)
        ax8.text(0.13,0.24, '$M_{*}$ = ' + str(format(M_star, '.2f')) + ' (+' + str(format(M_star_max_err, '.2f')) + ', -' + str(format(M_star_min_err, '.2f')) +') $M_{\odot}$', transform=ax8.transAxes, fontsize=9)
        ax8.text(0.13,0.16, '$R_{p}/R_{*}$ = ' + str(format(rp_ml, '.3f')) + ';  $a/R_{*}$ = ' + str(format(a_ml, '.2f')), transform=ax8.transAxes, fontsize=9)
        ax8.text(0.13,0.08, '$R_{p}$ = ' + str(format(rp_ml*R_star*109.2, '.2f')) + '$\,R_{\oplus}$' + ' = ' + str(format(rp_ml*R_star*10.973, '.2f')) + '$\,R_{jup}$', transform=ax8.transAxes, fontsize=9)
        ax8.text(0.13,0.00, 'odd/even mismatch = ' + str(format(tls_results.odd_even_mismatch, '.2f')) + ' $\sigma$', transform=ax8.transAxes, fontsize=9)
        
        # Make certain table variables explicit
        Rp_Rs = rp_ml
        Rp_j = rp_ml*R_star*10.973
        odev = tls_results.odd_even_mismatch
        
        # Save plots and results to combined vetting sheet for this planet candidate
        plt.tight_layout()   
        fname = [os.getcwd() + '/TLS_Vetting_sheets/' + Target + '_pc_' + str(pc) +'.svg']
        plt.savefig(fname[0], overwrite=True, bbox_inches='tight')
        plt.close('all')
        
        # Write data/results to the summary table for this target/planet
        fout.write(fmt1.format(                       
                     a = Target,
                     b = TOI,
                     c = Mission,
                     d = pc, 
                     e = Bitmask,
                     f = R_star,
                     g = M_star,
                     h = per_ml,
                     i = t_depth_ppm,
                     j = t_dur_days*24,
                     k = edepth,
                     l = inc_ml,
                     m = Rp_Rs,
                     n = Rp_j,
                     o = a_ml,
                     p = ecc_ml,
                     q = tls_results.SDE,
                     r = tls_results.snr,
                     s = odev))
        
# Write final note to summary table
fout.write(textwrap.dedent("""\
 
_____________________________                   
Notes:
    For a detection to be significant, we need a minimum S/N of 7.1 (Heller et al. 2020), 
    which was shown to yield only one false positive threshold-crossing event over the 
    entire Kepler mission (Jenkins et al. 2002). For the SDE_TLS value we apply a minimum
    detection threshold of 9, which results in a false-positive rate < 10^{-4} in the
    limiting case of white noise (Hippke & Heller 2019).
       
 """))    

fout.close()

raise SystemExit