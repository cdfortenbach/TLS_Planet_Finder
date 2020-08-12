#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 07:34:08 2020

@author: cfvm: Charles Fortenbach

Title: TLS_Planet_Finder

Version:  1.0.0

This project is licensed under the GNU GPLv3 License.

"""

import numpy as np
import matplotlib.pyplot as plt
import os
import datetime
import time
import textwrap
import pandas as pd
import lightkurve as lk
from transitleastsquares import transitleastsquares
from transitleastsquares import transit_mask
from transitleastsquares import cleaned_array
from transitleastsquares import catalog_info
import statsmodels.api as sm
lowess = sm.nonparametric.lowess


# Write summary table header to txt file in working directory
now = datetime.datetime.now()
timestr = time.strftime("%Y%m%d-%H%M%S")

fout = open(os.getcwd() + '/TLS_Smry_tables/' +'TLS_Smry_'
            + timestr + '.txt', "w", newline='')

fout.write('TLS_Planet_Finder run ID:  ' + str(now) + ' PT' + '\n' + '\n')

# Read in TLS config parameters:
print('\n') 
print('Reading in TLS config parameters . . . ')
TLS_config = np.genfromtxt('TLS_config.txt', dtype=str, delimiter=', ', skip_header=6)  # all str
window_length = int(TLS_config [0, 1])          # for lc corrector: A. Mayo used window_length=501
polyorder = int(TLS_config [1, 1])              # for lc corrector: A. Mayo used default 2 
sigma_upper = int(TLS_config [2, 1])            # ***Note that the upper and lower bounds are reversed in Lightkurve (bug?)  
sigma_lower = int(TLS_config [3, 1])            # ***Note that the upper and lower bounds are reversed in Lightkurve (bug?)
npps = int(TLS_config [4, 1])                   # planets to consider per system
lowess_frac_prim = float(TLS_config [5, 1])     # lowess filter factor for primary transit smoothing curve
flux_lim_prim = float(TLS_config [6, 1])        # flux limiter for primary transit smoothing curve
lowess_frac_sec = float(TLS_config [7, 1])      # lowess filter factor for secondary eclipse smoothing curve
flux_lim_sec = float(TLS_config [8, 1])         # flux limiter for secondary eclipse smoothing curve
osf = int(TLS_config [9, 1])                    # TLS oversampling factor
dgs = float(TLS_config [10, 1])                 # TLS transit grid step
tdm = float(TLS_config [11, 1])                 # TLS transit depth min

# Write config parameters to summary table
fout.write(' Configuration parameters for this run: ' + '\n')
fout.write('   window_length: ' 
           + str(window_length) + '\n')
fout.write('   polyorder: ' 
           + str(polyorder) + '\n')
fout.write('   sigma_upper: ' 
           + str(sigma_upper) + '\n')
fout.write('   sigma_lower: ' 
           + str(sigma_lower) + '\n')                  
fout.write('   planets per system limit: ' 
           + str(npps) + '\n') 
fout.write('   lowess_frac_prim: ' 
           + str(lowess_frac_prim) + '\n') 
fout.write('   flux_lim_prim: ' 
           + str(flux_lim_prim) + '\n') 
fout.write('   lowess_frac_sec: ' 
           + str(lowess_frac_sec) + '\n') 
fout.write('   flux_lim_sec: ' 
           + str(flux_lim_sec) + '\n') 
fout.write('   TLS oversampling factor: ' 
           + str(osf) + '\n') 
fout.write('   TLS transit grid step: ' 
           + str(dgs) + '\n') 
fout.write('   TLS transit depth min: ' 
           + str(tdm) + '\n' + '\n') 

# Define headers and formats for summary table output
head1 = ("|--------------target info-----------------|-----------------------------------results------------------------------------|\n")
head2 = (" Target          TOI    Msn   pc   Bitmask      Per   tdepth     tdur    SDE    SNR      R*      M*   Rp/R*     Rp   od/ev\n")
head3 = ("                                                (d)    (ppm)     (hr)                (Rsun)  (Msun) 	      (Rj)   (sig)\n")
head4 = (" ------          ---    ---   --   -------      ---   ------   ------   ----  -----    ----    ----   -----   ----     ---\n")
fmt1 =  (" {Target:15s} {TOI:5s} {Mission:>4s} {pc:4d} {Bitmask:>9s} {Period:8.4f} {tdepth:8.0f} {tdur:8.4f} {SDE:6.1f} {snr:6.1f} {R_star:7.2f} {M_star:7.2f} {Rp_Rs:7.3f} {Rp_j:6.2f} {odev:7.2f} \n")

fout.write(head1)
fout.write(head2)
fout.write(head3)
fout.write(head4)

# Define candidate planet names
names = ['nan', '1st', '2nd', '3rd', '4th']

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
    print('Downloading the light curve for TOI '+ TOI + ' . . .')
    print('This may take a while . . .')
    print('')
    search_result = lk.search_targetpixelfile(Target, mission=Mission)
    tpfs = search_result.download_all(quality_bitmask=Bitmask)
    lcs = lk.LightCurveCollection([tpf.to_lightcurve() for tpf in tpfs])
    lc = lcs.stitch()
    print('')
    print(search_result)
    raw_time = lc.time
    raw_flux = lc.flux
      
    # Extract BTJD base from lc file
    BTJD_base = int(lc.astropy_time.jd[0] - lc.time[0])
  
    # Gather priors for stellar properties
    ab, M_star, M_star_min_err, M_star_max_err, R_star, R_star_min_err, R_star_max_err = catalog_info(TIC_ID=int(TIC))
    
    # Start the planet loop
    for j in range(1, npps+1):
        
        pc = j  # jth planet candidate 
        
        # Plot raw light curve
        fignum = npps*n +int(pc)
        plt.figure(fignum, figsize = (8, 10.5))
        ax1 = plt.subplot2grid((6, 2), (0, 0), rowspan=1, colspan=1) 
        ax1 = plt.gca()
        ax1.xaxis.set_visible(False)
        ax1.set_ylabel("Raw flux")
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
                corrected_lc = lc.normalize().flatten(window_length, polyorder)
                clean_lc = iter_clip(corrected_lc, sigma_lower, sigma_upper)  
                return clean_lc
            
            lc = lcs.stitch(corrector_func=my_custom_corrector_func)
            
            t = lc.time
            y = lc.flux
            results = 0
            t_dur = 0
   
        else:          
            # Now, we use the transit_mask function to mask the in-transit points of each planet in turn.
            intransit = transit_mask(t, results.period, 2*t_dur, results.T0)
            y = y[~intransit]
            t = t[~intransit]
            t, y = cleaned_array(t, y)
          
        # Plot detrended and masked light curve
        ax2 = plt.subplot2grid((6, 2), (1, 0), rowspan=1, colspan=1)
        ax2 = plt.gca()
        ax2.set_xlim(lc.time.min(), lc.time.max())
        ax2.set_xlabel('Time [BTJD - ' + str(BTJD_base) +']')
        ax2.set_ylabel('Flux (ppm)')
        ax2.text(0.75, 0.15, 'Detrended', transform=ax2.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=0.8))
        stdev = np.std((y - 1.0)*10**6)  # Convert stdev to ppm units
        ax2.text(0.03, 0.15, '$\sigma = $' +str(format(stdev, '.0f')) + '$\,$ppm', color='red', transform=ax2.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=0.8))
        ax2.plot(t, (y - 1.0)*10**6, "k.", markersize=1)  # Convert flux to ppm units
        
        # Now run TLS on the detrended and masked light curve
        print(' ') 
        print(' ') 
        print(' ') 
        print('Search for the', names[j], 'planet using TLS:') 
        model = transitleastsquares(t, y)
        results = model.power(u=ab, oversampling_factor=osf,\
                              duration_grid_step=dgs, transit_depth_min=tdm)
        
        # Now use lowess to smooth and recover signal in noise for primary transit
        lowess_frac = lowess_frac_prim     # key toggle parameter for lowess filter
        filtered_primary = lowess(results.folded_y, (results.folded_phase - 0.5)*results.period, lowess_frac)
        
        # if flux level is greater than flux_lim_prim, then force to 1.0
        for i in range(0, len(filtered_primary[:, 1])):
            if filtered_primary[i, 1] > flux_lim_prim:
               filtered_primary[i, 1] = 1.0
        
        # Estimate transit depth and duration from filtered, phase-folded curve
        for i in range(int(0.4*len(filtered_primary[:, 1])),int(0.6*len(filtered_primary[:, 1]))):
            if filtered_primary[i, 1] < 1.0:
                transit_start = results.folded_phase[i]*results.period
                i_transit_start = i
                break
        for i in range(i_transit_start+1, int(0.6*len(filtered_primary[:, 1]))+4):
            if filtered_primary[i, 1] == 1.0 and filtered_primary[i+1, 1] == 1.0 and filtered_primary[i+2, 1] == 1.0 and filtered_primary[i+3, 1] == 1.0 and filtered_primary[i+4, 1] == 1.0:
                transit_end = results.folded_phase[i]*results.period
                break
        t_dur = (transit_end - transit_start)  # est. of transit duration, T14, in days
        t_depth_norm = np.min(filtered_primary[:, 1])  # est. of transit depth, norm 
        
        # Plot phase folded transit light curve vs time-from-mid-transit
        ax5 = plt.subplot2grid((6, 2), (2, 0), rowspan=1, colspan=1)
        ax5.scatter((results.folded_phase - 0.5)*results.period, (results.folded_y - 1.0)*10**6, color='darkgrey', s=1, alpha=0.5, zorder=2)
        ax5.plot((results.folded_phase - 0.5)*results.period, (filtered_primary[:, 1]-1)*10**6, color='red')
        ax5.set_xlim(-0.15, 0.15)
        ax5.set_xlabel('Time from mid-transit (days)')
        ax5.set_ylabel('Flux (ppm)')
        ax5.text(0.03, 0.30, 'Primary', transform=ax5.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=1.0))
        ax5.text(0.03, 0.15, str(names[j]) + ' planet', transform=ax5.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=1.0))
        
        # Now use lowess to smooth and recover signal in noise for secondary eclipse
        lowess_frac = lowess_frac_sec     # key toggle parameter for lowess filter
        filtered_secondary = lowess(results.folded_y, (results.folded_phase - 0.5)*results.period, lowess_frac)
        
        # if smoothed flux level is greater than flux_limit_sec, then force to 1.0
        for i in range(0, len(filtered_secondary[:, 1])):
            if filtered_secondary[i, 1] > flux_lim_sec:
                filtered_secondary[i, 1] = 1.0  
                  
        # Plot phase folded transit light curve vs time-from-mid-eclipse
        ax51 = plt.subplot2grid((6, 2), (2, 1), rowspan=1, colspan=1)
        ax51.scatter((results.folded_phase - 1.0)*results.period, (results.folded_y - 1.0)*10**6, color='darkgrey', s=1, alpha=0.5, zorder=2)
        ax51.scatter((results.folded_phase)*results.period, (results.folded_y - 1.0)*10**6, color='darkgrey', s=1, alpha=0.5, zorder=2)
        timeL = (results.folded_phase -1.0)*results.period
        fluxL = (filtered_secondary[:, 1] - 1.0)*10**6
        timeR = (results.folded_phase)*results.period
        fluxR = (filtered_secondary[:, 1] - 1.0)*10**6
        time = np.append(timeL, timeR)
        flux = np.append(fluxL, fluxR)
        ax51.set_xlim(-0.15, 0.15)
        ax51.set_xlabel('Time from mid-eclipse (days)')
        ax51.text(0.03, 0.30, 'Secondary', transform=ax51.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=1.0))
        ax51.text(0.03, 0.15, str(names[j]) + ' planet', transform=ax51.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=1.0))
        ax51.set_ylim(-0.06*(1-t_depth_norm)*10**6, 0.015*(1-t_depth_norm)*10**6)
        if (0.05*(1-t_depth_norm)*10**6) < (1-flux_lim_sec)*10**6:
            ax51.text(0.10, 0.55, 'Eclipse depth < instrument precision', transform=ax51.transAxes, color='red', fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=0.8))    
        else:
            ax51.plot(time, flux, 'r--', lw=2)
        
        # Plot phase folded transit light curve vs phase
        ax6 = plt.subplot2grid((6, 2), (3, 0), rowspan=1, colspan=2)
        ax6.plot((results.folded_phase[0:int(0.9*len(results.folded_phase))] - 0.5), (filtered_primary[0:int(0.9*len(results.folded_phase)), 1] - 1.0)*10**6, color='red')
        ax6.scatter((results.folded_phase - 0.5), (results.folded_y - 1.0)*10**6, color='darkgrey', s=1, alpha=0.5, zorder=2)
        ax6.scatter((results.folded_phase + 0.5), (results.folded_y - 1.0)*10**6, color='darkgrey', s=1, alpha=0.5, zorder=2)
        ax6.set_xlim(-0.2, 0.8)
        ax6.set_xlabel('Phase')
        ax6.set_ylabel('Flux (ppm)')
        ax6.text(0.01, 0.15, str(names[j]) + ' planet', transform=ax6.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=1.0))
           
        # Plot the TLS periodogram:
        ax4 = plt.subplot2grid((6, 2), (5, 0), rowspan=1, colspan=2)
        ax4 = plt.gca()
        ax4.axvline(results.period, alpha=0.4, lw=3)
        ax4.set_xlim(np.min(results.periods), np.max(results.periods))
        # Highlight the harmonics of the peak period
        for n in range(2, 10):
            ax4.axvline(n*results.period, alpha=0.4, lw=1, linestyle="dashed")
            ax4.axvline(results.period / n, alpha=0.4, lw=1, linestyle="dashed")
        ax4.set_ylabel(r'SDE')
        ax4.set_xlabel('Period (days)')
        ax4.text(0.60, 0.85,'TLS Periodogram, ' + str(names[j]) + ' planet', transform=ax4.transAxes, fontsize=9)
        ax4.text(0.65, 0.70,'Period (d) = ' + str(format(results.period, '.4f')) + ';  SDE = ' + str(format(results.SDE, '.1f')), transform=ax4.transAxes, fontsize=9)
        ax4.plot(results.periods, results.power, color='black', lw=0.5)
        ax4.set_xlim(0, max(results.periods))
                
        # Plot the full light curve, together with the best transit model (in red)
        ax7 = plt.subplot2grid((6, 2), (4, 0), rowspan=1, colspan=2)
        in_transit = transit_mask(
            t,
            results.period,
            t_dur,
            results.T0)
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
            results.model_lightcurve_time,
            (results.model_lightcurve_model - 1.0)*10**6, alpha=0.5, color='red', zorder=1)
        ax7.set_xlim(min(lc.time), max(lc.time))
        ax7.set_xlabel('Time [BTJD - ' + str(BTJD_base) +']')
        ax7.set_ylabel('Flux (ppm)')
        ax7.text(0.01, 0.15, ' Full light curve', transform=ax7.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=0.8))
        ax7.text(0.81, 0.15, str(names[j]) + ' planet transits', transform=ax7.transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='square', edgecolor='white', facecolor='white', alpha=0.8))
          
        # Print some statistics to the console:
        print(' ')     
        print('Period (d):', format(results.period, '.5f'))
        print('Transit depth (ppm):', format((1-t_depth_norm)*10**6, '.0f'))
        print('Transit duration (hrs):', format((t_dur)*24, '.4f'))
        print('Rp/Rs:', format(np.sqrt(1 - t_depth_norm),'.5f'))
        print('Signal detection eff.(SDE):', format(results.SDE,'.1f'))
        print('SNR:', format(results.snr,'.1f'))
        print(' ') 
          
        # Create subplot with stats only, no axes
        ax3 = plt.subplot2grid((6, 2), (0, 1), rowspan=2, colspan=1)
        ax3.axis('off')
        ax3.text(0.10,0.93, Target +'; TOI-' + TOI + '.xx', transform=ax3.transAxes)
        ax3.text(0.13,0.85, 'Bitmask = ' + str(Bitmask), transform=ax3.transAxes, fontsize=9)
        ax3.text(0.13,0.77, 'Planet candidate = ' + str(pc), transform=ax3.transAxes, fontsize=9)
        ax3.text(0.13,0.69, 'Period (days) = ' + str(format(results.period, '.5f')), transform=ax3.transAxes, fontsize=9)
        ax3.text(0.13,0.61, 'Transit depth (ppm) = ' + str(format((1 - t_depth_norm)*10**6, '.0f')), transform=ax3.transAxes, fontsize=9)
        ax3.text(0.13,0.53, 'Transit duration (hr) = ' + str(format(t_dur*24, '.4f')), transform=ax3.transAxes, fontsize=9)
        ax3.text(0.13,0.45, 'SDE = ' + str(format(results.SDE,'.1f')) + ';  SNR = ' + str(format(results.snr,'.1f')), transform=ax3.transAxes, fontsize=9)
        ax3.text(0.13,0.37, '$R_{*}$ = ' + str(format(R_star, '.2f')) + ' (+' + str(format(R_star_max_err, '.2f')) + ', -' + str(format(R_star_min_err, '.2f')) +') $R_{\odot}$', transform=ax3.transAxes, fontsize=9)
        ax3.text(0.13,0.29, '$M_{*}$ = ' + str(format(M_star, '.2f')) + ' (+' + str(format(M_star_max_err, '.2f')) + ', -' + str(format(M_star_min_err, '.2f')) +') $R_{\odot}$', transform=ax3.transAxes, fontsize=9)
        ax3.text(0.13,0.21, '$R_{p}/R_{*}$ = ' + str(format(np.sqrt(1 - t_depth_norm), '.3f')), transform=ax3.transAxes, fontsize=9)
        ax3.text(0.13,0.13, '$R_{p}$ = ' + str(format(np.sqrt(1 - t_depth_norm)*R_star*109.2, '.2f')) + '$\,R_{\oplus}$' + ' = ' + str(format(np.sqrt(1 - t_depth_norm)*R_star*10.973, '.2f')) + '$\,R_{jup}$', transform=ax3.transAxes, fontsize=9)
        ax3.text(0.13,0.05, 'odd/even mismatch = ' + str(format(results.odd_even_mismatch, '.2f')) + ' $\sigma$', transform=ax3.transAxes, fontsize=9)
        
        # Make certain table variables explicit
        Rp_Rs = np.sqrt(1 - t_depth_norm)
        Rp_j = np.sqrt(1 - t_depth_norm)*R_star*10.973
        odev = results.odd_even_mismatch
        
        plt.tight_layout()   
                
        # Save plots and results to combined vetting sheet for this planet candidate
        fname = [os.getcwd() + '/TLS_Vetting_sheets/' + Target + '_pc_' + str(pc) +'.svg']
        plt.savefig(fname[0], overwrite=True, bbox_inches='tight')
        plt.close('all')
        
        # Write data/results to the summary table for this target/planet
        fout.write(fmt1.format(                       
                 Target    = Target,
                 TOI       = TOI,
                 Mission   = Mission,
                 pc        = pc, 
                 Bitmask   = Bitmask,
                 Period    = results.period,
                 tdepth    = (1-t_depth_norm)*10**6,
                 tdur      = t_dur*24,
                 SDE       = results.SDE,
                 snr       = results.snr,
                 R_star    = R_star,
                 M_star    = M_star,
                 Rp_Rs     = Rp_Rs,
                 Rp_j      = Rp_j,
                 odev      = odev))
        
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