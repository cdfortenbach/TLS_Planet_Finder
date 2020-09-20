# TLS_Planet_Finder
Computer tool to analyze TESS transit light curves for new planet candidates.

## Description
This tool is an implementation of the Transit Least Squares (TLS) algorithm (and python package) described in Hippke & Heller (2019, A&A 623, A39). Here, we have built supporting code around a TLS core to read in a list of (TESS Objects of Interest) TOI's, determine the associated TIC (TESS Input Catalog, host star), and then search for up to four planets around each TIC.  For each planet candidate the code generates a vetting sheet with plots of various transit parameters.  The vetting sheets are based on similar sheets described in Heller et al. (2019, A&A, 627, A66). Upon completion of a run, a summary table is produced that captures all of the data for all of the TOI/TIC's. 

## Getting Started

### Hardware Configuration
* The TLS_Planet_Finder code was developed and tested on a modest desktop computer (see Table 1).  In addition, all production runs to date have been executed on this machine.  

 Table 1. Hardware for Development and Testing

| Component         | Details                |
|:------------------|:-----------------------|
| Processor         | Intel Core i9-9900     |
| Clock             | 3.10 GHz               |
| CPU's             | 8                      |
| Threads/CPU       | 2                      |
| RAM               | 32 GB                  |
| Internal Storage  | 3 TB SSD               |
| External Storage  | 3 TB                   |
| Cooling           | Air                    |


* The full TLS_Planet_Finder code has a very small memory storage footprint of less than 100KB.  A production run could generate a significant memory storage requirement depending on the size of the TIC/TOI data-sets.  The vetting sheet files (.svg) could be as large as 150MB each, but typically are more like 35MB. In execution, again the RAM requirements are dependent on the size of the TIC/TOI data-set that you are trying to analyze.  Normally 4GB of free RAM on top of any other system requirements should be adequate unless you are trying to analyze a very large data-set. 

### Software Configuration/Installation
* The desktop computer described in the previous section runs the Windows 10 Pro operating system.  We chose to develop the TLS_Planet_Finder to run in a Linux environment.  We used Oracle VirtualBox software (6.1.2) to run a Linux guest OS (specifically Ubuntu 18.04) under the Windows 10 host.  Of course a dedicated Linux machine should work just fine.  The installation and execution instructions provided here assume that the user will be using a Linux OS.

* The VirtualBox software allows files to be shared between the two OS's, which is very useful, but otherwise it keeps the two systems separate.  

* The TLS_Planet_Finder code was written in Python 3.7.6.  It is highly recommended that you use Python 3.7.6 to run the TLS_Planet_Finder

* TLS_Planet_Finder has a number of Python dependencies (see Table 2); however, many of them will be satisfied by an up-to-date scientific Python installation.  Unless you actively maintain your own scientific Python distribution, we recommend installing the latest Anaconda Linux distribution.

* It may make sense to use virtualenv (https://pypi.org/project/virtualenv/) to create an isolated Python environment for this project.

* For those packages that are not included in the basic conda distribution, you may need to enable an additional distribution channel with:

 ```conda config --add channels conda-forge```
 
 
 This only needs to be done once, not for each package.

* It is recommended that you install the version/build of the packages specified in Table 2. There are various complex dependencies and this build has been verified.  

* Many of the required packages (unless they have < pip > in the build column of Table 2) can then be installed from the Linux command line with:

 ```conda install package_name=version=build_string```
 
 For example:
 
 ```conda install joblib=0.11=py36_0```
 
* For those packages with a non-conda channel, you should use the command:

 ```conda install -c channel package_name=version=build_string```

 For example:

 ```conda install -c conda-forge astropy=3.2.3```
 

* For the remaining packages, you will need to use the pip package manager.  For example, to install numpy use the command:

 ```pip install numpy==1.18.5``` 

* When this process has been completed for all of the required packages, use the Linux command ```conda list``` to verify that all packages and modules shown in Table 2 are installed:

 Table 2. Python Packages/Modules for TLS_Planet_Finder Installation

| Name                  | Type              | Version/Build                          | Channel                             |
|:----------------------|:------------------|:---------------------------------------|:------------------------------------|
| batman-package        | Python pkg.       | 2.4.6  / pypi_0                        | pypi                                |
| lightkurve            | Python pkg.       | 1.9.0  / py_0                          | conda-forge                         |
| matplotlib            | Python pkg.       | 3.2.2  / py37hef1b27d_0                |                                     |
| numpy                 | Python pkg.       | 1.18.5 / py37ha1c710e_0                | <pip>                               |
| pandas                | Python pkg.       | 1.0.5  / py37h0573a6f_0                |                                     |
| scipy                 | Python pkg.       | 1.5.0  / py37h0b6359f_0                |                                     |
| transitleastsquares   | Python pkg.       | 1.0.25 / pypi_0                        | <pip>                               |


### Installing the TLS_Planet_Finder Code

1. Now, choose/create a working directory for the TLS_Planet_Finder code within your Linux system file structure.

2. Open a terminal in your Linux Download directory.  Then, download the TLS_Planet_Finder repository tar file from GitHub with the command: 

 ```curl -L https://github.com/cdfortenbach/TLS_Planet_Finder/tarball/master > master``` 

   You should see a folder/directory called **master** appear in the Download directory.

3. Now, unpack the **master** tar file with the command: 

 ```tar -xvzf master``` 
 
   you should now see a new folder/directory called **cdfortenbach-TLS_Planet_Finder-*commitID#***.   

4. Open this folder/directory and move the contents to your working directory.  

5. The downloaded tar file (**master**), and the now empty folder/directory **cdfortenbach-TLS_Planet_Finder-*commitID#*** can be deleted.  

6. You should verify that you have all files and folders listed in Table 3.

 Table 3. Files and Folders in the TLS_Planet_Finder Working Directory

| Files/Folders               | Details                                                                                  |
|:----------------------------|:-----------------------------------------------------------------------------------------|
| LICENSE                     | Open-source license GNU GPLv3                                                            |
| README.md                   | README for TLS_Planet_Finder in markdown format                                          |
| TLS_Planet_Finder.23.py     | Python program                                                                           |
| TLS_Smry_tables             | empty folder/directory to collect Summary tables                                         |
| TLS_Vetting_sheets          | empty folder/directory to collect Vetting sheets                                         |
| TLS_config.txt              | file to set certain configuration parameters for the TLS_Planet_Finder                   |
| TOI_list.txt                | list of TOI's that you want to run                                                       |

 Your working directory should now have all of the necessary elements of the TLS_Planet_Finder.  
 
7. That’s it!  You have completed the TLS_Planet_Finder installation.


### Executing the program

The following step-by-step process should guide you through making a run with the TLS_Planet_Finder code:

1.  First, using a GUI (e.g., Ubuntu, etc.) or with a Linux terminal, navigate to the TLS_Planet_Finder working directory.

2.  If you want to start fresh, then the folders **TLS_Smry_tables**, **TLS_Vetting_sheets**, should be cleaned out for a new run.  Any files in these folders/directories can be saved elsewhere as appropriate. The dummy **README** file in each of these folder/directories can be left in place.

3.  Using a .txt editor (or development environment), open the **TLS_config.txt** file and make the appropriate changes to the input parameters.  The following shows an example set-up for a run that looks for two planets for one TIC/TOI.  Only the boldface parameters should be changed.  Only make changes to User Input values in rows 7 - 16, and cols 33 - 40!!   Do not change the syntax of any item.  This is a comma delimited file.  Commas should only be included where indicated.

   TLS_config.txt file:

| Parameter Description                                    | User Input                |
|:---------------------------------------------------------|:--------------------------|
| window_length                                            |, **901**                  |  
| ndist_factor                                             |, **8**                    |    
| outlier_threshold                                        |, **1.65**                 |
| stdev_threshold                                          |, **2.75**                 | 
| sigma_lower                                              |, **10**                   | 
| sigma_upper                                              |, **3**                    | 
| planets_per_system:                                      |, **2**                    | 
| oversampling_factor:                                     |, **3**                    | 
| duration_grid_step:                                      |, **1.1**                  | 
| transit_depth_min:                                       |, **0.00001**              | 
  
   
4. Next, make sure that the **TOI_list.txt** file is properly configured. The comments in the **TOI_list.txt** file provide suggested values for these parameters.

   TOI_list.txt file:

|  TOI     | Mission  |   Bitmask     |                                          |
|:-------------------------------------------------------------------------------|
| #736     |    TESS  |      hard     | # 2 planet system with 1 sector of data  |
| #1260    |    TESS  |      hard     | # 2 planet system with 1 sector of data  |
|  270     |    TESS  |      hard     | # 3 planet system with 3 sectors of data |
   
  * For our example run we will use a file set up for one TOI/TIC.  We will consider TOI 270.  It has three confirmed planets.    In the interest of time and space in this document we will try to demonstrate recovery of two of these planets. Two other TOI's are listed in the **TOI_list.txt** file for example purposes, but they have been commented out.

5. Next, execute the program.  From the Linux command line (from the working directory) enter the following:
 ```
 python3 TLS_Planet_Finder.23.py
 ```

   * After a few seconds, you should see . . . 

```dos
 TLS Planet Finder:
 
 
 Reading in TLS config parameters . . . 
 
 Reading in TOI list . . . 
 
 
 New Target:
 Downloading light curve for the TOI 270 system . . .
 This may take a while . . .
 
 Warning: 31% (6168/19692) of the cadences will be ignored due to the quality mask (quality_bitmask=7407).
 
 SearchResult containing 3 data products.
 
  #   observation  target_name                     productFilename                     distance
 --- ------------- ----------- ------------------------------------------------------- --------
   0 TESS Sector 3   259377017 tess2018263035959-s0003-0000000259377017-0123-s_tp.fits      0.0
   1 TESS Sector 4   259377017 tess2018292075959-s0004-0000000259377017-0124-s_tp.fits      0.0
   2 TESS Sector 5   259377017 tess2018319095959-s0005-0000000259377017-0125-s_tp.fits      0.0
 
 Host star properties: 
   R_star: 0.374  +0.011 -0.011 solar radii
   M_star: 0.362  +0.020 -0.020 solar masses
   Limb darkening model: quadratic;  with coeffs: (0.1604, 0.4325)
 
 Standard cadence: 2.00 (min)
 Total time of dataset incl. gaps: 77.99 (d)
 Fraction of time covered by data: 83.33 (%)
 
 
 ** removed anomalous transient event from stitched/filtered light curve **
  
 Search for the 1st planet using TLS:
   0%|          | 0/8260 periods | 00:00<?Transit Least Squares TLS 1.0.25 (04 June 2020)
 Creating model cache for 45 durations
 Searching 46791 data points, 8260 periods from 0.602 to 38.995 days
 Using all 8 CPU threads
 100%|██████████| 8260/8260 periods | 01:10<00:00
   0%|          | 1/7345 [00:00<18:40,  6.55it/s]Searching for best T0 for period 5.66120 days
 100%|██████████| 7345/7345 [00:05<00:00, 1413.02it/s]
 
 
 Maximum likelihood result (primary transit):
               t0  =  -0.00000 (guess:  0.00000) days
               per =   5.66127 (guess:  5.66120) days
               rp  =   0.05780 (guess:  0.06002) stellar radii
               a   =  25.32472 (guess: 25.45542) stellar radii
               inc =  88.99002 (guess: 88.87451) deg
               ecc =   0.02083 (guess:  0.00000)
               w   =  -0.00000 (guess:  0.00000) deg
               f   =   0.00005 (guess:  0.00001)
         
               t_dur (hrs)     = 1.64
               t_depth (ppm)   = 3737  
 
 
 Maximum likelihood result (secondary eclipse):
               fp  =   0.00000 (guess:  0.00000)
               ts  =   0.00000 (guess:  0.00000) days
         
               ecl_depth (ppm) = 2
 
 Signal detection eff.(SDE): 31.5
 SNR: 56.7
  
  
 Search for the 2nd planet using TLS:
   0%|          | 0/8260 periods | 00:00<?Transit Least Squares TLS 1.0.25 (04 June 2020)
 Creating model cache for 45 durations
 Searching 45531 data points, 8260 periods from 0.602 to 38.995 days
 Using all 8 CPU threads
 100%|██████████| 8260/8260 periods | 01:00<00:00
   1%|▏         | 158/11826 [00:00<00:07, 1578.00it/s]Searching for best T0 for period 11.37609 days
 100%|██████████| 11826/11826 [00:07<00:00, 1687.78it/s]
 

 Maximum likelihood result (primary transit):
               t0  =   0.00000 (guess:  0.00000) days
               per =  11.37609 (guess: 11.37609) days
               rp  =   0.05028 (guess:  0.05226) stellar radii
               a   =  40.17029 (guess: 40.53559) stellar radii
               inc =  89.42016 (guess: 89.29325) deg
               ecc =   0.03809 (guess:  0.00000)
               w   =   0.00000 (guess:  0.00000) deg
               f   =   0.00005 (guess:  0.00001)
         
               t_dur (hrs)     = 2.10
               t_depth (ppm)   = 2841


 Maximum likelihood result (secondary eclipse):
               fp  =   0.00000 (guess:  0.00000)
               ts  =   0.00000 (guess:  0.00000) days
         
               ecl_depth (ppm) = 1
 
 Signal detection eff.(SDE): 30.8
 SNR: 34.2
```

6.  Now we can examine the products of the run.  A summary table should have been deposited into the **TLS_Smry_tables** folder/directory.  It is time and date stamped, and named something like: **TLS_Smry_20200920-094315.txt**.

 * We can also examine the vetting sheets that have been deposited into the **TLS_Vetting_sheets** folder/directory.   For our test run you should see the files listed in Table 4.  
	
 Table 4. Files deposited into the **TLS_Vetting_sheets** folder associated with the test run:

| File name                       | Description                                                          |
|:--------------------------------|:---------------------------------------------------------------------|
| TIC 259377017_pc_1.svg          | full target vetting sheet for TOI 270, planet candidate 1            |
| TIC 259377017_pc_2.svg          | full target vetting sheet for TOI 270, planet candidate 2            |

 
You have completed the example run.  Congratulations!

 
## Author
Charles Fortenbach

## Citations
When publishing results based on usage of the TLS_Planet_Finder please cite:

transitleastsquares:  
Hippke & Heller (2019, A&A 623, A39)  
Heller et al. (2019, A&A, 627, A66)  
https://github.com/hippke/tls

lightkurve:  
Lightkurve Collaboration, (2018)

batman:
Kreidberg, Laura. "batman: BAsic transit model cAlculatioN in Python." Publications of the Astronomical Society of the Pacific 127, no. 957 (2015): 1161.

pandas:  
Data structures for statistical computing in python, McKinney, Proceedings of the 9th Python in Science Conference, Volume 445, 2010.  
10.5281/zenodo.3509134

matplotlib:   
Hunter, J. D. 2007, Computing In Science & Engineering, 9, 90

numpy:   
Oliphant, T. E. 2015, Guide to NumPy, 2nd edn. (USA: CreateSpace Independent Publishing Platform)

scipy: 
Jones, E., Oliphant, T., Peterson, P., et al. 2001, SciPy: Open source scientific tools for Python


## Version History
* 1.0.0 (2020-08-12)
    Initial public release

* 1.1.0 (2020-09-20)
    Extensive update of TLS_Planet_Finder. Replaced light curve smoothing routines with maximum likelihood 'batman' transit model fit for primary and secondary transits. Replaced standard Lightkurve flattening routines with scipy.medfilt and custom outlier removal routines. Added various nan checks and continuation routines to avoid hard crashes in long multi-TOI runs. Added statistics on cadences and overall time covered by data.  Edited various vetting subplots and summary table.

## License
&copy; 2020 Charles Fortenbach

This project is licensed under the GNU GPLv3 License - see the LICENSE.md file for details.

## Acknowledgments
The code discussed here was developed as part of a group project at the University of California, Berkeley, under the guidance of Prof. Courtney Dressing. I would like to thank Andrew Mayo (UCB graduate student) for certain code elements related to downloading and detrending of the TESS lightcurves, and for his many helpful suggestions.  I would also like to thank Shishir Dholakia and Ziyi Lu (UCB students) for their helpful comments.
