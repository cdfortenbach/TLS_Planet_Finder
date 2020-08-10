# TLS_Planet_Finder
Computer tool to analyze TESS transit light curves for possible planets.

## Description
This tool is an implementation of the Transit Least Squares (TLS) algorithm (and python package) described in Hippke & Heller (2019, A&A 623, A39). Here, we have built supporting code around a TLS core to read in a list of TOI's, determine the associated TIC (host star), and then search for up to four planets around each TIC.  For each planet candidate the code generates a vetting sheet with plots of various transit parameters.  The vetting sheets are based on similar sheets discussed in Heller et al. (2019, A&A, 627, A66). Upon completion of a run, a summary table is produced that captures all of the data for all of the TOI/TIC's. 

## Getting Started

### Hardware Configuration
* The TLS Planet Finder code was developed and tested on a modest desktop computer (see Table 1).  In addition, all production runs to date have been executed on this machine.  

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


* The full TLS Planet Finder code has a very small memory storage footprint of less than 100KB.  A production run could generate a significant memory storage requirement depending on the size of the TIC/TOI data-sets.  The vetting sheet files (.svg) could be as large as 150MB each, but typically are more like 35MB. In execution, again the RAM requirements are dependent on the size of the TIC/TOI data-set that you are trying to analyze.  Normally 4GB of free RAM on top of any other system requirements should be adequate unless you are trying to analyze a very large data-set. 

### Software Configuration/Installation
* The desktop computer described in the previous section runs the Windows 10 Pro operating system.  We chose to develop the TLS_Planet_Finder to run in a Linux environment.  We used Oracle VirtualBox software (6.1.2) to run a Linux guest OS (specifically Ubuntu 18.04) under the Windows 10 host.  Of course a dedicated Linux machine should work just fine.  The installation and execution instructions provided here assume that the user will be using a Linux OS.

* The VirtualBox software allows files to be shared between the two OS's, which is very useful, but otherwise it keeps the two systems separate.  

* The TLS Planet Finder code was written in Python 3.7.6.  It is highly recommended that you use Python 3.7.6 to run the TLS Planet Finder

* TLS Planet Finder has a number of Python dependencies (see Table 2); however, many of them will be satisfied by an up-to-date scientific Python installation.  Unless you actively maintain your own scientific Python distribution, we recommend installing the latest Anaconda Linux distribution.

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

* When this process has been completed for all of the required packages, again use the Linux command ```conda list``` to verify that all packages and modules shown in Table 2 are installed:

 Table 2. Python Packages/Modules for TLS_Planet_Finder Installation

| Name                  | Type              | Version/Build                          | Channel                             |
|:----------------------|:------------------|:---------------------------------------|:------------------------------------|
| matplotlib            | Python pkg.       | 3.2.2  / py37hef1b27d_0                |                                     |
| numpy                 | Python pkg.       | 1.18.5 / py37ha1c710e_0                | <pip>                               |
| pandas                | Python pkg.       | 1.0.5  / py37h0573a6f_0                |                                     |
| lightkurve            | Python pkg.       | 1.9.0  / py_0                          | conda-forge                         |
| transitleastsquares   | Python pkg.       | 1.0.25 / <pip>                         | pypi                                |
| statsmodels           | Python pkg.       | 0.11.1 / py37h8f50634_2                | conda-forge                         |


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
| TLS_Planet_Finder.2.py      | Python program                                                                           |
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

3.  Using a .txt editor (or development environment), open the **TLS_config.txt** file and make the appropriate changes to the input parameters.  The following shows an example set-up for a run that looks for two planets for each TIC.  Only the boldface parameters should be changed.  Only make changes to User Input values in rows 7 - 18, and cols 33 - 40!!   Do not change the syntax of any item.  This is a comma delimited file.  Commas should only be included where indicated.

   TLS_config.txt file:

| Parameter Description                                    | User Input                |
|:---------------------------------------------------------|:--------------------------|
| window_length                                            |, **501**                  |   
| polyorder                                                |, **3**                    | 
| sigma_lower                                              |, **6**                    | 
| sigma_upper                                              |, **3**                    | 
| planets_per_system:                                      |, **2**                    | 
| lowess_frac_prim:                                        |, **0.010**                | 
| flux_lim_prim:                                           |, **0.999770**             | 
| lowess_frac_sec:                                         |, **0.050**                |
| flux_lim_sec:                                            |, **0.999950**             |
| oversampling_factor:                                     |, **3**                    | 
| duration_grid_step:                                      |, **1.1**                  | 
| transit_depth_min:                                       |, **0.00001**              | 
  
  * The comments in the **TLS_config.txt** file provide suggested values for these parameters.  
4. Next, make sure that the **TOI_list.txt** file is properly configured. The comments in the **TOI_list.txt** file provide suggested values for these parameters.
   
  * For our example run we will use a file set up for one TOI/TIC.  We will consider TOI 736.  It has two confirmed planets that we should be able to recover. 

5. Next, execute the program.  From the Linux command line (from the working directory) enter the following:
 ```
 python3 TLS_Planet_Finder.2.py
 ```

   * After a few seconds, you should see . . . 

```dos
 Reading in TLS config parameters . . . 
 
 Reading in TOI list . . . 
 
 Downloading the light curve for TOI 736 . . .
 This may take a while . . .
 
 
 SearchResult containing 1 data products.
 
 #   observation  target_name                     productFilename                     distance
 --- ------------- ----------- ------------------------------------------------------- --------
   0 TESS Sector 9   181804752 tess2019058134432-s0009-0000000181804752-0139-s_tp.fits      0.0
   0%|          | 0/2174 periods | 00:00<? 
  
 Search for the 1st planet using TLS:
 Transit Least Squares TLS 1.0.25 (04 June 2020)
 Creating model cache for 37 durations
 Searching 15489 data points, 2174 periods from 0.601 to 12.171 days
 Using all 8 CPU threads
 100%|██████████| 2174/2174 periods | 00:08<00:00
   4%|▎         | 482/12907 [00:00<00:02, 4810.39it/s]Searching for best T0 for period 4.98920 days
 100%|██████████| 12907/12907 [00:02<00:00, 4368.46it/s]
  
 Period (d): 4.98920
 Transit depth (ppm): 12071
 Transit duration (hrs): 0.5467
 Rp/Rs: 0.10987
 Signal detection eff.(SDE): 20.2
 SNR: 15.2
  
   0%|          | 0/2174 periods | 00:00<? 
  
 Search for the 2nd planet using TLS:
 Transit Least Squares TLS 1.0.25 (04 June 2020)
 Creating model cache for 37 durations
 Searching 15358 data points, 2174 periods from 0.601 to 12.171 days
 Using all 8 CPU threads
 100%|██████████| 2174/2174 periods | 00:07<00:00
  
 Period (d): 0.94774
 Transit depth (ppm): 5128
 Transit duration (hrs): 0.2063
 Rp/Rs: 0.07161
 Signal detection eff.(SDE): 14.0
 SNR: 7.5
```

6.  Now we can examine the products of the run.  A summary table should have been deposited into the **TLS_Smry_tables** folder/directory.  It is time and date stamped, and named something like: **TLS_Smry_20200810-145210.txt**.

 We can also examine the vetting sheets that have been deposited into the **TLS_Vetting_sheets** folder/directory.   For our test run you should see the files listed in Table 4.  
	
 Table 4. Files deposited into the **TLS_Vetting_sheets** folder associated with the test run:

| File name                       | Description                                                          |
|:--------------------------------|:---------------------------------------------------------------------|
| TIC 181804752_pc_1.svg          | full target vetting sheet for TOI 736, planet candidate 1            |
| TIC 181804752_pc_2.svg          | full target vetting sheet for TOI 736, planet candidate 2            |

 You have completed the example run.  Congratulations!

 

## Author
Charles Fortenbach

## Citations
When publishing results based on usage of the TLS_Planet_Finder please cite:

transitleastsquares:
Hippke & Heller (2019, A&A 623, A39)

lightkurve:
Lightkurve Collaboration, (2018)

statsmodels:
Seabold, Skipper, and Josef Perktold. “statsmodels: Econometric and statistical modeling with python.” Proceedings of the 9th Python in Science Conference. 2010.

pandas:
Data structures for statistical computing in python, McKinney, Proceedings of the 9th Python in Science Conference, Volume 445, 2010.
10.5281/zenodo.3509134

matplotlib: 
Hunter, J. D. 2007, Computing In Science & Engineering, 9, 90

numpy: 
Oliphant, T. E. 2015, Guide to NumPy, 2nd edn. (USA: CreateSpace Independent Publishing Platform)

## Version History
* 1.0.0 (2020-08-10)
    Initial public release

## License
&copy; 2020 Charles Fortenbach

This project is licensed under the GNU GPLv3 License - see the LICENSE.md file for details.

## Acknowledgments
The code discussed here was developed as part of a group project at the University of California, Berkeley, under the guidance of Prof. Courtney Dressing. I would like to thank Andrew Mayo (UCB graduate student) for certain code elements related to detrending of the TESS lightcurves, and for his many helpful suggestions.  I would also like to thank Shishir Dholakia (UCB student) for his helpful comments.
