# Installation

Things necessary to compile the codes succesfully in Python 3.6.3

The python version was install in a linux machine with the anaconda shell script
available at: https://www.anaconda.com/download/ (Last checked: 09/02/2018)

Installation:
bash Anaconda-latest-Linux-x86\_64.sh

## Packages

numpy 1.13.3
matplotlib 2.1.0
scipy 0.19.1
mcerp3 1.0.0
rpy2 2.9.0
R version 3.4.2

Install dtw and proxy packages of R in the right directories.  This requires
copy the packages folders 'dtw' and 'proxy' in the library directory of R.

# Scripts

Installing the package will provide you with the following scripts:

## Velocity picking (data on the bench)

 * simple\_harmonics: test the DTW picking with simple functions (input folder: Dummy)
 * pick\_vp\_dtw: do picking of P-waves using Dynamic Time Warping (input
   folder: NM11\_2087\_4A)
 * pick\_vs\_dtw: do picking of S-waves using Dynamic Time Warping (input
   folder: NM11\_2087\_4A)
 * cross\_corr\_p: calculate the time lag through cross-correlation for the
   P-wave (input folder: NM11\_2087\_4A)
 * cross\_corr\_s: calculate the time lag through cross-correlation for the
   S-wave (input folder: NM11\_2087\_4A)

## Velocity picking (data on the vessel)

 * vp\_dtw\_saturated\_loop: back-track P-wave picking points between waveforms
   from high to low pressures. (input folder: NM11-2087-4B\_sat500) 
 * vs\_dtw\_saturated\_loop: back-track S-wave picking points between waveforms
   from high to low pressures. (input folder: NM11-2087-4B\_sat500) 

## Densities and porosities

Depending on the denominator used in the Archimedes calculations the results
can be different.  Densities\_porosity2.py relies on the weights of the
saturated samples in water which are more accurate to measure.  Both
subroutines are included for comparison.

 * densities\_porosity: Calculate densities and porosities from the Archimedes
   measurements (input folder: NM11\_2087\_4A)
 * densities\_porosity2: Calculate densities and porosities from the Archimedes
   measurements (input folder: NM11\_2087\_4A)

## Miscellaneous subroutines

Calculate elastic constants and densities from the pycnometer files:

 * elastic\_constants
 * elastic\_constants\_pyc
 * pycnometer\_densities 
