# Poropyck

## Installation

This will guide you through the steps needed to install ``poropyck``.

### Installing the PyPI package

This ``poropyck`` package is available on PyPI.

    pip install poropyck

Installing via ``pip`` should get most of the Python dependencies you need.

The following is a list of packages used during development:

 * python 3.6.3
 * numpy 1.13.3
 * matplotlib 2.1.0
 * scipy 0.19.1
 * mcerp3 1.0.0

## Scripts

Installing the package will provide you with the following scripts:

### Velocity picking (data on the bench)

 * *simple\_harmonics*: test the DTW picking with simple functions
 * *pick\_vp\_dtw*: do picking of P-waves using Dynamic Time Warping
 * *pick\_vs\_dtw*: do picking of S-waves using Dynamic Time Warping
 * *cross\_corr\_p*: calculate the time lag through cross-correlation for the
   P-wave
 * *cross\_corr\_s*: calculate the time lag through cross-correlation for the
   S-wave

### Velocity picking (data on the vessel)

 * *vp\_dtw\_saturated\_loop*: back-track P-wave picking points between waveforms
   from high to low pressures
 * *vs\_dtw\_saturated\_loop*: back-track S-wave picking points between waveforms
   from high to low pressures

### Densities and porosities

Depending on the denominator used in the Archimedes calculations the results
can be different.  ``densities\_porosity2`` relies on the weights of the
saturated samples in water which are more accurate to measure.  Both
subroutines are included for comparison.

 * *densities\_porosity*: Calculate densities and porosities from the Archimedes
   measurements
 * *densities\_porosity2*: Calculate densities and porosities from the Archimedes
   measurements

### Samples

Sample data can be installed using the following script:

 * *install_poropyck_samples*: installs sample data into the specified location

### Miscellaneous subroutines

Calculate elastic constants and densities from the pycnometer files:

 * *elastic\_constants*
 * *elastic\_constants\_pyc*
 * *pycnometer\_densities*
