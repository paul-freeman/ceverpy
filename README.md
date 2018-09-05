# Poropyck

## Installation

### Option 1: conda

The ``poropyck`` package is available on Anaconda Cloud

    conda install poropyck -c freemapa

This should automatically get the necessary dependencies, when possible, and
is generally the easiest way to go.

### Option 2: pip

The ``poropyck`` package is also available on PyPI.

    pip install poropyck

Installing via ``pip`` should get most of the Python dependencies you need.

### Option 3 (advanced): GitHub

If the other option do no suit your needs, the package source is available on
GitHub.

The main script is located at ``poropyck/pick_dtw.py`` and sample data and
scripts can be found in ``demo/*``.


For reference, the following is a list of packages used during development:

 * python 3.6.3
 * numpy 1.13.3
 * matplotlib 2.1.0
 * scipy 0.19.1

## Execution

As of version 1.4, ``poropyck`` is provided as a library. So you can simply
import it into Python.

The following code sample demonstrates the usage:

    import poropyck

    lengths = [5.256, 5.25, 5.254, 5.254, 5.252, 5.252, 5.258, 5.265, 5.255, 5.252]
    dry_data = 'NM11_2087_4A_dry.csv'
    sat_data = 'NM11_2087_4A_sat.csv'

    dtw = poropyck.DTW(dry_data, sat_data, lengths)
    dry, sat = dtw.pick()

The resulting ``dry`` and ``sat`` variables in this example will contain
dictionary values for *distance*, *time*, and *velocity*, based on the user
picks. These variables will be *ufloat* values from the ``uncertainties``
package and so encapsulate the propogated error. For example, the saturated
velocity mean value could be viewed using ``sat['velocity'].n``. Or the dry
time deviation could be viewed using ``dry['time'].s``.


## Input data

It should be obvious that the ``lengths`` input into the code is just a list
of length measurements. If the list contains only 1 measurement, the length
variable will have no uncertainty, but the code will still work.

The 2 signal files used as input to ``poropyck`` are both CSV files.
**NOTE:** *The data is not expected to begin until line 22 of these files, so
data preprocessing may be necessary to accommodate this.*

For reference, your signal files should follow this format:

    # 21 lines ignored
    ...
    -7.2000e-07,-0.0261719  # line 22
    -7.1600e-07,-0.0267969
    ...
    3.9276e-05,-0.0310156
    # end of file

If you want to test your signal file, you can use the following Python code
to read the signal data (replace ``SIGNAL_FILE`` with your filename):

    import numpy
    print(numpy.loadtxt(SIGNAL_FILE, delimiter=',', skiprows=21).T[:2])
