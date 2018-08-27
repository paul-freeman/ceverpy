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

The main script is located at ``poropyck/pick_dtw.py`` and sample data can be
found in ``poropyck/demo/*``.


For reference, the following is a list of packages used during development:

 * python 3.6.3
 * numpy 1.13.3
 * matplotlib 2.1.0
 * scipy 0.19.1
 * mcerp3 1.0.0

## Execution

After installation, you should be able to run ``poropyck`` from the command-line.

    poropyck

This will execute using sample data. To use your own data, you must specify a
*template signal*, *query signal*, and *metadata*.

    poropyck -t TEMPLATE_SIGNAL -q QUERY_SIGNAL -m METADATA

## Data files

The 3 files used as input to ``poropyck`` are all simple text files. The
signal files are standard CSV files. **NOTE:** *Signal data is not expected
to begin until line 22 of these files, so data preprocessing may be necessary
to accommodate this.*

Your signal files should look like this:

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

The metadata should contain JSON length data at key location
``['length']['raw']``. The value at this key location should be a list of one
or more decimal length measurements. These values are used to compute an
uncertain length measurement. Any other JSON data is ignored.

Your metadata file should look like this:

    {
        "length": {
            "raw": [
                5.256,
                5.250,
                5.254,
                5.254,
                5.252,
                5.252,
                5.258,
                5.265,
                5.255,
                5.252
            ]
        }
    }

To test your metadata file, you can use this Python code (replace
METADATA_FILE with your filename):

    import json
    with open(METADATA_FILE) as dat:
        metadata = json.load(dat)
    print(metadata['length']['raw'])
