# Poropyck

## Changes in version 1.6

* input files now expect *microsecond* time units, as opposed to *seconds*
* all demo files updated to microsecond time units
* improved handling of calculations with zero error
* number of skipped rows is now defined at the top of the source code
* *time* and *velocity* result plots for template and query now share an x-axis.
* title for summary plot contains better instructions
* added support for single click on summary plot (for *no uncertainty* calculations)
* renamed ``plot_clicks`` to ``plot_time`` (this is more accurate)
* ``±`` no longer prints on results plot when there is no error
* moved *Execution* section of README to come after *Input data* section

## Installation

### Prerequisite - Dynamic Time Warping

You will need to install the [Dynamic Time Warping
package](https://github.com/paul-freeman/dtw) before you will be able to use
``poropyck``. This package is a Cython package and must therefore be compiled
for your local hardware (and NumPy) configuration. Once this is installed,
you will be able to run ``poropyck``.

### Option 1: conda

The ``poropyck`` package is available on Anaconda Cloud

    conda install poropyck -c freemapa

This should automatically get the necessary dependencies, when possible, and
is generally the easiest way to go.

### Option 2: pip

The ``poropyck`` package is also available on PyPI.

    pip install poropyck

Installing via ``pip`` may not get other dependencies, so you should ensure
you have the ``uncertainties`` package installeds, as well as some standards:
``matplotlib``, ``numpy``, and ``scipy``.

### Option 3 (advanced): GitHub

If the other option do no suit your needs, the package source is available on
GitHub.

Installation should be possible using the included ``setup.py`` script:

    python setup.py install

The main script is located at ``poropyck/pick_dtw.py`` and sample data and
scripts can be found in ``demo/*``.

For reference, the following is a list of packages used during development:

 * python 3.6.3
 * numpy 1.13.3
 * matplotlib 2.1.0
 * scipy 0.19.1

This is not to say that these package versions are required, but if you
encounter problems, this may be a place to start your search.

## Input data

It should be obvious that the ``lengths`` input into the code is just a list
of length measurements. If the list contains only 1 measurement, the length
variable will have no uncertainty, but ``poropyck`` should still work.

The signal files used as input to ``poropyck`` should be CSV files.

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

## Execution

As of version 1.4, ``poropyck`` is provided as a library. So you can simply
import it into Python.

Scripts similar to the ones shown in this section, as well as sample data are
provided in the GitHub repository.

### Basic use

The following code sample demonstrates the basic usage:

```python
import poropyck

lengths = [5.256, 5.25, 5.254, 5.254, 5.252, 5.252, 5.258, 5.265, 5.255, 5.252]
template = 'NM11_2087_4A_dry.csv'
query = 'NM11_2087_4A_sat.csv'

dtw = poropyck.DTW(template, query, lengths)
query_results = dtw.pick()
```

The ``query_results`` variable is a dictionary and will contain many values
describing the user picks. It will look something like this:

```json
{
    "file": "NM11_2087_4A_sat.csv",
    "window_start": 14.07,
    "window_end": 15.9,
    "distance": 5.2548,
    "distance_error": 0.004044749683231326,
    "time": 14.985,
    "time_error": 0.4575,
    "velocity": 3506.706706706707,
    "velocity_error": 107.0956363829241,
    "template": {
        "file": "NM11_2087_4A_dry.csv",
        "window_start": 14.836,
        "window_end": 16.64,
        "time": 15.738,
        "time_error": 0.45100000000000007,
        "velocity": 3338.9248951582163,
        "velocity_error": 95.71726030753038
    }
}
```

The template picks are included for reference (and in case the picks need to be repeated).

### Looping example

Here is a longer example showing how you might construct a script to loop over many files:

```python
import json
import matplotlib.pyplot as plt
import poropyck

# read input from JSON
with open('sample2_input.json') as jsonfile:
    data = json.load(jsonfile)

# put files into (template, query) tuples for DTW
templates = data['waves'][:-1]  # last file is never used as template
queries = data['waves'][1:]  # first file is never used as query
picks = zip(templates, queries)

# loop over each (template, query) tuple
data['picks'] = []
for template, query in picks:
    # run poropyck
    dtw = poropyck.DTW(template['file'], query['file'], data['lengths'])
    data['picks'].append(dtw.pick())

# write output to JSON
with open('sample2_output.json', 'w') as jsonfile:
    json.dump(data, jsonfile, indent=2)

# plot velocities with error
plt.errorbar(
    [pick['velocity'] for pick in data['picks']],
    [w['pressure'] for w in data['waves'][1:]],
    xerr=[pick['velocity_error'] for pick in data['picks']]
)
plt.show()
```

The sample JSON file used for input in this example would look like this:

```json
{
"lengths": [
    5.256,
    5.25,
    5.254,
    5.254,
    5.252,
    5.252,
    5.258,
    5.265,
    5.255,
    5.252
],
"wave_type": "P",
"waves": [
    {
    "file": "NM8A-2087-4B_8000_PP_sat500_u1.csv",
    "pressure": 8000
    },
    {
    "file": "NM8A-2087-4B_7000_PP_sat500_u1.csv",
    "pressure": 7000
    },
    {
    "file": "NM8A-2087-4B_6000_PP_sat500_u1.csv",
    "pressure": 6000
    },
    {
    "file": "NM8A-2087-4B_5000_PP_sat500_u1.csv",
    "pressure": 5000
    },
    {
    "file": "NM8A-2087-4B_4000_PP_sat500_u1.csv",
    "pressure": 4000
    },
    {
    "file": "NM8A-2087-4B_3000_PP_sat500_u1.csv",
    "pressure": 3000
    },
    {
    "file": "NM8A-2087-4B_2000_PP_sat500_u1.csv",
    "pressure": 2000
    },
    {
    "file": "NM8A-2087-4B_1000_PP_sat500_u1.csv",
    "pressure": 1000
    }
]
}
```

Obviously your needs may differ from this example, but if you are unfamiliar
with Python programming, this should get you started.
