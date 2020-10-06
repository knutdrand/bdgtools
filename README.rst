==============
Bedgraph Tools
==============


.. image:: https://img.shields.io/pypi/v/bdgtools.svg
        :target: https://pypi.python.org/pypi/bdgtools

.. image:: https://img.shields.io/travis/knutdrand/bdgtools.svg
        :target: https://travis-ci.com/knutdrand/bdgtools

.. image:: https://readthedocs.org/projects/bdgtools/badge/?version=latest
        :target: https://bdgtools.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


.. image:: https://pyup.io/repos/github/knutdrand/bdgtools/shield.svg
     :target: https://pyup.io/repos/github/knutdrand/bdgtools/
     :alt: Updates



Tools to work efficiently with bedgraph files.


* Free software: MIT license
* Documentation: https://bdgtools.readthedocs.io.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install bdgtools.

```bash
pip install bdgtools
```

## Usage
### Command Line
Create signal plots of tss+-500 :
```bash
bdgplot tss CTCF_treat_pileup.bdg genes.bed --rs 1000
```
or a V-plot over entire gene bodies:
```bash
bdgplot v CTCF_treat_pileup.bdg genes.bed
```

Create average signal plot over peak regions
```bash
bdgplot average CTCF_treat_pileup.bdg CTCF_peaks.narrowPeak
```
or a heatmap over all peak regions
```bash
bdgplot average CTCF_treat_pileup.bdg CTCF_peaks.narrowPeak
```

Save figure to a png file
```bash
bdgplot v CTCF_treat_pileup.bdg CTCF_peaks.narrowPeak -o CTCF_vplot.png
```
or save figure data to pickle

```bash
bdgplot v CTCF_treat_pileup.bdg CTCF_peaks.narrowPeak -od CTCF_vplot_data.pkl
```

### Python
Read bedgraph and bedfile from file and show a vplot: 

```python
from bdgtools import read_bedgraph, read_bedfile, VPlot, plot
bedgraph = read_bedgraph("CTCF_treat_pileup.bdg")
regions = read_bedfile("CTCF_peaks.narrowPeak")
max_peak_size = max(np.max(r.ends-r.starts) for r in regions)
plotter = VPlot(figure_width=1000)
df = plotter(bedgraph, regions)
plot(df, plotter, show=True)
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.


Features
--------

* TODO

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
