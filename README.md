# Dusty Galaxy Analyzer
Dusty Galaxy Analyzer is a tool to extract dusty galaxies from `SIMBA` simulations and analyze them in Python

## Get

You can download it via
```bash
git clone https://bitbucket.org/lq3552/dust_gal_analyzer
```

## Build

DGalA requires:

 * python(2.7)
 * numpy(>=1.12.0)
 * matplotlib(>=2.0.0)
 * [yt](http://yt-project.org)(>=3.4.1)
 * [caesar](http://caesar.readthedocs.io/en/latest/index.html)

You can install it by setup.py:
```bash
python setup.py install
```

or, if you need root access,
```bash
sudo python setup.py install
```

## Basic usage

To extract dusty galaxies from GIZMO snapshots:
```python
from dust_gal_analyzer import DustyGalaxyExtractor 
dge = DustyGalaxyExtractor(fname,replace=0)
dge.savec()
glist = dge.gal_extract()
```

To analyze dusty galaxy properties or statistics, such as the dust-mass function:
```python
from dust_gal_analyzer import DustyGalaxy
gal =  DustyGalaxy(glist)
gal.plot_dmf()
```

## Tools

You can find a tool to predict dust-to-gas ratio in `machine-learning-tools`, using a regressor trained by simulated data. 
Please define which input parameters are active following `param.list`, and set values of input parameters following `tab.txt`.
Then you can run the tool by
```bash
python dgr_extrarandomtree.py param.list tab.txt [DGR_output_list]
```
