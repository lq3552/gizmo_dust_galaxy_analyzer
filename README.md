# DGalA
DGalA (dusty galaxy analyzer) is a tool to extract dusty galaxies from GIZMO simulations and analyze them in Python

## Get DGalA

Currently this package is private.

Internal members can get it via
```bash
git clone https://bitbucket.org/lq3552/dust_gal_analyzer
```

## Build DGalA

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
## Usage

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
