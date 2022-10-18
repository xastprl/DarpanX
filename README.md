# DarpanX: A Python Package for Modeling X-ray Reflectivity of Multilayer Mirrors
<a href="https://ascl.net/2101.015"><img src="https://img.shields.io/badge/ascl-2101.015-blue.svg?colorB=262255" alt="ascl:2101.015" /></a>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7220886.svg)](https://doi.org/10.5281/zenodo.7220886)

DarpanX is a Python-3 package that provides functionality to compute reflectivity and other specular optical functions of multilayer/single-layer mirror for different energy and angles as well as to fit the XRR measurements of the mirrors. An  API  is  provided  for  this and it is also implemented as a local model in PyXspec for  fitting  the experimental XRR measurements. 

For more information regarding the DarpanX see:
* [Mondal, Biswajit et al. 2021, Astronomy and Computing.](https://doi.org/10.1016/j.ascom.2020.100446)


## Pre-requisites

Following python modules are required for the installation of darpanx :

```
 setuptools
 numpy
 scipy
 multiprocessing
 matplotlib
 tabulate
```

For the experimental data fitting, [PyXSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/index.html) (python version of [XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/)) is required.

## Installation
The DarpanX package can be downloaded from the Github link-
https://github.com/biswajitmb/DarpanX.git 

After downloading it, go to the DarpanX directory:

```
cd DarpanX
```

To install DarpanX package system-wide:

```
sudo python3.x setup.py install
```

## Usage

Open a python3.x shell like ipython

```
import darpanx as drp
import numpy as np
```

To create a multilayer object

```
m=drp.Multilayer(MultilayerType="SingleLayer",SubstrateMaterial="SiO2",LayerMaterial=["Pt"],Period=80)
```

Define the angle and energy values for which the optical functions are to be computed:

```
Energy=[8.0] # Incident beam energy in KeV
Theta=np.arange(start=0,stop=5,step=0.01) #Define Grazing angles in degree.
```

Compute reflectivity

```
m.get_optical_func(Theta=Theta,Energy=Energy,AllOpticalFun ="yes")
```

Plot outputs:

```
m.plot(ylog="no", Comp=["Ra","Ta","Aa"], AllComp="oplot", Scale="yes")
```

All available functionality with examples are described in the DarpanX_UserManual.pdf

Use of DarpanX should be acknowledged by referencing the paper:
[Mondal, B. et al. 2021, Astronomy and Computing.](https://doi.org/10.1016/j.ascom.2020.100446)






