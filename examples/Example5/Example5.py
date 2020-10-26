'''
An example, to compute the optical funcions for 
a Depth-Graded ML mirror for incident energy of 
8 keV over 0-5 degree anglular range.


Biswajit, 19-Jul-2020
'''

import darpanx as drp
import numpy as np

m=drp.Multilayer(MultilayerType="DepthGraded",SubstrateMaterial="SiO2",LayerMaterial=["Pt","C"],Repetition=150,D_max=128.1,D_min=31.7,C=0.245,GammaTop=0.7,Gamma=0.45,SigmaValues=[4.5])


Theta=[0.08] # Incident grazing angle in degree.
Energy=10**(np.arange(start=np.log10(0.1),stop=np.log10(80.0),step=0.01)) # 0.1-80.0 keV Energy in linear grid.

m.get_optical_func(Theta=Theta,Energy=Energy,AllOpticalFun ="yes")

m.plot(ylog="no",Comp=["Ra","Ta","Aa"],AllComp="oplot",OutFile="example3_Pt_C_DepthGraded")
