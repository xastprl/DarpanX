'''
An example, to compute the optical funcions for 
an arbitrary ML for incident energy of 8 keV over 
0-5 degree anglular range.


Biswajit, 19-Jul-2020
'''

import darpanx as drp
import numpy as np

#Define some random thickness profile:
expo=np.arange(7)-3
Z_Array=np.exp(-expo)+5.0

#Make multilayer structure
m=drp.Multilayer(MultilayerType="UserDefined",SubstrateMaterial="SiO2",LayerNum=7,LayerMaterial=["W","Si","W","C","Pt","C","Ni"],Z_Array=Z_Array)

Energy=[8.0] # Incident beam energy in KeV
Theta=np.arange(start=0,stop=5,step=0.01) # Define theta range in degree .

m.get_optical_func(Theta=Theta,Energy=Energy,AllOpticalFun ="yes")

m.plot(ylog="no",Comp=["Ra","Ta","Aa"],AllComp="oplot",OutFile="example3_userdefined")
