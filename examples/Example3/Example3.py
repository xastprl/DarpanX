'''
An example, to calculate the reflectivity of a Pt/SiC ML 
of 40 repetition with an extra PtO2 layer on the top formed 
on a SiO2 substrate for incident energy of 8 keV over 0-10 
degree anglular range.


Biswajit, 19-Jul-2020
'''

import darpanx as drp
import numpy as np

m=drp.Multilayer(MultilayerType="UserDefinedML",SubstrateMaterial="SiO2",TopLayerNum=1,EachRepetitionLayerNum=2,LayerMaterial=["PtO2","Pt","SiC"],TopLayerThick=[10],EachRepetitionLayerThick=[20,40],Repetition=40,SigmaValues=[1.0])

Energy=[8.0] # Incident beam energy in KeV
Theta=np.arange(start=0,stop=10,step=0.01) # Define theta range in degree

m.get_optical_func(Theta=Theta,Energy=Energy)

m.plot(ylog="yes",Comp=["Ra"],OutFile="Example2_Pt_SiC")
