'''
An example, to calculate the reflectivity, transmitivity, and absorbtance
of a Pt single layer on SiO2 substrate for incident energy of 8 keV over 
0-5 degree angle range.


Biswajit, 19-Jul-2020
'''
import darpanx as drp 
import numpy as np 

m=drp.Multilayer(MultilayerType="SingleLayer",SubstrateMaterial="SiO2",LayerMaterial=["Pt"],Period=80) # defines the parameters of the multilayer structure.

Energy=[8.0] # Incident beam energy in KeV
Theta=np.arange(start=0,stop=5,step=0.01) # Define Grazing angles in degree .

m.get_optical_func(Theta=Theta,Energy=Energy,AllOpticalFun ="yes") #compute the optical functions.

m.plot(ylog="no", Comp=["Ra","Ta","Aa"], AllComp="oplot",OutFile="Example1_Single_Pt_80", Scale="yes") # Show the outputs
