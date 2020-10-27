#***********************************************************************#
#          DarpanX: A Python Package for Modeling X-ray                 #
#              Reflectivity of Multilayer Mirrors.                      #
#***********************************************************************#
# Copyright (C) 2019-2020. Author: Biswajit Mondal                      #
#                                                                       #
# This program is free software: you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# This program is distributed in the hope that it will be useful,       #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
# GNU General Public License for more details.                          #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with this program. If not, see <http://www.gnu.org/licenses/>.  #
#***********************************************************************#

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
