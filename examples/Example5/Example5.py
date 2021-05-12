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

m.plot(ylog="no",Comp=["Ra","Ta","Aa"],AllComp="oplot",OutFile="example5_Pt_C_DepthGraded",Struc='yes')
