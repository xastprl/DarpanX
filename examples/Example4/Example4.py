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

m.plot(ylog="no",Comp=["Ra","Ta","Aa"],AllComp="oplot",OutFile="example4_userdefined",Struc='yes')
