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

import darpanx as drp
import os,glob

NK_dir=os.path.join(drp.get_dir(),'..','nk_data','nist')

Elements=['SiO2','PtO','PtO2','NiO','CH2','SiC','B4C']
Density=[2.65,14.90,10.2,6.67,0.9,3.21 ,2.52]

n=len(Elements)

for i in range(n):
    m=drp.nkcal(NK_dir=NK_dir,Formula=Elements[i])
    m.get_nk(OutDir=os.path.join(NK_dir,"nk"),Density=Density[i])
    m.plot_sfnk(Comp=['SF','NK'],OutFile=os.path.join(NK_dir,'plots',Elements[i]),ShowPlot=False)

'''
SF_files=sorted(glob.glob(NK_dir+"/sf/*.sf"))
Elements=[]
for i in range(len(SF_files)):
    name=SF_files[i].split('/')[-1]
    Elements=Elements+[name.split('.')[0]]


n=len(Elements)

for i in range(n):
    m=drp.nkcal(NK_dir=NK_dir,Formula=Elements[i])
    m.get_nk(OutDir=os.path.join(NK_dir,"nk"),Density=None)
    m.plot_sfnk(Comp=['SF','NK'],OutFile=os.path.join(NK_dir,'plots',Elements[i]),ShowPlot=False)
'''



