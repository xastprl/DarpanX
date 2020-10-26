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
#Purpose:
  This scrpt will download the X-ray form factor, ttenuation and scattering tables data tables provided by NIST (available on: 
  https://www.nist.gov/pml/x-ray-form-factor-attenuation-and-scattering-tables) for the elements Z=1 to Z=92 and each elements 
  scatteing factors are corrected for various factors as mentioned in the Ref: "Detailed Tabulation of Atomic Form Factors, 
  Photoelectric Absorption and Scattering Cross Section, and Mass Attenuation Coefficients in the Vicinity of Absorption Edges 
  in the Soft X-Ray (Z=30–36, Z=60–89, E=0.1 keV–10 keV), Addressing Convergence Issues of Earlier Work, DOI: 10.1063/1.1321055, 
  C. T. Chandler, American Institute of Physics(2001)"

#Last Modified:
  Biswajit, 2nd-Jul-2020
'''

import urllib.request
import numpy as np
import os
import datetime

nk_dir=os.path.join('..','nk_data','nist','nist_sf','test')

def get_html(Z):
    link = 'http://physics.nist.gov/cgi-bin/ffast/ffast.pl?Z=%d&Formula=&gtype=4&lower=&upper=&density=&frames=no' % Z
    http = urllib.request.urlopen(link)
    return http.readlines()

def find_str(file, string):
    for no,line in enumerate(file):
        if line.find(string) != -1:
            return line,no
    return 0,0

for Z in range(1,93):

    file=get_html(Z)

    l,no=find_str(file, b'(Z =')
    l=l.decode('utf-8')
    element=l.split('<b>')[1].split('&#160;')[0].strip()
    print('Downloading data for Element : ',element)

    l,no=find_str(file, b'Form Factors, Attenuation and Scattering Cross-sections')
    l=l.decode('utf-8')
    de=l.split('=')[2].split('keV')[0]
    print('Energy range : ',de,' keV')

    l,no=find_str(file, b'Nominal density')
    l=l.decode('utf-8')
    A=float(l.split(' ')[0])
    
    l,no=find_str(file, b'Nominal density')
    l=l.decode('utf-8')
    Rho=float(l.split('=')[1])
    print('Z, A(g/mol), Rho(g/cm3) = ',Z,', ',A,', ',Rho)

    l,no=find_str(file, b'Relativistic correction')
    l=l.decode('utf-8')
    f_rel=float(l.split('=')[1].split(',')[1].split('<i>e</i>')[0].strip())
    print('Relativistic correction estimate (3/5CL): ',f_rel,' e/atom')

    l,no=find_str(file, b'Nuclear Thomson')
    l=l.decode('utf-8')
    f_NT=float(l.split('=')[1].split('<i>e</i>')[0].strip())
    print('Nuclear Thomson correction facor: ',f_NT,' e/atom')
    print('-------------------------------------------------------------')
    print(' ')
    current_time = datetime.datetime.now()
    f = open(nk_dir+"/"+element+".sf", 'w')
    f.write("#This file is generted by using get_NISTData.py (available within DarpanX package)  by "+str(os.uname()[1])+" on "+str(current_time.year)+"-"+str(current_time.month)+"-"+str(current_time.day)+";\n")
    f.write("#\n")
    f.write("#Element="+str(element)+"\n")
    f.write("#Z="+str(Z)+"\n")
    f.write("#A(g/mol)="+str(A)+"\n")
    f.write("#Rho(g/cm3)="+str(Rho)+"\n")
    f.write("#Energy range="+str(de)+" keV\n")
    f.write("#\n")
    f.write('# E (keV) \t f1 (e/atom) \t f2 (e/atom)\n')

    predata, no = find_str(file, b'Form Factors, Attenuation and Scattering Cross-sections')
    data=file[no+3:-1]
    data_row=[]
    for i in range(len(data)):
        data_row=data_row+[data[i].decode('utf-8')]
    data=np.fromstring(' '.join(data_row), sep = ' ').reshape((len(data_row), 8))
    data_table=np.c_[data[:,0], data[:,1] + f_rel + f_NT, data[:,2]]

    np.savetxt(f, data_table)
    f.close()

