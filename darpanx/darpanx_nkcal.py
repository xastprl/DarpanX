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
Purpose:
    Calculate SF and NK for compound elements.
Biswajit, Jul.04.2020

'''
import datetime
import os,sys
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.constants as sc
from darpanx.darpanx_utilities import*


class nkcal(object):
    def __init__(self,NK_dir=None,Formula=None,Energy=None):
        self.NK_dir=os.path.join(NK_dir,"nk")
        self.SF_dir=os.path.join(NK_dir,"sf")
        self.Formula=Formula
        #self.Density=Density
        self.Energy=Energy


    def get_indiv_element(self,Formula=None):
        element=Formula
        ele=re.findall('[A-Z][^A-Z]*', element)
        no_ele=len(ele)
        elements=[]
        weight=[]
        ch_len=0
        wth_len=0
        for i in range(no_ele):
            try:
                wth1=re.compile("([a-zA-Z]+)([0-9]+)")
                wth = wth1.match(ele[i]).groups()
                ele[i]=wth[0]
                wth = int(wth[1])
                wth_len=wth_len+len(str(wth))
            except: wth = 1
            elements=elements+[ele[i]]
            weight=weight+[wth]
            ch_len=ch_len+len(ele[i])
        if ch_len+wth_len != len(element):raise Exception("%% DarpanX_Error :Incorrect material name. The first letter of each element of a compound material should be in the capital, e.g, for Silicon Dioxide it should be 'SiO2' not 'siO2' etc.")
        return elements,weight

    def find_str(self,File, string):
        for no,line in enumerate(File):
            if line.find(string) != -1:
                return line,no
        return None,None

    def read_file(self,Element=None,File_dir=None):
        SF_dir=File_dir
        if os.path.isfile(os.path.join(SF_dir,Element)) is False:
            raise Exception("%% DarpanX_Error : File not exist- "+os.path.join(SF_dir,Element))
        f=open(os.path.join(SF_dir,Element),'r')
        try:
            l=f.readlines()[0:9]
            f.close()
            A=float(self.find_str(l,'#A(g/mol)')[0].split('=')[1])
            Rho=float(self.find_str(l,'#Rho(g/cm3)=')[0].split('=')[1])
            nist_data=np.loadtxt(SF_dir+'/'+Element)
        except: raise Exception("%% DarpanX_Error : Invalid format of file '"+Element+"' inside NK_dir")
        return A,Rho,nist_data[:,0], nist_data[:,1], nist_data[:,2]

    def get_compound_sf(self,Formula=None,Energy=None):
        SF_dir=self.SF_dir
        ele,wt=self.get_indiv_element(Formula=Formula)
        no_ele=len(ele)
        self.f1=0
        self.f2=0
        total_wt=0
        f1=0
        f2=0
        A=0
        if Energy is None:
            minE=[]
            maxE=[]
            for i in range(no_ele):
                total_wt=total_wt+wt[i]
                globals()['A%s'%str(i)],globals()['Rho%s'%str(i)],globals()['e%s'%str(i)],globals()['f1%s'%str(i)],globals()['f2%s'%str(i)]=self.read_file(Element=ele[i]+'.sf',File_dir=SF_dir)
                minE=minE+[min(globals()['e%s'%str(i)])]
                maxE=maxE+[max(globals()['e%s'%str(i)])]
            # findout common energie range to all data sets for interpolation.
            minE_ind=minE.index(max(minE))
            #maxE_ind=maxE.index(min(maxE))
            minE=max(minE)
            maxE=min(maxE)
            common_e=np.array(globals()['e%s'%str(minE_ind)])
            if common_e[-1] > maxE:
                ind=np.nonzero(common_e <= maxE)
                common_e=common_e[0:max(ind)]
            for i in range(no_ele):
                globals()['F1%s'%str(i)]=interpolate.interp1d(np.log10(globals()['e%s'%str(i)]), (globals()['f1%s'%str(i)]), kind='linear')
                globals()['F2%s'%str(i)]=interpolate.interp1d(np.log10(globals()['e%s'%str(i)]), (globals()['f2%s'%str(i)]), kind='linear')
                globals()['f1%s'%str(i)]=globals()['F1%s'%str(i)](np.log10(common_e))
                globals()['f2%s'%str(i)]=globals()['F2%s'%str(i)](np.log10(common_e))
                A=A+(globals()['A%s'%str(i)])*wt[i]
                f1=f1+(wt[i]*(globals()['f1%s'%str(i)]))
                f2=f2+(wt[i]*(globals()['f2%s'%str(i)]))
                del globals()['F1%s'%str(i)],globals()['A%s'%str(i)],globals()['f1%s'%str(i)]
                del globals()['F2%s'%str(i)],globals()['e%s'%str(i)],globals()['f2%s'%str(i)]
            Energy=common_e
        elif Energy is not None:
            for i in range(no_ele):
                total_wt=total_wt+wt[i]
                globals()['A%s'%str(i)],globals()['Rho%s'%str(i)],globals()['e%s'%str(i)],globals()['f1%s'%str(i)],globals()['f2%s'%str(i)]=self.read_file(Element=ele[i]+'.sf',File_dir=SF_dir)
                globals()['F1%s'%str(i)]=interpolate.interp1d(np.log10(globals()['e%s'%str(i)]), globals()['f1%s'%str(i)], kind='linear')
                globals()['F2%s'%str(i)]=interpolate.interp1d(np.log10(globals()['e%s'%str(i)]), globals()['f2%s'%str(i)], kind='linear')
                globals()['f1%s'%str(i)]=globals()['F1%s'%str(i)](np.log10(Energy))
                globals()['f2%s'%str(i)]=globals()['F2%s'%str(i)](np.log10(Energy))
                A=A+(globals()['A%s'%str(i)])*wt[i]
                f1=f1+(wt[i]*globals()['f1%s'%str(i)])
                f2=f2+(wt[i]*globals()['f2%s'%str(i)])
                del globals()['F1%s'%str(i)],globals()['A%s'%str(i)],globals()['e%s'%str(i)]
                del globals()['F2%s'%str(i)],globals()['f1%s'%str(i)], globals()['f2%s'%str(i)]
        if no_ele > 1: Rho=None
        else: Rho=Rho0
        #f1=f1/total_wt
        #f2=f2/total_wt
        self.f1=f1
        self.f2=f2
        return A,Rho,Energy,f1,f2

    #def get_nk(self,Formula=Formula,Density=None,Energy=None):
    def get_nk(self,OutDir=None,Density=None,FindNK=False):
        '''
        * If Density is given the NK will calcualte for given density
        * If FindNK=True then it will search for existing .nk file of
          the compound. if not exist then calculate from SF data.

        It will return 5 quantity:
        (1) Rho in g/cm3
        (2) Energy array in Angstrom
        (3) db, where db=(Na*r_e*Lambda^2 / 2pi)*(f1-i f2) and RI=n=1-(Na*r_e*Lambda^2 / 2pi)*Rho*(f1-i f2)
        (4) N, Real part of RI
        (5) K, Imaginary part of NK
        '''
        self.A=0
        self.Rho=0
        self.E=0
        self.N=0
        self.K=0
        self.db=0
        NK_dir=self.NK_dir
        SF_dir=self.SF_dir
        Formula=self.Formula
        #Density=self.Density
        Energy=self.Energy
        density_stat=1
        if FindNK is True and os.path.isfile(os.path.join(NK_dir,Formula+'.nk')) is True:
            A,Rho,e_angs,n,k=self.read_file(Element=Formula+'.nk',File_dir=NK_dir)
            self.A=A
            self.Rho=Rho
            self.E=a2kev(e_angs)
            self.N=n
            self.K=k
            db=(1-(n+ 1.0j*k))/Rho
            #self.db=db
            return Rho,e_angs,n,k,db,density_stat
        else:
            #Na= 6.6022e23 #/mol
            Na= 6.6022 # in 1.0e23 /mol
            #r_e=2.82e-13 # cm
            r_e=2.82 #in 1.0e-13 cm
            A,Rho,e,f1,f2=self.get_compound_sf(Formula=Formula,Energy=Energy)
            density_stat=1
            if Density is not None: Rho=Density
            elif Rho is None and Density is None:
                #print("%% DarpanX_Setting : Density is not getting from header of "+Formula+" data file. Set Density=1 (g/cm3). To define density see option < DensityCorrection >. Ignore if < DensityCorrection > = 'yes' is defined")
                print("%% DarpanX_Setting : nkcal is not getting the density from header of "+Formula+" data file. Set Density=1 (g/cm3) in nkcal.")
                Rho=1.0
                density_stat=None
            de=str(min(e))+"-"+str(max(e))
            #cons=(r_e*Rho*Na*1.0e-6)/(2.0*3.14*A)
            cons=(r_e*Na*1.0e-6)/(2.0*3.14*A)
            e_angs=kev2a(e)
            #n=1-(cons*(e_angs**2)*(f1 - 1.0J*f2))
            db=(cons*(e_angs**2)*(f1 - 1.0J*f2))
            n=1-(Rho*db)
            if OutDir != None:
                current_time = datetime.datetime.now()
                f = open(os.path.join(NK_dir,Formula+".nk"), 'w')
                f.write("#This file is generted by using darpanx_nkcal.get_nk() (available within DarpanX package)  by "+str(os.uname()[1])+" on "+str(current_time.year)+"-"+str(current_time.month)+"-"+str(current_time.day)+";\n")
                f.write("#\n")
                f.write("#Formula="+str(Formula)+"\n")
                f.write("#A(g/mol)="+str(A)+"\n")
                f.write("#Rho(g/cm3)="+str(Rho)+"\n")
                f.write("#Energy range="+str(de)+" keV\n")
                f.write("#\n")
                f.write('# Wavelength (A) \t n  \t k \n')
                data_table=np.c_[e_angs, n.real , n.imag]
                np.savetxt(f, data_table)
                f.close()
                self.A=A
                self.Rho=Rho
                self.E=e
                self.N=n.real
                self.K=n.imag
                #self.db=db
                return print("%% DarpanX_message : NK data is saved in- "+os.path.join(NK_dir,Formula+".nk"))
            self.A=A
            self.Rho=Rho
            self.E=e
            #self.db=db
            self.N=n.real
            self.K=n.imag
            return Rho,e_angs,n.real,n.imag,db,density_stat

    def plot_sfnk(self,Comp=None,OutFile=None,ylog='yes',xlog='yes',title=None,OutFileFormat='pdf',xlim=None,ylim=None,ShowPlot=True):
        try:
            n=len(Comp)
        except:
            raise Exception('%% DarpanX_Error : No Componenet is selelted for plot. Select Comp=["'"SF"'", "'"NK"'"]')
        x=self.E
        xlab='Energy (keV)'
        for i in range(n):
            
            fig, ax = plt.subplots(num=None, figsize=(10, 8), dpi=80, facecolor='w', edgecolor='k')
            plt.xticks(size = 20)
            plt.yticks(size = 20)
            if xlim != None:ax.set_xlim(xlim[0],xlim[1])
            if ylim != None:ax.set_ylim(ylim[0],ylim[1])
            if xlog == 'yes':plt.xscale('log')
            if ylog == 'yes':plt.yscale('log')
            if Comp[i] == 'SF':
                if title is None: plt.title("DarpanX Output: Element->"+self.Formula,fontsize=18,color='k',loc='left')
                plt.ylabel('SF (e/atom)', fontsize=20)
                plt.xlabel(xlab, fontsize=20)
                ax.set_xlim(min(self.E),max(self.E))
                plt.grid(True)
                ax.plot(x,self.f1,lw=3.0,label='f1')
                ax.plot(x,self.f2,lw=3.0,label='f2')
                plt.legend(fontsize=24)
                ax.axhline(1.0, color='darkgreen',ls='--', lw=1.5, alpha=0.5)
                if OutFile is not None:
                    plt.savefig(OutFile+'_sf.'+OutFileFormat)
            elif Comp[i] == 'NK':
                if title is None: plt.title("DarpanX Output: Element->"+self.Formula+", Density->"+str(self.Rho)+"$g/cm^{3}$",fontsize=18,color='k',loc='left')
                plt.ylabel('RI', fontsize=20)
                plt.xlabel(xlab, fontsize=20)
                ax.set_xlim(min(self.E),max(self.E))
                plt.grid(True)
                ax.plot(x,self.N,lw=3.0,label='Re(RI)')
                ax.plot(x,self.K,lw=3.0,label='Imag(RI)')
                plt.legend(fontsize=24)
                ax.axhline(1.0, color='darkgreen',ls='--', lw=1.5, alpha=0.5)
                if OutFile is not None:
                    plt.savefig(OutFile+'_nk.'+OutFileFormat)
            else:
                print('%% DarpanX_message : Unknown Comp = '+Comp[i]+' ....Skipping')
        if ShowPlot is True: plt.show(block=False)
