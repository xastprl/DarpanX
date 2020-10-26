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

from darpanx.darpanx_frsnl import*
from darpanx.darpanx_utilities import*
import matplotlib.pyplot as plt
import numpy as np
import sys
import multiprocessing as mulp
from darpanx.get_dir import*

class Multilayer(object):

    def __init__(self,NK_dir=None, LayerMaterial=None, Repetition=None, NumStack=None, MultilayerType=None,
            EachRepetitionLayerNum=None, LayerNum=None, Period=None,Gamma=None,D_max=None,D_min=None,C=None,GammaTop=None, SigmaValues=[0.0],
            TopLayerThick=None, EachRepetitionLayerThick=None, Z_Array=None,NumCore=None,AmbientDensity=None,SubstrateDensity=None,
            SubstrateMaterial='Vacuum', TopLayerNum=0,AmbientMaterial='Vacuum',DensityCorrection='no', InterfaceProf='ErrorFunction', 
            SigmaCorrMethod='NevotCroce', BeamPolarization=0,DetectorPA=1, ProjEffCorr='no', SampleSize=100,IncidentBeamDia=0.1, 
            LayerDensity=[], InsRes=0, ShowPar = 'yes'):
     
        #intro()

        try:
            NumCore=int(NumCore)
            if NumCore > 0:
                print("%% DarpanX_status: Parallel processing is using with no.cores = "+str(NumCore))
                self.NumCore=NumCore
            else: print("%% DarpanX_Error: < NumCore > should be an integer, like 1,2,3...")
        except:self.NumCore=None
        
        if NK_dir == None:
            print("%% DarpanX_Settings: Undefined < NK_dir >. Set default direcrory (nk_data/nist)")
            dir1=get_dir()
            NK_dir=os.path.join(str(dir1),"..","nk_data","nist")

        
        self.NK_dir=NK_dir
        self.MultilayerType=MultilayerType
        self.LayerMaterial=LayerMaterial
        self.SubstrateMaterial=SubstrateMaterial
        self.Repetition=Repetition
        self.EachRepetitionLayerNum=EachRepetitionLayerNum
        self.TopLayerNum=TopLayerNum
        self.AmbientMaterial=AmbientMaterial
        self.Period=Period
        self.Gamma=Gamma
        self.D_max=D_max
        self.D_min=D_min
        self.C=C
        self.SigmaValues=SigmaValues
        self.GammaTop=GammaTop
        self.TopLayerThick=TopLayerThick
        self.EachRepetitionLayerThick=EachRepetitionLayerThick
        self.LayerNum=LayerNum
        self.NumStack=NumStack
        self.Z_Array=Z_Array
        self.DensityCorrection=DensityCorrection
        self.InterfaceProf=InterfaceProf
        self.SigmaCorrMethod=SigmaCorrMethod
        self.BeamPolarization=BeamPolarization
        self.DetectorPA=DetectorPA
        self.ProjEffCorr=ProjEffCorr
        self.SampleSize=SampleSize
        self.IncidentBeamDia=IncidentBeamDia
        #self.OriginalDensity=OriginalDensity
        self.Density=LayerDensity
        self.InsRes=InsRes
        #self.NK_Array=NK_Array
        self.AmbientDensity=AmbientDensity
        self.SubstrateDensity=SubstrateDensity
        self.ShowPar = ShowPar

        if self.ShowPar == 'yes':specific(MultilayerType =MultilayerType ,LayerMaterial=LayerMaterial, Repetition=Repetition, NumStack=NumStack,
        EachRepetitionLayerNum=EachRepetitionLayerNum, LayerNum=LayerNum, Period=Period,Gamma=Gamma,D_max=D_max,
        D_min=D_min,C=C,GammaTop=GammaTop, SigmaValues=SigmaValues,TopLayerThick=TopLayerThick, EachRepetitionLayerThick=EachRepetitionLayerThick,
        Z_Array=Z_Array,NumCore=NumCore,AmbientDensity=AmbientDensity,SubstrateDensity=SubstrateDensity,SubstrateMaterial=SubstrateMaterial, 
        TopLayerNum=TopLayerNum,AmbientMaterial=AmbientMaterial,DensityCorrection=DensityCorrection,InterfaceProf=InterfaceProf,SigmaCorrMethod=SigmaCorrMethod, 
        BeamPolarization=BeamPolarization,DetectorPA=DetectorPA, ProjEffCorr=ProjEffCorr,
        SampleSize=SampleSize,IncidentBeamDia=IncidentBeamDia,LayerDensity=LayerDensity, InsRes=InsRes)

    def angular_scan(self,i):
        Rs,Rp,Ra=self.mul.get_reflectivity(self.Z_Array,self.NK_Array,self.sigma,self.Theta_in[i],self.Energy_in,IncidentBeamDia=self.IncidentBeamDia)
        return Rs,Rp,Ra
    def energy_scan(self,i):
        ncls=self.mul.Materials_NK(self.Energy_in[i])
        NK_Array=self.mul.refractiveindexarray(ncls,LayerDensity=self.Density,AmbientDensity=self.AmbientDensity,SubstrateDensity=self.SubstrateDensity)
        Rs,Rp,Ra=self.mul.get_reflectivity(self.Z_Array,NK_Array,self.sigma,self.Theta_in,self.Energy_in[i],IncidentBeamDia=self.IncidentBeamDia)
        return Rs,Rp,Ra

    def angular_scan_all(self,i):
        Rs,Rp,Ra,Ts,Tp,Ta,As,Ap,Aa=self.mul.get_refl_tran_absp(self.Z_Array,self.NK_Array,self.sigma,self.Theta_in[i],self.Energy_in,IncidentBeamDia=self.IncidentBeamDia)
        return Rs,Rp,Ra,Ts,Tp,Ta,As,Ap,Aa
    def energy_scan_all(self,i):
        ncls=self.mul.Materials_NK(self.Energy_in[i])
        NK_Array=self.mul.refractiveindexarray(ncls,LayerDensity=self.Density,AmbientDensity=self.AmbientDensity,SubstrateDensity=self.SubstrateDensity)
        Rs,Rp,Ra,Ts,Tp,Ta,As,Ap,Aa=self.mul.get_refl_tran_absp(self.Z_Array,NK_Array,self.sigma,self.Theta_in,self.Energy_in[i],IncidentBeamDia=self.IncidentBeamDia)
        return Rs,Rp,Ra,Ts,Tp,Ta,As,Ap,Aa

    def get_optical_func(self,Theta=None,Energy=None,AllOpticalFun='no'):
        #global mul,NK_Array,sigma,ncls
        if self.ShowPar == 'yes': indepen_var(Theta=Theta,Energy=Energy)
        Theta=np.array(Theta)
        Energy=np.array(Energy)
        self.Theta=Theta
        self.Energy=Energy
        self.AllOpticalFun=AllOpticalFun

        self.Rs=None
        self.Rp=None
        self.Ra=None
        self.Ts=None
        self.Tp=None
        self.Ta=None
        self.As=None
        self.Ap=None
        self.Aa=None
        del self.Rs,self.Rp,self.Ra,self.Ts,self.Tp,self.Ta,self.As,self.Ap,self.Aa

        if Theta is None or Energy is None:
            raise Exception("%% DarpanX_Error: Define independent varoables, < Theta > or < Energy >")
        if AllOpticalFun != 'yes':
            if len(Energy) == 1 and len(Theta) > 1:
                self.Rs=np.empty(len(Theta))
                self.Rp=np.empty(len(Theta))
                self.Ra=np.empty(len(Theta))

            elif len(Energy) > 1 and len(Theta) == 1:
                self.Rs=np.empty(len(Energy))
                self.Rp=np.empty(len(Energy))
                self.Ra=np.empty(len(Energy))
            else:
                raise Exception("%% DarpanX_Error: Define either < Theta > or < energy > as an array of independent parameters ")

        elif AllOpticalFun == 'yes':
            if len(Energy) == 1 and len(Theta) > 1:
                self.Rs=np.empty(len(Theta))
                self.Rp=np.empty(len(Theta))
                self.Ra=np.empty(len(Theta))
                self.Ts=np.empty(len(Theta))
                self.Tp=np.empty(len(Theta))
                self.Ta=np.empty(len(Theta))
                self.As=np.empty(len(Theta))
                self.Ap=np.empty(len(Theta))
                self.Aa=np.empty(len(Theta))
            elif len(Energy) > 1 and len(Theta) == 1:
                self.Rs=np.empty(len(Energy))
                self.Rp=np.empty(len(Energy))
                self.Ra=np.empty(len(Energy))
                self.Ts=np.empty(len(Energy))
                self.Tp=np.empty(len(Energy))
                self.Ta=np.empty(len(Energy))
                self.As=np.empty(len(Energy))
                self.Ap=np.empty(len(Energy))
                self.Aa=np.empty(len(Energy))
            else:
                raise Exception("%% DarpanX_Error: Define either < Theta > or < energy > as an array of independent parameters ")

        self.mul=ML_Structure(NK_dir=self.NK_dir, LayerMaterial=self.LayerMaterial, Repetition=self.Repetition, 
                    NumStack=self.NumStack, MultilayerType=self.MultilayerType,EachRepetitionLayerNum=self.EachRepetitionLayerNum, 
                    LayerNum=self.LayerNum, SubstrateMaterial=self.SubstrateMaterial, TopLayerNum=self.TopLayerNum,AmbientMaterial=
                    self.AmbientMaterial,DensityCorrection=self.DensityCorrection, InterfaceProf=self.InterfaceProf, SigmaCorrMethod
                    =self.SigmaCorrMethod, BeamPolarization=self.BeamPolarization,DetectorPA=self.DetectorPA, ProjEffCorr=self.ProjEffCorr, 
                    SampleSize=self.SampleSize, InsRes=self.InsRes)

        Theta=deg2rad(Theta)
        Energy=kev2a(Energy)
        self.Theta_in=Theta
        self.Energy_in=Energy
        self.Z_Array=self.mul.LStructure(Period=self.Period,Gamma=self.Gamma,D_max=self.D_max,D_min=self.D_min,C=self.C,GammaTop=self.GammaTop, 
                TopLayerThick=self.TopLayerThick, EachRepetitionLayerThick=self.EachRepetitionLayerThick, Z_Array=self.Z_Array)
        #ncls=mul.Materials_NK(self.Energy)
        #NK_Array=mul.refractiveindexarray(ncls)
        self.sigma=self.mul.roughnessarray(self.SigmaValues)
        #self.sigma=sigma

        if AllOpticalFun != 'yes':
            if len(Energy) == 1 and len(Theta) > 1:
                ncls=self.mul.Materials_NK(Energy)
                self.NK_Array=self.mul.refractiveindexarray(ncls,LayerDensity=self.Density,AmbientDensity=self.AmbientDensity,SubstrateDensity=self.SubstrateDensity)
                if self.NumCore is None:
                    for i in range(len(Theta)): 
                        self.Rs[i],self.Rp[i],self.Ra[i]=self.mul.get_reflectivity(self.Z_Array,self.NK_Array,self.sigma,Theta[i],Energy,IncidentBeamDia=self.IncidentBeamDia)
                else:
                    pool = mulp.Pool(self.NumCore)
                    R = pool.map(self.angular_scan, range(len(self.Theta_in)))
                    pool.close()
                    R=np.array(R)
                    self.Rs=R[:,0]
                    self.Rp=R[:,1]
                    self.Ra=R[:,2]

            elif len(Energy) > 1 and len(Theta) == 1:
                if self.NumCore is None:
                    for i in range(len(Energy)):
                        ncls=self.mul.Materials_NK(Energy[i])
                        self.NK_Array=self.mul.refractiveindexarray(ncls,LayerDensity=self.Density,AmbientDensity=self.AmbientDensity,SubstrateDensity=self.SubstrateDensity)
                        self.Rs[i],self.Rp[i],self.Ra[i]=self.mul.get_reflectivity(self.Z_Array,self.NK_Array,self.sigma,Theta,Energy[i],IncidentBeamDia=self.IncidentBeamDia)
                else:
                    ncls_ForShowPar=self.mul.Materials_NK(Energy[0])
                    NK_Array_ForShowPar=self.mul.refractiveindexarray(ncls_ForShowPar,LayerDensity=self.Density,AmbientDensity=self.AmbientDensity,SubstrateDensity=self.SubstrateDensity)
                    pool = mulp.Pool(self.NumCore)
                    R = pool.map(self.energy_scan, range(len(self.Energy_in)))
                    pool.close()
                    R=np.array(R)
                    self.Rs=R[:,0]
                    self.Rp=R[:,1]
                    self.Ra=R[:,2]

            #self.Theta=rad2deg(self.Theta)
            #self.Energy=a2kev(self.Energy)
            #return self.Rs,self.Rp,self.Ra  

        elif self.AllOpticalFun == 'yes':
            np.seterr(divide='ignore', invalid='ignore')
            if len(Energy) == 1 and len(Theta) > 1:
                ncls=self.mul.Materials_NK(Energy)
                self.NK_Array=self.mul.refractiveindexarray(ncls,LayerDensity=self.Density,AmbientDensity=self.AmbientDensity,SubstrateDensity=self.SubstrateDensity)
                if self.NumCore is None:
                    for i in range(len(Theta)):
                        self.Rs[i],self.Rp[i],self.Ra[i],self.Ts[i],self.Tp[i],self.Ta[i],self.As[i],self.Ap[i],self.Aa[i]=self.mul.get_refl_tran_absp(self.Z_Array,self.NK_Array,self.sigma,Theta[i],Energy,IncidentBeamDia=self.IncidentBeamDia)
                else:
                    pool = mulp.Pool(self.NumCore)
                    R = pool.map(self.angular_scan_all, range(len(self.Theta_in)))
                    pool.close()
                    R=np.array(R)
                    self.Rs=R[:,0]
                    self.Rp=R[:,1]
                    self.Ra=R[:,2]
                    self.Ts=R[:,3]
                    self.Tp=R[:,4]
                    self.Ta=R[:,5]
                    self.As=R[:,6]
                    self.Ap=R[:,7]
                    self.Aa=R[:,8]
                # Show the densities uses in the calculation    

            elif len(Energy) > 1 and len(Theta) == 1:
                if self.NumCore is None:
                    for i in range(len(Energy)):
                        ncls=self.mul.Materials_NK(Energy[i])
                        self.NK_Array=self.mul.refractiveindexarray(ncls,LayerDensity=self.Density,AmbientDensity=self.AmbientDensity,SubstrateDensity=self.SubstrateDensity)
                        self.Rs[i],self.Rp[i],self.Ra[i],self.Ts[i],self.Tp[i],self.Ta[i],self.As[i],self.Ap[i],self.Aa[i]=self.mul.get_refl_tran_absp(self.Z_Array,self.NK_Array,self.sigma,Theta,Energy[i],IncidentBeamDia=self.IncidentBeamDia)
                else:
                    ncls_ForShowPar=self.mul.Materials_NK(Energy[0])
                    NK_Array_ForShowPar=self.mul.refractiveindexarray(ncls_ForShowPar,LayerDensity=self.Density,AmbientDensity=self.AmbientDensity,SubstrateDensity=self.SubstrateDensity)
                    pool = mulp.Pool(self.NumCore)
                    R = pool.map(self.energy_scan_all, range(len(self.Energy_in)))
                    pool.close()
                    R=np.array(R)
                    self.Rs=R[:,0]
                    self.Rp=R[:,1]
                    self.Ra=R[:,2]
                    self.Ts=R[:,3]
                    self.Tp=R[:,4]
                    self.Ta=R[:,5]
                    self.As=R[:,6]
                    self.Ap=R[:,7]
                    self.Aa=R[:,8]
            #self.Theta=rad2deg(self.Theta)
            #self.Energy=a2kev(self.Energy)
            #return self.Rs,self.Rp,self.Ra,self.Ts,self.Tp,self.Ta,self.As,self.Ap,self.Aa
        if self.ShowPar == 'yes': material_prop(LayerDensity=self.mul.LayerDensity,AmbientDensity=self.mul.AmbientDensity,SubstrateDensity=self.mul.SubstrateDensity)

    def plot(self,OpFun='yes',Struc='no',Scale ='no',Comp=None,AllComp=None,OutFile=None,ylog='yes',xlog='no',title='DarpanX Output',OutFileFormat='pdf',xlim=None,ylim=None):
        #np.seterr(divide='ignore', invalid='ignore')
        try:
            n=len(Comp)
        except:
            raise Exception('%% DarpanX_Error : No Componenet is selelted for plot. Select Comp=["'"Ra"'", "'"Rs"'", etc]')
        if len(self.Energy) == 1 and len(self.Theta) > 1:
            x=self.Theta
            xlab='Theta (degree)'
            xsave='Theta'
        elif len(self.Energy) > 1 and len(self.Theta) == 1:
            x=self.Energy
            xlab='Energy (keV)'
            xsave='Energy'
        if OpFun == 'yes':
            if AllComp == 'oplot':
                fig, ax = plt.subplots(num=None, figsize=(10, 8), dpi=80, facecolor='w', edgecolor='k')
                plt.title(title,fontsize=18,color='k',loc='left')
                plt.xticks(size = 20)
                plt.yticks(size = 20)
                if xlim != None:ax.set_xlim(xlim[0],xlim[1])
                if ylim != None:ax.set_ylim(ylim[0],ylim[1])
                if xlog == 'yes':plt.xscale('log')
                if ylog == 'yes':plt.yscale('log')
                plt.ylabel('Optical Functions', fontsize=20)
                plt.xlabel(xlab, fontsize=20)
                #plt.legend(fontsize=24)
            name=''
            func_status=0
            for i in range(n):
                if self.AllOpticalFun != 'yes':
                    if Comp[i] == 'Ts':
                        print('%% DarpanX_message : "'"Ts"'" is not calculated by "'"get_optical_func()"'" for AllOpticalFun="'"no"'"')
                        func_status=1
                    if Comp[i] == 'Tp':
                        print('%% DarpanX_message : "'"Tp"'" is not calculated by "'"get_optical_func()"'" for AllOpticalFun="'"no"'"')
                        func_status=1
                    if Comp[i] == 'Ta':
                        print('%% DarpanX_message : "'"Ta"'" is not calculated by "'"get_optical_func()"'" for AllOpticalFun="'"no"'"')
                        func_status=1
                    if Comp[i] == 'As':
                        print('%% DarpanX_message : "'"As"'" is not calculated by "'"get_optical_func()"'" for AllOpticalFun="'"no"'"')
                        func_status=1
                    if Comp[i] == 'Ap':
                        print('%% DarpanX_message : "'"Ap"'" is not calculated by "'"get_optical_func()"'" for AllOpticalFun="'"no"'"')
                        func_status=1
                    if Comp[i] == 'Aa':
                        print('%% DarpanX_message : "'"Aa"'" is not calculated by "'"get_optical_func()"'" for AllOpticalFun="'"no"'"')
                        func_status=1
                
                if AllComp != 'oplot' and func_status != 1:
                    fig, ax = plt.subplots(num=None, figsize=(10, 8), dpi=80, facecolor='w', edgecolor='k')
                    plt.title(title,fontsize=18,color='k',loc='left')
                    plt.xticks(size = 20)
                    plt.yticks(size = 20)
                    plt.xlabel(xlab, fontsize=20)
                    if xlim != None:ax.set_xlim(xlim[0],xlim[1])
                    if ylim != None:ax.set_ylim(ylim[0],ylim[1])
                    if xlog == 'yes':plt.xscale('log')
                    if ylog == 'yes':plt.yscale('log')
                    #plt.legend(fontsize=24)
                if Comp[i] == 'Rs':
                    #ax.set_xlim(min(self.Theta),max(self.Theta))
                    plt.grid(True)
                    plt.plot(x,self.Rs,label='Reflectivity-s')
                    if AllComp != 'oplot':
                        plt.ylabel('Reflectivity', fontsize=20)
                        plt.legend(fontsize=24)
                    if OutFile is not None and AllComp != 'oplot':
                        plt.savefig(OutFile+'_'+xsave+'_Rs.'+OutFileFormat)
                elif Comp[i] == 'Rp':
                    #ax.set_xlim(min(self.Theta),max(self.Theta))
                    plt.grid(True)
                    plt.plot(x,self.Rp,label='Reflectivity-p')
                    if AllComp != 'oplot':
                        plt.ylabel('Reflectivity', fontsize=20)
                        plt.legend(fontsize=24)
                    if OutFile is not None and AllComp != 'oplot':
                        plt.savefig(OutFile+'_'+xsave+'_Rp.'+OutFileFormat)
                elif Comp[i] == 'Ra':
                    #ax.set_xlim(min(self.Theta),max(self.Theta))
                    plt.grid(True)
                    plt.plot(x,self.Ra,label='Reflectivity-a')
                    if AllComp != 'oplot':
                        plt.ylabel('Reflectivity', fontsize=20)
                        plt.legend(fontsize=24)
                    if OutFile is not None and AllComp != 'oplot':
                        plt.savefig(OutFile+'_'+xsave+'_Ra.'+OutFileFormat)
                elif Comp[i] == 'Ts' and func_status != 1:
                    #ax.set_xlim(min(self.Theta),max(self.Theta))
                    plt.grid(True)
                    plt.plot(x,self.Ts,label='Transmittivity-s')
                    if AllComp != 'oplot':
                        plt.ylabel('Transmittivity', fontsize=20)
                        plt.legend(fontsize=24)
                    if OutFile is not None and AllComp != 'oplot':
                        plt.savefig(OutFile+'_'+xsave+'_Ts.'+OutFileFormat)
                elif Comp[i] == 'Tp' and func_status != 1:
                    #ax.set_xlim(min(self.Theta),max(self.Theta))
                    plt.grid(True)
                    plt.plot(x,self.Tp,label='Transmittivity-p')
                    if AllComp != 'oplot':
                        plt.ylabel('Transmittivity', fontsize=20)
                        plt.legend(fontsize=24)
                    if OutFile is not None and AllComp != 'oplot':
                        plt.savefig(OutFile+'_'+xsave+'_Tp.'+OutFileFormat)
                elif Comp[i] == 'Ta' and func_status != 1:
                    #ax.set_xlim(min(self.Theta),max(self.Theta))
                    plt.grid(True)
                    plt.plot(x,self.Ta,label='Transmittivity-a')
                    if AllComp != 'oplot':
                        plt.ylabel('Transmittivity', fontsize=20)
                        plt.legend(fontsize=24)
                    if OutFile is not None and AllComp != 'oplot':
                        plt.savefig(OutFile+'_'+xsave+'_Ta.'+OutFileFormat)

                elif Comp[i] == 'As' and func_status != 1:
                    #ax.set_xlim(min(self.Theta),max(self.Theta))
                    plt.grid(True)
                    plt.plot(x,self.As,label='Absorptance-s')
                    if AllComp != 'oplot':
                        plt.ylabel('Absorptance', fontsize=20)
                        plt.legend(fontsize=24)
                    if OutFile is not None and AllComp != 'oplot':
                        plt.savefig(OutFile+'_'+xsave+'_As.'+OutFileFormat)
                elif Comp[i] == 'Ap' and func_status != 1:
                    #ax.set_xlim(min(self.Theta),max(self.Theta))
                    plt.grid(True)
                    plt.plot(x,self.Ap,label='Absorptance-p')
                    if AllComp != 'oplot':
                        plt.legend(fontsize=24)
                        plt.ylabel('Absorptance', fontsize=20)
                    if OutFile is not None and AllComp != 'oplot':
                        plt.savefig(OutFile+'_'+xsave+'_Ap.'+OutFileFormat)
                elif Comp[i] == 'Aa' and func_status != 1:
                    #ax.set_xlim(min(self.Theta),max(self.Theta))
                    plt.grid(True)
                    plt.plot(x,self.Aa,label='Absorptance-a')
                    if AllComp != 'oplot':
                        plt.ylabel('Absorptance', fontsize=20)
                        plt.legend(fontsize=24)
                    if OutFile is not None and AllComp != 'oplot':
                        plt.savefig(OutFile+'_'+xsave+'_Aa.'+OutFileFormat)
                else:
                    print('%% DarpanX_message : Unknown Comp = '+Comp[i]+' ....Skipping')
                name=name+Comp[i]
            if AllComp == 'oplot':plt.legend(fontsize=24)
            if AllComp == 'oplot' and OutFile is not None:
                plt.savefig(OutFile+'_'+xsave+'_'+name+'.'+OutFileFormat)


            '''
            fig, ax = plt.subplots(num=None, figsize=(10, 8), dpi=80, facecolor='w', edgecolor='k')
            plt.title(title,fontsize=18,color='k',loc='left')
            plt.xticks(size = 20)
            plt.yticks(size = 20)
            if xlim != None:ax.set_xlim(xlim[0],xlim[1])
            if ylim != None:ax.set_ylim(ylim[0],ylim[1])
            if xlog == 'yes':plt.xscale('log')
            if ylog == 'yes':plt.yscale('log')
            if self.AllOpticalFun != 'yes':
                if len(self.Energy) == 1 and len(self.Theta) > 1:
                    plt.ylabel('Reflectivity', fontsize=20)
                    plt.xlabel('Theta (degree)', fontsize=20)
                    ax.set_xlim(min(self.Theta),max(self.Theta))
                    plt.grid(True)
                    if AllComp != 'oplot':
                        plt.plot(self.Theta,self.Ra,'-r',label='Reflectivity')
                    if AllComp == 'oplot':
                        plt.plot(self.Theta,self.Ra,'-r',label='Average R')
                        plt.plot(self.Theta,self.Rs,'-b',label='S-comp R')
                        plt.plot(self.Theta,self.Rp,'--k',label='P-comp R')
                    plt.legend(fontsize=24)
                    if OutFile is not None:
                        plt.savefig(OutFile+'.pdf')
                elif len(self.Energy) > 1 and len(self.Theta) == 1:
                    plt.ylabel('Reflectivity', fontsize=20)
                    plt.xlabel('Energy (keV)', fontsize=20)
                    ax.set_xlim(min(self.Energy),max(self.Energy))
                    plt.grid(True)
                    plt.plot(self.Energy,self.Ra,'-r',label='Average R')
                    if AllComp == 'oplot':
                        plt.plot(self.Energy,self.Rs,'-b',label='S-comp R')
                        plt.plot(self.Energy,self.Rp,'--k',label='P-comp R')
                    plt.legend(fontsize=24)
                    if OutFile is not None:
                        plt.savefig(OutFile+'.pdf')
            elif self.AllOpticalFun == 'yes':
                if len(self.Energy) == 1 and len(self.Theta) > 1:
                    plt.ylabel('Optical functions', fontsize=20)
                    plt.xlabel('Theta (degree)', fontsize=20)
                    ax.set_xlim(min(self.Theta),max(self.Theta))
                    plt.grid(True)
                    plt.plot(self.Theta,self.Ra,'-r',label="Reflactivity")
                    plt.plot(self.Theta,self.Ta,'--b',label="Transmittivity")
                    plt.plot(self.Theta,self.Aa,'-g',label="Absorptance")
                    plt.legend(fontsize=24)
                    if OutFile is not None:
                        plt.savefig(OutFile+'_Ra.pdf')
                elif len(self.Energy) > 1 and len(self.Theta) == 1:
                    plt.ylabel('Optical functions', fontsize=20)
                    plt.xlabel('Energy (keV)', fontsize=20)
                    ax.set_xlim(min(self.Energy),max(self.Energy))
                    plt.grid(True)
                    plt.plot(self.Energy,self.Ra,'-r',label="Reflactivity")
                    plt.plot(self.Energy,self.Ta,'--b',label="Transmittivity")
                    plt.plot(self.Energy,self.Aa,'-g',label="Absorptance")
                    plt.legend(fontsize=24)
                    if OutFile is not None:
                        plt.savefig(OutFile+'.'+OutFileFormat)
        '''
        def substrate_fill():
            totalthick=(3.0*sum(self.Z_Array)/100.0)#1.8
            plt.fill_between([0,1.1], y1=-totalthick,y2=0.1, color='k',alpha=1.0)
            plt.annotate(self.SubstrateMaterial, (0.50,-totalthick),color='white',fontsize=15.0)
        if Struc == 'yes':
            fig1, ax1 = plt.subplots(num=None, figsize=(6, 8), dpi=80, facecolor='w', edgecolor='k')
            #fig1.tight_layout(pad=0)
            plt.axis('off')
            ax1.set_xlim(0.0,1.1)
            #plt.yticks([])
            cl=['red','blue','darkgreen','yellow','fuchsia','grey','khaki','lime','orange']
            cl1=['green','purple','cyan','aquamarine','chartreuse']
            plt.subplots_adjust(left=0.01,right=0.99, top=0.98,bottom=0.02 )
            plt.title(title,fontsize=18,color='k',position=[0.18,0.975])
            if self.MultilayerType == 'UserDefinedML':
                substrate_fill()
                #totalthick=sum(self.EachRepetitionLayerThick)
                #plt.fill_between([0,1.1], y1=-totalthick,y2=0.1, color='k',alpha=1.0)
                #plt.annotate(self.SubstrateMaterial, (0.46,-totalthick),color='white')
                y1=0.1
                for i in range(self.Repetition):
                    y2=y1
                    for j in range(self.EachRepetitionLayerNum-1,-1,-1):
                        y2=y2+self.EachRepetitionLayerThick[j]
                        if i == 0:plt.fill_between([0,1.1], y1=y1,y2=y2, color=cl[j],alpha=1.0,label=self.LayerMaterial[self.TopLayerNum+j])
                        else: plt.fill_between([0,1.1], y1=y1,y2=y2, color=cl[j],alpha=1.0)
                        y1=y2
                if self.TopLayerNum > 0:
                    for k in range(self.TopLayerNum-1,-1,-1):
                        y2=y2+self.TopLayerThick[k]
                        plt.fill_between([0,1.1], y1=y1,y2=y2, color=cl1[k],alpha=1.0,label=self.LayerMaterial[k])
                        y1=y2

                ax1.legend(loc='lower left', bbox_to_anchor= (0.02, -0.02), ncol=len(self.LayerMaterial),borderaxespad=0, frameon=False,fontsize=18) 

            elif self.MultilayerType == 'BiLayer' or self.MultilayerType == 'DepthGraded' or self.MultilayerType == 'ClusterGraded' or self.MultilayerType == 'SingleLayer':
                substrate_fill()
                #if self.MultilayerType == 'DepthGraded':
                #    totalthick=(self.D_max+self.D_min)/2.0
                #elif self.MultilayerType == 'ClusterGraded':
                #    totalthick=self.Period[0]
                #else: totalthick=self.Period
                #plt.fill_between([0,1.1], y1=-totalthick,y2=0.1, color='k',alpha=1.0)
                #plt.annotate(self.SubstrateMaterial, (0.46,-totalthick),color='white')
                y1=0.1
                for i in range(len(self.Z_Array)-1,-1,-1):
                    if (i % 2.0 == 0):
                        clbil='red'
                        mat=self.LayerMaterial[0]
                    else:
                        clbil='blue'
                        mat=self.LayerMaterial[1]
                    y2=y1+self.Z_Array[i]
                    if i > 1:plt.fill_between([0,1.1], y1=y1,y2=y2, color=clbil,alpha=1.0)
                    elif i == 1: plt.fill_between([0,1.1], y1=y1,y2=y2, color=clbil,alpha=1.0,label=mat)
                    elif i ==0: plt.fill_between([0,1.1], y1=y1,y2=y2, color=clbil,alpha=1.0,label=mat)
                    y1=y2
                ax1.legend(loc='lower left', bbox_to_anchor= (0.02, -0.02), ncol=len(self.LayerMaterial),borderaxespad=0, frameon=False,fontsize=18) 

            elif self.MultilayerType == 'UserDefined':
                substrate_fill()
                y1=0.1
                lm=np.array(self.LayerMaterial)
                dup_mat=find_duplicate(self.LayerMaterial)
                cl_usd=lm.astype('U256')
                
                for i in range(len(dup_mat)):
                    c=np.where(lm == dup_mat[i])[0]
                    cl_usd[c]=cl[i]
                    c35=0
                    for k in range(len(c)):
                        c35=c35+1
                        if c[k] == 0:
                            y1=0.1
                            y2=y1+self.Z_Array[c[k]]
                            #print(c[k],y1,y2)
                        else:
                            y1=0.1+sum(self.Z_Array[0:c[k]])
                            y2=y1+self.Z_Array[c[k]]
                            #print(c[k],y1,y2)
                        #print(y1,y2)
                        if c35 == 1:
                            plt.fill_between([0,1.1], y1=y1,y2=y2, color=cl_usd[c[k]],alpha=1.0,label=self.LayerMaterial[c[k]])
                        else:plt.fill_between([0,1.1], y1=y1,y2=y2, color=cl_usd[c[k]],alpha=1.0)
                
                ax1.legend(loc='lower left', bbox_to_anchor= (0.02, -0.02), ncol=len(dup_mat),borderaxespad=0, frameon=False,fontsize=18)
            if Scale == 'yes':
                TZ=sum(self.Z_Array)
                TZ1="%0.2f"%TZ
                plt.annotate("", xy=(1.0, 0.1), xytext=(1.0, 0.1+TZ),arrowprops=dict(arrowstyle="<->", connectionstyle="arc3",color='white',lw=2.0))
                plt.text(1.0, (0.1+TZ/2.0), str(TZ1)+'$\AA$',{'color': 'black', 'fontsize': 24, 'ha': 'center', 'va': 'center','bbox': dict(boxstyle="round", fc="white", ec="black", pad=0.2)},rotation=90)
                TZ=0
            if OutFile is not None:
                plt.savefig(OutFile+'_structure.'+OutFileFormat) 
        plt.show(block=False)


    '''
    def save(self,OutFile='test_darpanx',SaveAll='yes'):
        current_time = datetime.datetime.now()
        with open(OutFile+".py", "w") as f1:
            f1.write("# This file is generted by DarpanX on: "+str(os.uname()[1])+" at "+str(current_time.year)+"-"+str(current_time.month)+"-"+str(current_time.day)+". Don't edit it:")
            f1.write("NK_dir="+str(self.NK_dir)+'\n')
            f1.write("LayerMaterial="+str(self.LayerMaterial)+'\n')
            f1.write("Repetition="+str(self.Repetition)+'\n')
            f1.write("NumStack="+str(self.NumStack)+'\n')
            f1.write("MultilayerType="+str(self.MultilayerType)+'\n')
            f1.write("EachRepetitionLayerNum="+str(self.EachRepetitionLayerNum)+'\n')
            f1.write("LayerNum="+str(self.LayerNum)+'\n')
            f1.write("Period="+str(self.Period)+'\n')
            f1.write("Gamma="+str(self.Gamma)+'\n')
            f1.write("D_max="+str(self.D_max)+'\n')
            f1.write("D_min="+str(self.D_min)+'\n')
            f1.write("C="+str(self.C)+'\n')
            f1.write("GammaTop="+str(self.GammaTop)+'\n')
            f1.write("SigmaValues="+str(self.SigmaValues)+'\n')
            f1.write("TopLayerThick="+str(self.TopLayerThick)+'\n')
            f1.write("EachRepetitionLayerThick="+str(self.EachRepetitionLayerThick)+'\n')
            f1.write("Z_Array="+str(self.Z_Array)+'\n')
            f1.write("NumCore="+str(self.NumCore)+'\n')
            #f1.write("NK_Array="+str(self.NK_Array)+'\n')
            f1.write("SubstrateMaterial="+str(self.SubstrateMaterial)+'\n')
            f1.write("TopLayerNum="+str(self.TopLayerNum)+'\n')
            f1.write("AmbientMaterial="+str(self.AmbientMaterial)+'\n')
            f1.write("DensityCorrection="+str(self.DensityCorrection)+'\n')
            f1.write("InterfaceProf="+str(self.InterfaceProf)+'\n')
            f1.write("SigmaCorrMethod="+str(self.SigmaCorrMethod)+'\n')
            f1.write("BeamPolarization="+str(self.BeamPolarization)+'\n')
            f1.write("DetectorPA="+str(self.DetectorPA)+'\n')
            f1.write("ProjEffCorr="+str(self.ProjEffCorr)+'\n')
            f1.write("SampleSize="+str(self.SampleSize)+'\n')
            f1.write("IncidentBeamDia="+str(self.IncidentBeamDia)+'\n')
            #f1.write("OriginalDensity="+str(self.OriginalDensity)+'\n')
            f1.write("Density="+str(self.Density)+'\n')
            f1.write("InsRes="+str(self.InsRes)+'\n')
    ''' 


