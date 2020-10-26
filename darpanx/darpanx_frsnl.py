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

import sys
import os
import numpy as np
from scipy import interpolate
from darpanx.darpanx_utilities import*
import scipy.constants as const
import darpanx.darpanx_nkcal as classnkcal 

class ML_Structure(object):

    def __init__(self, NK_dir=None, LayerMaterial=None, Repetition=None, NumStack=None, MultilayerType=None,
            EachRepetitionLayerNum=None, LayerNum=None, SubstrateMaterial='Vacuum', TopLayerNum=0,AmbientMaterial='Vacuum', 
            DensityCorrection='no', InterfaceProf='ErrorFunction', SigmaCorrMethod='NevotCroce', BeamPolarization=0, 
            DetectorPA=1, ProjEffCorr='no', SampleSize=100, InsRes=0):
        
        self.NK_dir = NK_dir
        self.MultilayerType = MultilayerType
        self.Repetition = Repetition
        self.LayerMaterial = LayerMaterial
        self.SubstrateMaterial = SubstrateMaterial
        self.AmbientMaterial = AmbientMaterial
        self.BeamPolarization = BeamPolarization
        self.DetectorPA = DetectorPA
        self.ProjEffCorr = ProjEffCorr 
        self.DensityCorrection = DensityCorrection
        #self.OriginalDensity = OriginalDensity
        #self.Density = Density
        self.InsRes = InsRes
        self.SampleSize=SampleSize

        dirstatus = os.path.isdir(NK_dir)
        if NK_dir == None:
            raise Exception("%% DarpanX_Error: Undefined < NK_dir >")
        else:
            dirstatus = os.path.isdir(NK_dir)
            while dirstatus != True:
                print("%% DarpanX_Input: NK_dir is not exist.")
                print('Give the < NK_dir > path:')
                NK_dir = str(input())
                dirstatus = os.path.isdir(NK_dir)
                


        if ProjEffCorr == 'yes':
            try:
                self.SampleSize = float(SampleSize)
            except:raise Exception('%% DarpanX_Error: Incorrect < SampleSize >')

        #if ScanAxis is None:
        #    raise Exception('%% DarpanX_Error: Undefined < ScanXAxis > --Options are; < "'"theta"'" >, < "'"ene"'" >')
        try:
            if LayerMaterial is None or len(LayerMaterial) < 1:
                raise Exception('%% DarpanX_Error: Undefined < LayerMaterial > -- It should be a string array, like ["'"Si"'"] or ["'"Pt"'","'"C"'"]')
            else: pass
        except:raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial > -- It should be a string array, like ["'"Si"'"] or ["'"Pt"'","'"C"'"]')
        if self.NK_dir is None:raise Exception('%% DarpanX_Error: Undefined < NK_dir >')

        if InterfaceProf == 'ErrorFunction':
            self.InterfaceProf = 0
        elif InterfaceProf == 'ExponentialFunction':
            self.InterfaceProf = 1
        elif InterfaceProf == 'LinearFunction':
            self.InterfaceProf = 2
        elif InterfaceProf == 'SinusoidalFunction':
            self.InterfaceProf = 3
        elif InterfaceProf == 'StepFunction':
            self.InterfaceProf = 4
        else:
            raise Exception('%% DarpanX_Error: Undefined < InterfaceProf > --options are; < "'"ErrorFunction"'" >, < "'"ExponentialFunction"'" >, < "'"LinearFunction"'" >, < "'"SinusoidalFunction"'" >, < "'"StepFunction"'" > .')
        if SigmaCorrMethod == 'NevotCroce': self.SigmaCorrMethod = 1
        elif SigmaCorrMethod == 'DebyeWaller': self.SigmaCorrMethod = 0
        elif SigmaCorrMethod == 'Avg': self.SigmaCorrMethod = 2
        else :
            raise Exception('%% DarpanX_Error: Undefined < SigmaCorrMethod > --options are; < "'"NevotCroce"'" >, < "'"DebyeWaller"'" >')


        density_stat=0
        if SubstrateMaterial != 'Vacuum':
            nk_class=classnkcal.nkcal(NK_dir=NK_dir,Formula=SubstrateMaterial,Energy=None)
            Rho,SubLambda, SubN, SubK, DB,density_statout=nk_class.get_nk(FindNK=True)
            self.SubstrateDensity=Rho
            self.density_stat_sub=density_statout
            #SubLambda, SubN, SubK = np.loadtxt(self.NK_dir + SubstrateMaterial+'.nk', comments=';', unpack=True)    
            #self.SubN_Array = interpolate.interp1d(SubLambda, SubN, kind='linear')
            #self.SubK_Array = interpolate.interp1d(SubLambda, SubK, kind='linear')
            self.SubN_Array = interpolate.interp1d(SubLambda, DB.real, kind='linear')
            self.SubK_Array = interpolate.interp1d(SubLambda, DB.imag, kind='linear')
            del SubLambda,SubN,SubK,nk_class,Rho,DB
        elif SubstrateMaterial == 'Vacuum':
            self.SubstrateDensity=1.0
            self.density_stat_sub=1
        if AmbientMaterial != 'Vacuum':
            nk_class=classnkcal.nkcal(NK_dir=NK_dir,Formula=AmbientMaterial,Energy=None)
            Rho,AmbLambda, AmbN, AmbK, DB,density_statout =nk_class.get_nk(FindNK=True)
            self.AmbientDensity=Rho
            self.density_stat_amb=density_statout
            #AmbLambda, AmbN, AmbK = np.loadtxt(self.NK_dir + AmbientMaterial+'.nk', comments=';', unpack=True)
            self.AmbN_Array = interpolate.interp1d(AmbLambda, DB.real, kind='linear')
            self.AmbK_Array = interpolate.interp1d(AmbLambda, DB.imag, kind='linear')
            del AmbLambda,AmbN,AmbK ,nk_class,Rho,DB
        elif AmbientMaterial == 'Vacuum': 
            self.AmbientDensity=1.0
            self.density_stat_amb=1
        if MultilayerType == 'UserDefinedML':
            try:
                self.TopLayerNum = int(TopLayerNum)
                self.EachRepetitionLayerNum = EachRepetitionLayerNum
                self.Repetition = int(Repetition)
                self.n = int((EachRepetitionLayerNum*Repetition) + TopLayerNum)
                self.Z_Array = np.empty(self.n)
                #if self.Repetition == 1:
                self.OriginalDensity=np.empty(len(self.LayerMaterial))
                self.Sigma_Array = np.zeros(self.n+1)
                #else:
                #    self.Sigma_Array = np.zeros(self.n+2)
                self.NK_Array = np.empty((self.n+2), dtype=np.complex_)
                #self.interface=np.zeros(self.n+1)+self.InterfaceProf
                self.interface = self.InterfaceProf
            except: raise Exception("%% DarpanX_Error: Incorrect or Undefined  < TopLayerNum >, < EachRepetitionLayerNum >, < Repetition > for UserDefinedML structure")
            if (self.TopLayerNum + self.EachRepetitionLayerNum) != len(self.LayerMaterial):
                raise Exception('%% DarpanX_Error: Inconsistent < TopLayerNum > and < EachRepetitionLayerNum > with the length(< LayerMaterial >)')
            if self.Repetition < 1 :
                raise Exception('%% DarpanX_Error: Incorrect < Repetition >. It should be an integer number, like, 1,2,3...etc')
            density_stat=[]
            for topnum in range(self.TopLayerNum+self.EachRepetitionLayerNum):
                nk_class=classnkcal.nkcal(NK_dir=NK_dir,Formula=LayerMaterial[topnum],Energy=None)
                Rho,LnLambda,LnN,LnK,DB,density_statout=nk_class.get_nk(FindNK=True) 
                self.OriginalDensity[topnum]=Rho
                density_stat=density_stat+[density_statout]
                #LnLambda,LnN,LnK = np.loadtxt(self.NK_dir+LayerMaterial[topnum]+'.nk',comments=';',unpack=True)
                globals()['N_Array%s' % str(topnum)] = interpolate.interp1d(LnLambda, DB.real, kind='linear')
                globals()['K_Array%s' % str(topnum)] = interpolate.interp1d(LnLambda, DB.imag, kind='linear')
                del nk_class,Rho,LnLambda,LnN,LnK,DB

        elif MultilayerType == 'SingleLayer': 
            try:
                if len(LayerMaterial) != 1:
                    raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial > -- should be a 1-element string array like ["'"Pt"'"]')
            except: raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial > -- should be a 1-element string array like ["'"Pt"'"]')
            self.n =1
            self.Z_Array = np.empty(self.n)
            self.Sigma_Array = np.zeros(self.n+1)
            self.NK_Array = np.empty((self.n+2), dtype=np.complex_)
            self.interface = self.InterfaceProf
            nk_class=classnkcal.nkcal(NK_dir=NK_dir,Formula=LayerMaterial[0],Energy=None)
            Rho,L1Lambda,L1N,L1K,DB1,density_statout=nk_class.get_nk(FindNK=True)
            #L1Lambda,L1N,L1K=np.loadtxt(self.NK_dir+self.LayerMaterial[0]+'.nk',comments=';',unpack=True)
            self.L1N_Array = interpolate.interp1d(L1Lambda, DB1.real, kind='linear')
            self.L1K_Array = interpolate.interp1d(L1Lambda, DB1.imag, kind='linear')
            self.OriginalDensity=[Rho]
            density_stat=[density_statout]
            del L1Lambda,L1N,L1K,DB1,Rho,nk_class

        elif MultilayerType == 'BiLayer':
            try:
                if len(LayerMaterial) != 2:
                    raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial > -- should be a 2-element array, like ["'"W"'","'"Si"'"].')
            except:raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial > -- should be a 2-element array, like ["'"W"'","'"Si"'"].')
            try:
                self.n = 2*Repetition
                self.Z_Array = np.empty(self.n)
                self.Sigma_Array = np.zeros(self.n+1)
                self.NK_Array = np.empty((self.n+2), dtype=np.complex_)
                #self.interface=np.zeros(self.n+1)+self.InterfaceProf
                self.OriginalDensity=[0,0]
                self.interface = self.InterfaceProf
            except: raise Exception("%% DarpanX_Error: < Repetition > is required for Bilayer Structure.")
            density_stat=[0,0]
            nk_class=classnkcal.nkcal(NK_dir=NK_dir,Formula=LayerMaterial[0],Energy=None)
            Rho,L1Lambda,L1N,L1K,DB1,density_statout=nk_class.get_nk(FindNK=True)
            self.OriginalDensity[0]=Rho
            density_stat[0]=density_statout
            nk_class=classnkcal.nkcal(NK_dir=NK_dir,Formula=LayerMaterial[1],Energy=None)
            Rho,L2Lambda,L2N,L2K,DB2,density_statout=nk_class.get_nk(FindNK=True)
            self.OriginalDensity[1]=Rho
            density_stat[1]=density_statout
            #L1Lambda,L1N,L1K=np.loadtxt(self.NK_dir+LayerMaterial[0]+'.nk',comments=';',unpack=True)
            #L2Lambda,L2N,L2K =np.loadtxt(self.NK_dir+LayerMaterial[1]+'.nk',comments=';',unpack=True)
            self.L1N_Array = interpolate.interp1d(L1Lambda, DB1.real, kind='linear')
            self.L1K_Array = interpolate.interp1d(L1Lambda, DB1.imag, kind='linear')
            self.L2N_Array = interpolate.interp1d(L2Lambda, DB2.real, kind='linear')
            self.L2K_Array = interpolate.interp1d(L2Lambda, DB2.imag, kind='linear')
            del L1Lambda,L1N,L1K,L2Lambda,L2N,L2K,Rho,nk_class,DB1,DB2


        elif MultilayerType == 'DepthGraded':
            if len(self.LayerMaterial) != 2:
                raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial > -- should be a 2-element array, like ["'"W"'","'"Si"'"].')
            try:
                self.n = 2*self.Repetition
                self.Z_Array = np.empty(self.n)
                self.Sigma_Array = np.zeros(self.n+1)
                self.NK_Array = np.empty((self.n+2), dtype=np.complex_)
                #self.interface=np.zeros(self.n+1)+self.InterfaceProf
                self.OriginalDensity=[0,0]
                self.interface = self.InterfaceProf
            except: raise Exception("%% DarpanX_Error: < Repetition > is required for DepthGraded Structure.")
            density_stat=[0,0]
            nk_class=classnkcal.nkcal(NK_dir=NK_dir,Formula=LayerMaterial[0],Energy=None)
            Rho,L1Lambda,L1N,L1K,DB1,density_statout=nk_class.get_nk(FindNK=True)
            self.OriginalDensity[0]=Rho
            density_stat[0]=density_statout
            nk_class=classnkcal.nkcal(NK_dir=NK_dir,Formula=LayerMaterial[1],Energy=None)
            Rho,L2Lambda,L2N,L2K,DB2,density_statout=nk_class.get_nk(FindNK=True)
            self.OriginalDensity[1]=Rho
            density_stat[1]=density_statout
            #L1Lambda,L1N,L1K=np.loadtxt(self.NK_dir+LayerMaterial[0]+'.nk',comments=';',unpack=True)
            #L2Lambda,L2N,L2K =np.loadtxt(self.NK_dir+LayerMaterial[1]+'.nk',comments=';',unpack=True)
            self.L1N_Array = interpolate.interp1d(L1Lambda, DB1.real, kind='linear')
            self.L1K_Array = interpolate.interp1d(L1Lambda, DB1.imag, kind='linear')
            self.L2N_Array = interpolate.interp1d(L2Lambda, DB2.real, kind='linear')
            self.L2K_Array = interpolate.interp1d(L2Lambda, DB2.imag, kind='linear')
            del L1Lambda,L1N,L1K,L2Lambda,L2N,L2K,Rho,nk_class,DB1,DB2

        elif MultilayerType == 'ClusterGraded':
            if len(self.LayerMaterial) != 2:
                raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial > -- should be a 2-element array, like ["'"W"'","'"Si"'"].')
            try:
                self.NumStack = int(NumStack)
                self.n = 2*self.Repetition
                self.Z_Array = np.empty(self.n)
                self.Sigma_Array = np.zeros(self.n+1)
                self.NK_Array = np.empty((self.n+2), dtype=np.complex_)
                #self.interface=np.zeros(self.n+1)+self.InterfaceProf
                self.OriginalDensity=[0,0]
                self.interface = self.InterfaceProf
            except: raise Exception("%% DarpanX_Error: < Repetition >, < NumStack > are required for ClusterGraded Structure.")
            density_stat=[0,0]
            nk_class=classnkcal.nkcal(NK_dir=NK_dir,Formula=LayerMaterial[0],Energy=None)
            Rho,L1Lambda,L1N,L1K,DB1,density_statout=nk_class.get_nk(FindNK=True)
            self.OriginalDensity[0]=Rho
            density_stat[0]=density_statout
            nk_class=classnkcal.nkcal(NK_dir=NK_dir,Formula=LayerMaterial[1],Energy=None)
            Rho,L2Lambda,L2N,L2K,DB2,density_statout=nk_class.get_nk(FindNK=True)
            self.OriginalDensity[1]=Rho
            density_stat[1]=density_statout
            #L1Lambda,L1N,L1K=np.loadtxt(self.NK_dir+LayerMaterial[0]+'.nk',comments=';',unpack=True)
            #L2Lambda,L2N,L2K =np.loadtxt(self.NK_dir+LayerMaterial[1]+'.nk',comments=';',unpack=True)
            self.L1N_Array = interpolate.interp1d(L1Lambda, DB1.real, kind='linear')
            self.L1K_Array = interpolate.interp1d(L1Lambda, DB1.imag, kind='linear')
            self.L2N_Array = interpolate.interp1d(L2Lambda, DB2.real, kind='linear')
            self.L2K_Array = interpolate.interp1d(L2Lambda, DB2.imag, kind='linear')
            del L1Lambda,L1N,L1K,L2Lambda,L2N,L2K,Rho,nk_class,DB1,DB2

        elif MultilayerType == 'UserDefined':
            try:
                self.LayerNum=int(LayerNum)
                self.Z_Array = np.empty(LayerNum)#np.array(self.LayerNum).astype(float)
                self.Sigma_Array = np.zeros(LayerNum+1)
                self.NK_Array = np.empty((LayerNum+2), dtype=np.complex_)
                #self.interface=np.zeros(self.n+1)+self.InterfaceProf
                self.OriginalDensity=np.empty(len(self.LayerMaterial))
                self.interface = self.InterfaceProf
                self.n = LayerNum
            except: raise Exception("%% DarpanX_Error: Undefined or incorrect  < LayerNum >.")
            if len(self.LayerMaterial) != LayerNum:
                raise Exception('%% DarpanX_Error: Incorrect dimension of  < LayerMaterial > -- len(LayerMaterial) should be a '+str(LayerNum))
            #if DensityCorrection == 'yes' and len(OriginalDensity) != LayerNum:
            #    raise Exception('%% DarpanX_Error: Incorrect dimension of < OriginalDensity >..It should be='+str(len(LayerMaterial)))
            density_stat=[]
            for topnum in range(self.LayerNum):
                nk_class=classnkcal.nkcal(NK_dir=NK_dir,Formula=LayerMaterial[topnum],Energy=None)
                Rho,LnLambda,LnN,LnK,DB,density_statout=nk_class.get_nk(FindNK=True)
                self.OriginalDensity[topnum]=Rho
                density_stat=density_stat+[density_statout]
                #LnLambda,LnN,LnK = np.loadtxt(self.NK_dir+LayerMaterial[topnum]+'.nk',comments=';',unpack=True)
                globals()['N_Array%s' % str(topnum)] = interpolate.interp1d(LnLambda, DB.real, kind='linear')
                globals()['K_Array%s' % str(topnum)] = interpolate.interp1d(LnLambda, DB.imag, kind='linear')
                del nk_class,Rho,LnLambda,LnN,LnK,DB

        else:
            raise Exception('%% DarpanX_Error: Undefined < MultilayerType > --options are < ["'"UserDefinedML"'","'"SingleLayer"'", "'"DepthGraded"'", "'"ClusterGraded"'", "'"UserDefined"'"] >.')
        self.density_stat=density_stat
    #--------------------- LStructure() --------------------
    #
    # Purpose: It will define the Multi-Layer Structure from
    #          input parameters.
    #          
    # Modification History:
    #                     June.27.2019,... Biswajit
    #                     May.27.2020,.... Biswajit
    #-------------------------------------------------------


    def LStructure(self,Period=None,Gamma=None,D_max=None,D_min=None,C=None, 
            GammaTop=None, TopLayerThick=None, EachRepetitionLayerThick=None, Z_Array=None):

        if EachRepetitionLayerThick is not None:
            try:
                err12=len(EachRepetitionLayerThick)
            except: raise Exception("%% DarpanX_Error: Incorrect < EachRepetitionLayerThick >. It should be a float array.")
        if TopLayerThick is not None:
            try:
                err12=len(TopLayerThick)
            except: raise Exception("%% DarpanX_Error: Incorrect < TopLayerThick >. It should be a float array.")
        if self.MultilayerType == 'UserDefinedML':
            if EachRepetitionLayerThick is None:
                raise Exception('%% DarpanX_Error: Undefined < EachRepetitionLayerThick >')
            elif len(EachRepetitionLayerThick) != self.EachRepetitionLayerNum :
                raise Exception('%% DarpanX_Error: len(EachRepetitionLayerThick){'+str(len(EachRepetitionLayerThick))+'} should be = < EachRepetitionLayerNum >{'+str(self.EachRepetitionLayerNum)+'}')
            else: 
                for i in range(self.Repetition):
                    self.Z_Array[(self.TopLayerNum+(self.EachRepetitionLayerNum*i)):(self.TopLayerNum+(self.EachRepetitionLayerNum*(i+1)))] = EachRepetitionLayerThick[:]
                
            if self.TopLayerNum > 0 and TopLayerThick is None:
                raise Exception('%% DarpanX_Error: Undefined < TopLayerThick >')
            elif self.TopLayerNum > 0 and len(TopLayerThick) != self.TopLayerNum:
                raise Exception('%% DarpanX_Error: len(TopLayerThick){'+str(len(TopLayerThick))+'} should be = < TopLayerNum >{'+str(self.TopLayerNum)+'}')
            elif self.TopLayerNum > 0: 
                self.Z_Array[0:self.TopLayerNum] = TopLayerThick[:]
            return self.Z_Array

        elif self.MultilayerType == 'SingleLayer':
            try:
                self.Z_Array[0] = Period
            except: raise Exception("%% DarpanX_Error: Incorrect or Undefined < Period > ... it is required for Single-layer Structure and should be a float.")
            return self.Z_Array
            
        elif self.MultilayerType == 'BiLayer':
            try:
                Z_L1 = Gamma*Period
                Z_L2 = (1.0-Gamma)*Period
            except: raise Exception("%% DarpanX_Error: Incorrect or Undefined < Period >, < Gamma > ... are required for Bilayer Structure.")    
            self.Z_Array[0::2]=Z_L1
            self.Z_Array[1::2]=Z_L2
            return self.Z_Array

        elif self.MultilayerType == 'DepthGraded':
            try:
                k_g = (D_min/D_max)**(1.0/C)
                b_g = (1.0-((self.Repetition)*k_g))/(k_g-1.0)
                a_g = D_min*(b_g+(self.Repetition))**C
                d1 = a_g/(b_g + np.arange(self.Repetition) + 1)**(C)
            except: raise Exception("%% DarpanX_Error: Incorrect or Undefined < D_max >,< D_min >,< C >, < Gamma > all are required for DepthGraded Structure.")
            if GammaTop is None:
                print("%% DarpanX_Setting:  GammaTop is not provided. Set GammaTop = Gamma ")
                GammaTop = Gamma
            self.Z_Array[0::2]=Gamma*d1
            self.Z_Array[1::2]=(1.0-Gamma)*d1
            self.Z_Array[0]=GammaTop*d1[0]
            self.Z_Array[1]=(1-GammaTop)*d1[0]
            return self.Z_Array

        elif self.MultilayerType == 'ClusterGraded':
            NumStack=self.NumStack
            LayerInEachStuck = (self.n)//self.NumStack
            try:
                for i in range(1,NumStack+1):
                    Z_L1=Gamma[i-1]*Period[i-1]
                    Z_L2=(1.0-Gamma[i-1])*Period[i-1]
                    self.Z_Array[((LayerInEachStuck)*(i-1)):((LayerInEachStuck)*i):2]=Z_L1
                    self.Z_Array[(((LayerInEachStuck)*(i-1))+1):((LayerInEachStuck)*i):2]=Z_L2
            except: raise Exception("%% DarpanX_Error: Undefined < Period >, < Gamma > or Incorrect dimension. len(Period) and len(Gamma) should be- "+str(NumStack))
            return self.Z_Array

        elif self.MultilayerType == 'UserDefined':
            try:
                self.Z_Array[:] = Z_Array
            except: raise Exception("%% DarpanX_Error: Undefined < Z_Array > or incorrect dimension. len(Z_Array) should be- "+str(self.LayerNum))
            return self.Z_Array
    
    
    #------------------ Materials_NK() ------------------
    #Modification History:
    #                    June.27.2019,... Biswajit
    #                    May.31.2020,.... Biswajit
    #    
    #---------------------------------------------------

    def Materials_NK(self,Energy):
       
        if self.SubstrateMaterial == 'Vacuum':
            ncs = complex(0.0, 0.0)
        else:
            RI_rs = self.SubN_Array(Energy)
            RI_is = self.SubK_Array(Energy)
            ncs = complex(RI_rs, RI_is)  # ;RI of Substrate
            del RI_rs,RI_is

        if self.AmbientMaterial == 'Vacuum':
            nca = complex(0.0, 0.0)
        else:
            RI_ra = self.AmbN_Array(Energy)
            RI_ia = self.AmbK_Array(Energy)
            nca = complex(RI_ra, RI_ia)  # ;RI of Ambient
            del RI_ra,RI_ia

        if self.MultilayerType =='UserDefinedML' :
            ncls = [nca]
            for topnum in range(self.TopLayerNum+self.EachRepetitionLayerNum):
                RI_rn = globals()['N_Array%s' % str(topnum)](Energy)
                RI_in = globals()['K_Array%s' % str(topnum)](Energy)
                ncls = ncls+[complex(RI_rn,RI_in)]
            ncls = ncls+[ncs]
            del RI_rn,RI_in
            return ncls

        if self.MultilayerType == 'SingleLayer':
            RI_r1 = self.L1N_Array(Energy)
            RI_i1 = self.L1K_Array(Energy)
            ncls = [nca,complex(RI_r1, RI_i1),ncs]
            del RI_r1,RI_i1
            return ncls

        elif self.MultilayerType =='BiLayer' or self.MultilayerType == 'DepthGraded' or self.MultilayerType == 'ClusterGraded':
            RI_r1 = self.L1N_Array(Energy)
            RI_i1 = self.L1K_Array(Energy)
            RI_r2 = self.L2N_Array(Energy)
            RI_i2 = self.L2K_Array(Energy)
            # ;RI of two layer at WL=lambda and substrate
            ncls = [nca,complex(RI_r1, RI_i1), complex(RI_r2, RI_i2),ncs]
            del RI_r1,RI_i1,RI_r2,RI_i2
            return ncls


        elif self.MultilayerType =='UserDefined':
            ncls = [nca]
            for topnum in range(self.LayerNum):
                RI_rn = globals()['N_Array%s' % str(topnum)](Energy)
                RI_in = globals()['K_Array%s' % str(topnum)](Energy)
                ncls = ncls+[complex(RI_rn,RI_in)]
            ncls = ncls+[ncs]
            del RI_rn,RI_in
            return ncls

    #--------------- refractiveindexarray() ----------------
    #
    # Purpose: It will give the RI Structure.
    # Modification History:
    #                     June.27.2019,... Biswajit
    #                     May.27.2020,.... Biswajit
    #-------------------------------------------------------
    def refractiveindexarray(self,ncls,LayerDensity=[],AmbientDensity=None,SubstrateDensity=None):
        if AmbientDensity is None :
            AmbientDensity=self.AmbientDensity
            if self.density_stat_amb is None: raise Exception('%% DarpanX_Error: Density of the ambient medium of '+self.AmbientMaterial+' is not getting from the header of the NK database or the NK file not exist. Give it by the keyword <AmbientDensity>')
        else: self.AmbientDensity=AmbientDensity
        if SubstrateDensity is None : 
            SubstrateDensity=self.SubstrateDensity
            if self.density_stat_sub is None: raise Exception('%% DarpanX_Error: Density of the '+self.SubstrateMaterial+' substrate is not getting from the header of the NK database or the NK file not exist. Give it by the keyword <SubstrateDensity>')
        else: self.SubstrateDensity=SubstrateDensity
        if self.DensityCorrection != 'yes': 
            LayerDensity = self.OriginalDensity
            if (None in self.density_stat) is True: raise Exception('%% DarpanX_Error: Density of the layer material/materials are not getting from the header of the NK database or the NK file not exist. Give it by the keyword < DensityCorrection > and < LayerDensity >')
        try:
            if self.DensityCorrection =='yes' and len(LayerDensity) != len(self.LayerMaterial):
                raise Exception('%% DarpanX_Error: Incorrect dimension of < LayerDensity >..It should be='+str(len(self.LayerMaterial)))
        except:raise Exception('%% DarpanX_Error: Incorrect dimension or value of < LayerDensity > -- It should be a float array, like [5.0] or [16.0,10.0]')
        #else: raise Exception("%% DarpanX_Error: Incorrect option for < DensityCorrection >. Available options are 'yes' or 'no' ")
        self.LayerDensity=LayerDensity
        n=self.n
        #ncls[0]=1.0-(ncls[0]*AmbientDensity)
        #ncls[-1]=1.0-(ncls[-1]*SubstrateDensity)
        ncls0=1.0-(ncls[0]*AmbientDensity)
        nclsn=1.0-(ncls[-1]*SubstrateDensity)
        if self.MultilayerType == 'SingleLayer':
            b=ncls[1]*(LayerDensity[0])
            b=1.0-b
            self.NK_Array[0]=ncls0#ncls[0]
            self.NK_Array[1]=b
            self.NK_Array[2]=nclsn#ncls[-1]
            del ncls
        elif self.MultilayerType =='BiLayer' or self.MultilayerType == 'DepthGraded' or self.MultilayerType == 'ClusterGraded':
            b=ncls[1]*(LayerDensity[0])
            b=1.0-b

            bb=ncls[2]*(LayerDensity[1])
            bb=1.0-bb
            
            self.NK_Array[0]=ncls0#ncls[0]
            self.NK_Array[n+1]=nclsn#ncls[-1]
            ncl=[b,bb]
            self.NK_Array[1:-1:2]=ncl[0]
            self.NK_Array[2:-1:2]=ncl[1]
            del ncls,ncl
        elif self.MultilayerType =='UserDefinedML':
            #nlayersno=self.TopLayerNum+self.EachRepetitionLayerNum
            #if self.DensityCorrection ==  'yes' and len(self.OriginalDensity) == nlayersno and len(Density) == nlayersno :
            #if self.DensityCorrection ==  'yes':
            ncl=[]
            for i in range(1,len(ncls)-1):
                b=ncls[i]*(LayerDensity[i-1])
                b=1.0-b
                ncl=ncl+[b]
                
            toplno=self.TopLayerNum+1
            self.NK_Array[0]=ncls0#ncls[0]
            self.NK_Array[n+1]=nclsn#ncls[-1]
            self.NK_Array[1:toplno]=ncl[0:(self.TopLayerNum)]
            for i in range(self.Repetition):
                self.NK_Array[(toplno+(self.EachRepetitionLayerNum*i)):(toplno+(self.EachRepetitionLayerNum*(i+1)))]=ncl[self.TopLayerNum::]

        elif self.MultilayerType =='UserDefined':
            ncl=[]
            for i in range(1,len(ncls)-1):
                b=ncls[i]*(LayerDensity[i-1])
                b=1.0-b
                ncl=ncl+[b]

            toplno=self.LayerNum+1
            self.NK_Array[0]=ncls0#ncls[0]
            self.NK_Array[n+1]=nclsn#ncls[-1]
            self.NK_Array[1:toplno]=ncl[0:(self.LayerNum)]

        return self.NK_Array
    
    #----------------- roughnessarray() ---------------------
    #
    # Purpose: It will give the Sigma array.
    # Modification History:
    #                     June.27.2019,... Biswajit
    #                     May.27.2020,.... Biswajit
    #-------------------------------------------------------
    def roughnessarray(self,SigmaValues):
        n=self.n
        SigmaValues=np.array(SigmaValues)
        try:
            err13=len(SigmaValues)
        except:raise Exception('%% DarpanX_Error: Incorrect < SigmaValues > ... It should be a float array.')
        if self.MultilayerType =='UserDefinedML':
            if len(SigmaValues) == 1:
                self.Sigma_Array[:]=SigmaValues[0]
            elif self.Repetition == 1 and len(SigmaValues) == (self.EachRepetitionLayerNum+self.TopLayerNum)+1 :
                self.Sigma_Array[0:self.TopLayerNum]=SigmaValues[0:self.TopLayerNum]
                self.Sigma_Array[self.TopLayerNum:(self.TopLayerNum+self.EachRepetitionLayerNum+1)]=SigmaValues[self.TopLayerNum::]
            elif self.Repetition > 1 and len(SigmaValues) == self.EachRepetitionLayerNum+self.TopLayerNum+2 :
                self.Sigma_Array[0:(self.TopLayerNum+1)]=SigmaValues[0:(self.TopLayerNum+1)]
                for i in range(self.Repetition):
                    self.Sigma_Array[(self.TopLayerNum+1+(self.EachRepetitionLayerNum*i)):(self.TopLayerNum+1+(self.EachRepetitionLayerNum*(i+1)))]=SigmaValues[self.TopLayerNum+1:-1]
                self.Sigma_Array[-1]=SigmaValues[-1]
            else:
                if self.Repetition == 1:
                    raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions.len(SigmaValues) should be- 1 or '+str(self.EachRepetitionLayerNum+self.TopLayerNum+1))
                else:
                    raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions. len(SigmaValues) should be- 1 or '+str(self.EachRepetitionLayerNum+self.TopLayerNum+2)+', or '+str(self.n+1))

        elif self.MultilayerType == 'SingleLayer':
            if len(SigmaValues) == 1:
                self.Sigma_Array[:]=SigmaValues[0]
            elif len(SigmaValues) == 2 :
                self.Sigma_Array[0]=SigmaValues[0]
                self.Sigma_Array[1]=SigmaValues[1]
            else:
                raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions. len(SigmaValues) should be- 1 or 2')

        elif self.MultilayerType =='BiLayer' or self.MultilayerType == 'DepthGraded' or self.MultilayerType == 'ClusterGraded':
            if n > 2:
                if len(SigmaValues) == 1:
                    self.Sigma_Array[:]=SigmaValues[0]
                elif (len(SigmaValues) == 4):
                    self.Sigma_Array[0]=SigmaValues[0]
                    #self.Sigma_Array[n]=SigmaValues[3]
                    #----------------------------------------------------------------
                    self.Sigma_Array[1:n:2] = SigmaValues[1]      #All even index elements starting from 1 to n-1 of sigma is replaced by sigma_l1_l2 
                    self.Sigma_Array[2:n:2] = SigmaValues[2]      #All odd index element of sigma is replaced by sigma_l1_l2
                    #---------------------------------------------------------------
                    self.Sigma_Array[n]=SigmaValues[3]
                elif (self.MultilayerType == 'ClusterGraded' and len(SigmaValues) == (self.NumStack*2+2)):
                    #---------------For the rughness correcrtion at individual stack in case of typ=4 --------------
                    layer_each_stuck=(n)//self.NumStack
                    for i in range(1,self.NumStack+1):
                        self.Sigma_Array[((layer_each_stuck)*(i-1))+1:((layer_each_stuck)*i):2]=SigmaValues[(2*i)-1]
                        self.Sigma_Array[(((layer_each_stuck)*(i-1))+2):((layer_each_stuck)*i+1):2]=SigmaValues[2*i]
                    #----------------------------------------------------------------------------------------------
                    self.Sigma_Array[0]=SigmaValues[0]
                    self.Sigma_Array[n]=SigmaValues[-1]

                elif (len(SigmaValues) == n+1):
                    self.Sigma_Array[:]=SigmaValues[:]
                else:
                    if (self.MultilayerType == 'ClusterGraded' and self.Repetition ==1): raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions len(SigmaValues) should be- 1 or 3 or '+str(self.n+1)+' or '+str((2*self.NumStack)+2))
                    elif (self.MultilayerType == 'ClusterGraded' and self.Repetition > 1): raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions len(SigmaValues) should be- 1 or 4 or '+str(self.n+1)+' or '+str((2*self.NumStack)+2))
                    elif (self.MultilayerType != 'ClusterGraded' and self.Repetition ==1) : raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions len(SigmaValues) should be- 1 or 3 or '+str(self.n+1))
                    else: raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions. len(SigmaValues) should be- 1 or 4 or '+str(self.n+1) )
            elif n==2:
                if len(SigmaValues) == 1:
                    self.Sigma_Array[:]=SigmaValues[0]
                elif (len(SigmaValues) == 3):
                    self.Sigma_Array[:]=SigmaValues[:]
                else:
                    raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions. len(SigmaValues) should be- 1 or 3.')
        
        elif self.MultilayerType =='UserDefined':
            if len(SigmaValues) == 1:
                self.Sigma_Array[:]=SigmaValues[0]
            elif len(SigmaValues) == self.n+1:
                self.Sigma_Array=SigmaValues
            else:raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions. len(SigmaValues) should be- 1 or '+str(self.n+1))
        return self.Sigma_Array

    #------- get_reflectivity(),get_refl_tran_absp() ------
    #
    # Purpose: Returns the reflectivity of the whole system
    #          of n number of layers for a single theta/ene
    #          rgy (theta_i/energy_i) value.
    # Modification History:
    #                     June.27.2019,... Biswajit
    #                     May.31.2020,.... Biswajit
    #-------------------------------------------------------
    def get_reflectivity(self,Z_Array,NK_Array,sigma,theta_i,energy_i,IncidentBeamDia=0.1):
        n=self.n
        theta0 = theta_i
        LAMBDA = energy_i
        
        rijs = rij_s(NK_Array[0:-1], NK_Array[1::], theta0, NK_Array[0])
        rijp = rij_p(NK_Array[0:-1], NK_Array[1::], theta0, NK_Array[0])
        if  self.SigmaCorrMethod == 1:
            q=(2.0*np.sqrt(q_z(LAMBDA, NK_Array[0:-1], theta0, NK_Array[0])*q_z(LAMBDA, NK_Array[1::], theta0, NK_Array[0])))
            corr = w(q, sigma, self.interface, theta0, NK_Array[0])
            rijs = rijs*corr
            rijp = rijp*corr
        elif self.SigmaCorrMethod == 0:
            q=(2.0*q_z(LAMBDA, NK_Array[0:-1], theta0, NK_Array[0]))
            corr = w(q, sigma, self.interface, theta0, NK_Array[0])
            rijs = rijs*corr
            rijp = rijp*corr
        elif self.SigmaCorrMethod == 2:
            q1 = q_z(LAMBDA, NK_Array[0:-1], theta0, NK_Array[0])
            q2 = q_z(LAMBDA, NK_Array[1::], theta0, NK_Array[0])
            corr = w(q1+q2, sigma, self.interface, theta0, NK_Array[0]) / w(q1-q2, sigma, self.interface, theta0, NK_Array[0])
            rijs = rijs*corr
            rijp = rijp*corr

        #phase[0] -> for 1st layer || phase[n-1] -> bottom most layer.
        phase = fi(LAMBDA, Z_Array, NK_Array[1:-1], theta0, NK_Array[0])
        #print("---Time in each Theta loop %s seconds ---" % (time.time() - start_time_th))
        dphi = np.exp(phase[n-1]*complex(0, 1))
        R_s = (rijs[n-1]+(rijs[n]*dphi))/(1.+(rijs[n-1]*rijs[n]*dphi))  # ;n-1 -> this is the last layer above the substrate
        R_p = (rijp[n-1]+(rijp[n]*dphi))/(1.+(rijp[n-1]*rijp[n]*dphi))

        if n == 2 :
            dphi = np.exp(phase[0]*complex(0, 1))
            R_s = (rijs[0]+(R_s*dphi))/(1.0+(rijs[0]*R_s*dphi))
            R_p = (rijp[0]+(R_p*dphi))/(1.0+(rijp[0]*R_p*dphi))
        elif n >= 3:
            for l in range(n-2, -1, -1):
                dphi = np.exp(phase[l]*complex(0, 1))
                R_s = (rijs[l]+(R_s*dphi))/(1.0+(rijs[l]*R_s*dphi))
                R_p = (rijp[l]+(R_p*dphi))/(1.0+(rijp[l]*R_p*dphi))

        projection_factor = projection(self.ProjEffCorr,IncidentBeamDia,self.SampleSize,theta0)
        Reflectivity_s = np.real(R_s*np.conj(R_s))*projection_factor
        Reflectivity_p = np.real(R_p*np.conj(R_p))*projection_factor
        Reflectivity_a = ((Reflectivity_s*self.DetectorPA*(1.0+self.BeamPolarization)) +
                      (Reflectivity_p*(1.0-self.BeamPolarization)))/((self.BeamPolarization*(self.DetectorPA-1.0)) + (self.DetectorPA+1.0))
        
        return Reflectivity_s,Reflectivity_p,Reflectivity_a


    def get_refl_tran_absp(self,Z_Array,NK_Array,sigma,theta_i,energy_i,IncidentBeamDia=0.1):
        n=self.n
        theta0 = theta_i
        LAMBDA = energy_i
        
        rijs = rij_s(NK_Array[0:-1], NK_Array[1::], theta0, NK_Array[0])
        rijp = rij_p(NK_Array[0:-1], NK_Array[1::], theta0, NK_Array[0])
        tijs = tij_s(NK_Array[0:-1], NK_Array[1::], theta0, NK_Array[0])
        tijp = tij_p(NK_Array[0:-1], NK_Array[1::], theta0, NK_Array[0])
        #q1 = q_z(LAMBDA, NK_Array[0:-1], theta0, NK_Array[0])
        #q2 = q_z(LAMBDA, NK_Array[1::], theta0, NK_Array[0])
        #corr_r = 1.0 / w(q1-q2, sigma, self.interface, theta0, NK_Array[0])
        if  self.SigmaCorrMethod == 1:
            q=(2.0*np.sqrt(q_z(LAMBDA, NK_Array[0:-1], theta0, NK_Array[0])*q_z(LAMBDA, NK_Array[1::], theta0, NK_Array[0])))
            corr = w(q, sigma, self.interface, theta0, NK_Array[0])
            corr_r = w(q, sigma, self.interface, theta0, NK_Array[0])
            rijs = rijs*corr
            rijp = rijp*corr
            tijs = tijs*corr_r
            tijp = tijp*corr_r
        elif self.SigmaCorrMethod == 0:
            q=(2.0*q_z(LAMBDA, NK_Array[0:-1], theta0, NK_Array[0]))
            corr = w(q, sigma, self.interface, theta0, NK_Array[0])
            corr_r = w(q, sigma, self.interface, theta0, NK_Array[0])
            rijs = rijs*corr
            rijp = rijp*corr
            tijs = tijs*corr_r
            tijp = tijp*corr_r
        elif self.SigmaCorrMethod == 2:
            q1 = q_z(LAMBDA, NK_Array[0:-1], theta0, NK_Array[0])
            q2 = q_z(LAMBDA, NK_Array[1::], theta0, NK_Array[0])
            corr = w(q1+q2, sigma, self.interface, theta0, NK_Array[0]) / w(q1-q2, sigma, self.interface, theta0, NK_Array[0])
            corr_r = 1.0 / w(q1-q2, sigma, self.interface, theta0, NK_Array[0])
            rijs = rijs*corr
            rijp = rijp*corr
            tijs = tijs*corr_r
            tijp = tijp*corr_r

        #phase[0] -> for 1st layer || phase[n-1] -> bottom most layer.
        phase = fi(LAMBDA, Z_Array, NK_Array[1:-1], theta0, NK_Array[0])
        #print("---Time in each Theta loop %s seconds ---" % (time.time() - start_time_th))
        ph_coff = phase[n-1]*complex(0, 1)
        dphir = np.exp(0.5*ph_coff)
        dphi = np.exp(ph_coff)
        R_s = (rijs[n-1]+(rijs[n]*dphi))/(1.+(rijs[n-1]*rijs[n]*dphi))
        R_p = (rijp[n-1]+(rijp[n]*dphi))/(1.+(rijp[n-1]*rijp[n]*dphi))
        T_s = (tijs[n-1]*(tijs[n]*dphir))/(1.+(rijs[n-1]*rijs[n]*dphi))
        T_p = (tijp[n-1]*(tijp[n]*dphir))/(1.+(rijp[n-1]*rijp[n]*dphi))
        if n == 2 :
            ph_coff = phase[0]*complex(0, 1)
            dphir = np.exp(0.5*ph_coff)
            dphi = np.exp(ph_coff)
            R_s = (rijs[0]+(R_s*dphi))/(1.0+(rijs[0]*R_s*dphi))
            R_p = (rijp[0]+(R_p*dphi))/(1.0+(rijp[0]*R_p*dphi))
            T_s = (tijs[0]*T_s*dphir)/(1.0+(rijs[0]*R_s*dphi))
            T_p = (tijp[0]*T_p*dphir)/(1.0+(rijp[0]*R_p*dphi))
        elif n >= 3:
            R_s_old=0
            R_p_old=0
            for l in range(n-2, -1, -1):
                ph_coff = phase[l]*complex(0, 1)
                dphir = np.exp(0.5*ph_coff)
                dphi = np.exp(ph_coff)
                R_s = (rijs[l]+(R_s*dphi))/(1.0+(rijs[l]*R_s*dphi))
                R_p = (rijp[l]+(R_p*dphi))/(1.0+(rijp[l]*R_p*dphi))
                T_s = (tijs[l]*T_s*dphir)/(1.0+(rijs[l]*R_s_old*dphi))
                T_p = (tijp[l]*T_p*dphir)/(1.0+(rijp[l]*R_p_old*dphi))
                R_s_old=R_s
                R_p_old=R_p
        #insres=gauss(Theta,self.InsRes)
        projection_factor = projection(self.ProjEffCorr,IncidentBeamDia,self.SampleSize,theta0)
        Reflectivity_s = np.real(R_s*np.conj(R_s))*projection_factor
        Reflectivity_p = np.real(R_p*np.conj(R_p))*projection_factor
        Reflectivity_a = ((Reflectivity_s*self.DetectorPA*(1.0+self.BeamPolarization)) +
                      (Reflectivity_p*(1.0-self.BeamPolarization)))/((self.BeamPolarization*(self.DetectorPA-1.0)) + (self.DetectorPA+1.0))

        proj_coeff=((NK_Array[n+1]*np.sqrt(1.0-(NK_Array[0]*np.cos(theta0)/NK_Array[n+1])**2))/(NK_Array[0]*np.sin(theta0)))
        #proj_coeff=np.real((NK_Array[n+1]*np.sqrt(1.0-(NK_Array[0]*np.cos(theta0)/NK_Array[n+1])**2))/(NK_Array[0]*np.sin(theta0)))
        #Transmittivity_s = np.real(T_s*np.conj(T_s))*proj_coeff
        #Transmittivity_p = np.real(T_p*np.conj(T_p))*proj_coeff
        Transmittivity_s = np.real((T_s*np.conj(T_s))*proj_coeff)
        Transmittivity_p = np.real((T_p*np.conj(T_p))*proj_coeff)
        Transmittivity_a = ((Transmittivity_s*self.DetectorPA*(1.0+self.BeamPolarization)) +
                      (Transmittivity_p*(1.0-self.BeamPolarization)))/((self.BeamPolarization*(self.DetectorPA-1.0)) + (self.DetectorPA+1.0))

        Absorptance_s = 1.0 - Reflectivity_s - Transmittivity_s
        Absorptance_p = 1.0 - Reflectivity_p - Transmittivity_p
        Absorptance_a = ((Absorptance_s*self.DetectorPA*(1.0+self.BeamPolarization)) +
                         (Absorptance_p*(1.0-self.BeamPolarization)))/((self.BeamPolarization*(self.DetectorPA-1.0)) + (self.DetectorPA+1.0))

        return Reflectivity_s,Reflectivity_p,Reflectivity_a,Transmittivity_s,Transmittivity_p,Transmittivity_a,Absorptance_s,Absorptance_p,Absorptance_a

    
