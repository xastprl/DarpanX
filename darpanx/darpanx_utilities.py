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

import numpy as np
import scipy.constants as const
from tabulate import tabulate
from os import path
from darpanx.get_dir import*
import sys
import multiprocessing as mlp
'''
#-----------------------------------------------
Perpose:
    Different functions to the Calculate, calling by main 
    DarpanX class.
Modification History:
    June.27.2019,... Biswajit
    May.31.2020,.... Biswajit

#-----------------------------------------------
'''

#;==============================================================================================

#;----------------------------------------------------------------------------------------------
def rij_s(ni,nj,theta0,n0):
    a=(n0*np.cos(theta0))
    x=ni*(np.sqrt(1.0-((a/ni)**2)))
    y=nj*(np.sqrt(1.0-((a/nj)**2)))
    return ((x-y)/(x+y))
#;-----------------------------------------------------------------------------------------------
def rij_p(ni,nj,theta0,n0):
    a=(n0*np.cos(theta0))
    x=nj*(np.sqrt(1.0-((a/ni)**2)))
    y=ni*(np.sqrt(1.0-((a/nj)**2)))
    return ((x-y)/(x+y))
#;-----------------------------------------------------------------------------------------------

def tij_s(ni,nj,theta0,n0):
    a=(n0*np.cos(theta0))
    x=ni*(np.sqrt(1.0-((a/ni)**2)))
    y=nj*(np.sqrt(1.0-((a/nj)**2)))
    return (2.0*x/(x+y))
#;-----------------------------------------------------------------------------------------------
def tij_p(ni,nj,theta0,n0):
    a=(n0*np.cos(theta0))
    x=nj*(np.sqrt(1.0-((a/ni)**2)))
    y=ni*(np.sqrt(1.0-((a/nj)**2)))
    return (2.0*x/(x+y))
#;-----------------------------------------------------------------------------------------------


#;=============================== Phase Function ================================================

def fi(LAMBDA,di,ni,theta0,n0):
    return (4.0*np.pi/LAMBDA)*di*ni*np.sqrt(1.0-(n0*np.cos(theta0)/ni)**2)

#--------------------------------------------------------------------------------------------
#;=============================== INTERFACE PROFILE FUNCTION=====================================

def w(s,sigma1,interface1,theta0,n0):
    #sigma1 should be an array.
    if (interface1 == 0) : return np.exp(-(s*sigma1)**2 / 2.0)
    elif (interface1 == 1) : return (1.0 / (1.0+(0.5*(s*sigma1)**2)))
    elif (interface1 == 2) : 
        a=s*sigma1*np.sqrt(3)
        a=np.array(a)
        b=np.array(np.sin(a))
        a_ind_0 = np.where(a == 0)
        a[a_ind_0] = 1
        b[a_ind_0] = 1
        return (b/a)
    elif (interface1 == 3) :
        sa=np.pi/np.sqrt(np.pi**2 - 8.0)

        a=(sa*sigma1*s)-(np.pi/2.0)
        b=(sa*sigma1*s)+(np.pi/2.0)
        return ((np.pi/4.0)*((np.sin(a)/(a))+(np.sin(b)/(b)))) 
    elif (interface1 == 4) : return np.cos(sigma1*s)
    #elif (sigma1 == 0) : return 1.0

def q_z(lambda1,ni,theta0,n0):
    return (2.0*np.pi/lambda1)*np.sqrt(1.0-(n0*np.cos(theta0)/ni)**2)


def projection(ProjEffCorr,IncidentBeamDia,SampleSize,theta):

    if ProjEffCorr == 'yes':
        slit_width=IncidentBeamDia
        Ai=slit_width/np.sin(theta)
        sample_area= SampleSize #; in mm in 1 Dimension
        if Ai < sample_area :
            Ar=Ai
        else:
            Ar=sample_area
        projection_factor=(Ar/Ai)
    else: projection_factor=1
    return projection_factor    

#Note: Can define another convolution function
#      for instrumental resolution...etc

#This Function converts angstrom units to kev
def a2kev(lambda1):
    ENERGY=np.array(lambda1)
    ENERGY=const.h*const.c/(lambda1*1.0e-10*const.e*1000.0)
    return ENERGY

#This Function converts KeV to angstrom unit.
def kev2a(lambda1):
    ENERGY=np.array(lambda1)
    LAMBDA=const.h*const.c/(ENERGY*1.0e-10*const.e*1000.0)
    return LAMBDA

def deg2rad(theta):
    return np.array(theta)*(const.pi/180.0)

def rad2deg(theta):
    return np.array(theta)*(180.0/const.pi)

def gauss(x,sd):
    s=2.0*sd**2
    gaus=exp(-x**2/s)/(np.sqrt(3.14*s))

def intro():
    print('=======================================================')
    print('  ')
    print('---------- Welcome to DarpanX Project Code ------------')
    print('------ Developed at Physical Research Laboratory-------')
    print('  ')
    print('=======================================================')
    print("  ")
    return None

def specific(MultilayerType=None ,LayerMaterial=None, Repetition=None, NumStack=None,
        EachRepetitionLayerNum=None, LayerNum=None, Period=None,Gamma=None,D_max=None,
        D_min=None,C=None,GammaTop=None, SigmaValues=None,TopLayerThick=None, EachRepetitionLayerThick=None, 
        Z_Array=None,NumCore=None,AmbientDensity=None,SubstrateDensity=None,
        SubstrateMaterial=None, TopLayerNum=None,AmbientMaterial=None,DensityCorrection=None, 
        InterfaceProf=None,SigmaCorrMethod=None, BeamPolarization=None,DetectorPA=None, ProjEffCorr=None, 
        SampleSize=None,IncidentBeamDia=None,LayerDensity=None, InsRes=None):

    pdtabulate=lambda df:tabulate(df,headers=['MultilayerType',MultilayerType],tablefmt='psql')
    if AmbientDensity is None: AmbientDensity='Default'
    if SubstrateDensity is None: SubstrateDensity='Default'

    if (MultilayerType == 'UserDefinedML'):
        if DensityCorrection == 'yes':
            list = [['Repetition',Repetition], ['TopLayerNum', str(TopLayerNum)], ['EachRepetitionLayerNum', str(EachRepetitionLayerNum)], ['EachRepetitionLayerThick', str(EachRepetitionLayerThick)], ['LayerMaterial',str(LayerMaterial)], ['DensityCorrection',DensityCorrection], ['LayerDensity',str(LayerDensity)], ['AmbientMaterial',AmbientMaterial], ['SubstrateMaterial', SubstrateMaterial],['AmbientDensity',AmbientDensity],['SubstrateDensity',SubstrateDensity],['SigmaValues',SigmaValues]]
        else: list = [['Repetition',Repetition], ['TopLayerNum', str(TopLayerNum)], ['EachRepetitionLayerNum', str(EachRepetitionLayerNum)], ['EachRepetitionLayerThick', str(EachRepetitionLayerThick)], ['LayerMaterial',str(LayerMaterial)], ['DensityCorrection',DensityCorrection], ['AmbientMaterial',AmbientMaterial], ['SubstrateMaterial', SubstrateMaterial],['SigmaValues',SigmaValues]]

    elif (MultilayerType == 'BiLayer'):
        if DensityCorrection == 'yes':
            list = [['Repetition',Repetition], ['Period', str(Period)], ['Gamma', str(Gamma)], ['LayerMaterial',str(LayerMaterial)], ['DensityCorrection',DensityCorrection], ['LayerDensity',str(LayerDensity)], ['AmbientMaterial',AmbientMaterial], ['SubstrateMaterial', SubstrateMaterial],['AmbientDensity',AmbientDensity],['SubstrateDensity',SubstrateDensity],['SigmaValues',SigmaValues]]
        else:list = [['Repetition',Repetition], ['Period', str(Period)], ['Gamma', str(Gamma)], ['LayerMaterial',str(LayerMaterial)], ['DensityCorrection',DensityCorrection], ['AmbientMaterial',AmbientMaterial], ['SubstrateMaterial', SubstrateMaterial],['AmbientDensity',AmbientDensity],['SubstrateDensity',SubstrateDensity],['SigmaValues',SigmaValues]]
 
    elif (MultilayerType == 'DepthGraded'):
        if DensityCorrection == 'yes':
            list = [['Repetition',Repetition],['D_max',str(D_max)],['D_min',str(D_min)],['C',str(C)], ['Period', str(Period)], ['GammaTop', str(GammaTop)], ['Gamma', str(Gamma)],['LayerMaterial',str(LayerMaterial)], ['DensityCorrection',DensityCorrection], ['LayerDensity',str(LayerDensity)], ['AmbientMaterial',AmbientMaterial], ['SubstrateMaterial', SubstrateMaterial],['AmbientDensity',AmbientDensity],['SubstrateDensity',SubstrateDensity],['SigmaValues',SigmaValues]]
        else:list = [['Repetition',Repetition],['D_max',str(D_max)],['D_min',str(D_min)],['C',str(C)], ['Period', str(Period)], ['GammaTop', str(GammaTop)], ['Gamma', str(Gamma)],['LayerMaterial',str(LayerMaterial)], ['DensityCorrection',DensityCorrection], ['AmbientMaterial',AmbientMaterial], ['SubstrateMaterial', SubstrateMaterial],['AmbientDensity',AmbientDensity],['SubstrateDensity',SubstrateDensity],['SigmaValues',SigmaValues]]

    elif (MultilayerType == 'ClusterGraded'):

        if DensityCorrection == 'yes':
            list = [['Repetition',Repetition],['NumStack',NumStack],['Period', str(Period)], ['Gamma', str(Gamma)],['LayerMaterial',str(LayerMaterial)], ['LayerDensity',str(LayerDensity)],['DensityCorrection',DensityCorrection], ['AmbientMaterial',AmbientMaterial], ['SubstrateMaterial', SubstrateMaterial],['AmbientDensity',AmbientDensity],['SubstrateDensity',SubstrateDensity],['SigmaValues',SigmaValues]]
        else: list = [['Repetition',Repetition],['NumStack',NumStack],['Period', str(Period)], ['Gamma', str(Gamma)],['LayerMaterial',str(LayerMaterial)], ['DensityCorrection',DensityCorrection], ['LayerDensity',str(LayerDensity)], ['AmbientMaterial',AmbientMaterial], ['SubstrateMaterial', SubstrateMaterial],['AmbientDensity',AmbientDensity],['SubstrateDensity',SubstrateDensity],['SigmaValues',SigmaValues]]

    elif (MultilayerType == 'SingleLayer'):
        if DensityCorrection == 'yes':
            list = [['Period', str(Period)], ['LayerMaterial',str(LayerMaterial)], ['DensityCorrection',DensityCorrection], ['LayerDensity',str(LayerDensity)], ['AmbientMaterial',AmbientMaterial], ['SubstrateMaterial', SubstrateMaterial],['AmbientDensity',AmbientDensity],['SubstrateDensity',SubstrateDensity],['SigmaValues',SigmaValues]]
        else:list = [['Period', str(Period)], ['LayerMaterial',str(LayerMaterial)], ['DensityCorrection',DensityCorrection], ['AmbientMaterial',AmbientMaterial], ['SubstrateMaterial', SubstrateMaterial],['AmbientDensity',AmbientDensity],['SubstrateDensity',SubstrateDensity],['SigmaValues',SigmaValues]]

    elif (MultilayerType == 'UserDefined'):
        if DensityCorrection == 'yes':
            list = [['LayerNum', str(LayerNum)],['Z_Array', str(Z_Array)], ['LayerMaterial',str(LayerMaterial)], ['DensityCorrection',DensityCorrection], ['LayerDensity',str(LayerDensity)], ['AmbientMaterial',AmbientMaterial], ['SubstrateMaterial', SubstrateMaterial],['AmbientDensity',AmbientDensity],['SubstrateDensity',SubstrateDensity],['SigmaValues',SigmaValues]]
        else:list = [['LayerNum', str(LayerNum)],['Z_Array', str(Z_Array)],['LayerMaterial',str(LayerMaterial)], ['DensityCorrection',DensityCorrection], ['AmbientMaterial',AmbientMaterial], ['SubstrateMaterial', SubstrateMaterial],['AmbientDensity',AmbientDensity],['SubstrateDensity',SubstrateDensity],['SigmaValues',SigmaValues]]
    else: list =[['0',0]]

    if ProjEffCorr == 'yes':list1 = [['InterfaceProf', InterfaceProf], ['SigmaCorrMethod',SigmaCorrMethod], ['ProjEffCorr',ProjEffCorr], ['SampleSize',str(SampleSize)], ['IncidentBeamDia', str(IncidentBeamDia)],['BeamPolarization',str(BeamPolarization)],['DetectorPA',str(DetectorPA)]]
    elif  ProjEffCorr != 'yes' and ProjEffCorr != None: list1 = [['InterfaceProf', InterfaceProf], ['SigmaCorrMethod',SigmaCorrMethod], ['ProjEffCorr',ProjEffCorr],['BeamPolarization',str(BeamPolarization)],['DetectorPA',str(DetectorPA)]]
    elif ProjEffCorr is None : list1=[['0',0]]
    list=list+list1
    print(' ')
    print(pdtabulate(list))
    del pdtabulate,list
    return None

def indepen_var(Xscan=None,Energy=None,Theta=None):

    pdtabulate=lambda df:tabulate(df,headers=['Independent Parameters',''],tablefmt='psql')
    if Xscan is None:
        if Energy is not None and len(Energy) == 1: Xscan ='Theta'
        elif Theta is not None and len(Theta) == 1: Xscan ='Energy'
    if Xscan =='Theta':list = [['Xscan',Xscan],['Energy (keV)',Energy]]
    elif Xscan =='Energy': list = [['Xscan',Xscan],['Theta (deg)',Theta]]
    else:list=[['0',0]]
    print(' ')
    print(pdtabulate(list))
    print(' ')
    del pdtabulate,list
    return None

def material_prop(LayerDensity=None,AmbientDensity=None,SubstrateDensity=None):
    pdtabulate=lambda df:tabulate(df,headers=['','Density Used in Calculation'],tablefmt='psql')
    list = [['AmbientDensity',AmbientDensity],['SubstrateDensity',SubstrateDensity],['LayerDensity',str(LayerDensity)]]
    print(' ')
    print(pdtabulate(list))
    print(' ')
    del pdtabulate,list
    return None

def check_error(MultilayerType=None ,LayerMaterial=None, Repetition=None, NumStack=None,
        EachRepetitionLayerNum=None, LayerNum=None, Period=None,Gamma=None,D_max=None,
        D_min=None,C=None,GammaTop=None, SigmaValues=None,TopLayerThick=None, EachRepetitionLayerThick=None,
        Z_Array=None,NumCore=None,AmbientDensity=None,SubstrateDensity=None,
        SubstrateMaterial=None, TopLayerNum=None,AmbientMaterial=None,DensityCorrection=None,
        InterfaceProf=None,SigmaCorrMethod=None, BeamPolarization=None,DetectorPA=None, ProjEffCorr=None,
        SampleSize=None,IncidentBeamDia=None,LayerDensity=None, InsRes=None):

    if ProjEffCorr == 'yes':
        try:
            SampleSize = float(SampleSize)
        except:raise Exception('%% DarpanX_Error: Incorrect < SampleSize >')
    try:
        if LayerMaterial is None or len(LayerMaterial) < 1:
            raise Exception('%% DarpanX_Error: Undefined < LayerMaterial > -- It should be a string array, like ["'"Si"'"] or ["'"Pt"'","'"C"'"]')
        else: pass
    except:raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial > -- It should be a string array, like ["'"Si"'"] or ["'"Pt"'","'"C"'"]')
    if InterfaceProf != 'ErrorFunction' and InterfaceProf != 'ExponentialFunction' and InterfaceProf != 'LinearFunction' and InterfaceProf != 'SinusoidalFunction' and InterfaceProf != 'StepFunction':
        raise Exception('%% DarpanX_Error: Undefined < InterfaceProf > --options are; < "'"ErrorFunction"'" >, < "'"ExponentialFunction"'" >, < "'"LinearFunction"'" >, < "'"SinusoidalFunction"'" >, < "'"StepFunction"'" > .')
    if SigmaCorrMethod != 'NevotCroce' and SigmaCorrMethod != 'DebyeWaller' and SigmaCorrMethod != 'Avg':
        raise Exception('%% DarpanX_Error: Undefined < SigmaCorrMethod > --options are; < "'"NevotCroce"'" >, < "'"DebyeWaller"'" >')
    
    if MultilayerType != None:
        if MultilayerType == 'UserDefinedML':
            try:
                TopLayerNum = int(TopLayerNum)
                EachRepetitionLayerNum = EachRepetitionLayerNum
                Repetition = int(Repetition)
            except: raise Exception("%% DarpanX_Error: Incorrect or Undefined  < TopLayerNum >, < EachRepetitionLayerNum >, < Repetition > for UserDefinedML structure")
            try:
                if (TopLayerNum + EachRepetitionLayerNum) != len(LayerMaterial):
                    raise Exception('%% DarpanX_Error: Inconsistent < TopLayerNum > and < EachRepetitionLayerNum > with the length(< LayerMaterial >)')
            except: raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial >')
            if Repetition < 1 :
                raise Exception('%% DarpanX_Error: Incorrect < Repetition >. It should be an integer number, like, 1,2,3...etc')
            if EachRepetitionLayerThick is None:
                raise Exception('%% DarpanX_Error: Undefined < EachRepetitionLayerThick >')
            elif len(EachRepetitionLayerThick) != EachRepetitionLayerNum :
                raise Exception('%% DarpanX_Error: len(EachRepetitionLayerThick){'+str(len(EachRepetitionLayerThick))+'} should be = < EachRepetitionLayerNum >{'+str(EachRepetitionLayerNum)+'}')
            if TopLayerNum > 0 and TopLayerThick is None:
                raise Exception('%% DarpanX_Error: Undefined < TopLayerThick >')
            elif TopLayerNum > 0 and len(TopLayerThick) != TopLayerNum:
                raise Exception('%% DarpanX_Error: len(TopLayerThick){'+str(len(TopLayerThick))+'} should be = < TopLayerNum >{'+str(TopLayerNum)+'}')
            if SigmaValues != None:
                try:
                    if len(SigmaValues) != 1 and len(SigmaValues) != (EachRepetitionLayerNum+TopLayerNum)+1 and len(SigmaValues) != EachRepetitionLayerNum+TopLayerNum+2:
                        if Repetition == 1:raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions.len(SigmaValues) should be- 1 or '+str(EachRepetitionLayerNum+TopLayerNum+1))
                        else:raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions. len(SigmaValues) should be- 1 or '+str(EachRepetitionLayerNum+TopLayerNum+2)+', or '+str(len(LayerMaterial)+1))
                except:raise Exception("%% DarpanX_Error: Incorect < SigmaValues >")

        elif MultilayerType == 'SingleLayer':
            try:
                if len(LayerMaterial) != 1:
                    raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial > -- should be a 1-element string array like ["'"Pt"'"]')
            except: raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial > -- should be a 1-element string array like ["'"Pt"'"]')
            try:
                if len(Period) != 1: raise Exception('%% DarpanX_Error: Incorrect or undefied < Period > ')
                Period = float(Period[0])
            except: raise Exception("%% DarpanX_Error: < Period > (float) is required for Single-layer Structure.")
            try:
                if len(SigmaValues) != 1 and len(SigmaValues) != 2:raise Exception('%% DarpanX_Error: Incorrect dimensions of < SigmaValues >. len(SigmaValues) should be- 1 or 2')
            except:raise Exception('%% DarpanX_Error: Undefined or incorrect < SigmaValues >. len(SigmaValues) should be- 1 or 2')

        elif MultilayerType == 'BiLayer':
            try:
                if len(LayerMaterial) != 2:
                    raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial > -- should be a 2-element array, like ["'"W"'","'"Si"'"].')
            except:raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial > -- should be a 2-element array, like ["'"W"'","'"Si"'"].')
            try:
                n = 2*Repetition
            except: raise Exception("%% DarpanX_Error: < Repetition > is required for Bilayer Structure.")
            try:
                Period = float(Period[0])
                Gamma = float(Gamma[0])
            except: raise Exception("%% DarpanX_Error: < Period >, < Gamma > are required for Bilayer Structure.")
            n=len(LayerMaterial) 
            if n > 2:
                if len(SigmaValues) != 1 and len(SigmaValues) != 4 and len(SigmaValues) != n+1:raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions len(SigmaValues) should be- 1 or 4 or '+str(n+1))
            elif n==2:
                if len(SigmaValues) != 1 and len(SigmaValues) != 3:raise Exception("%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions. len(SigmaValues) should be- 1 or 3")

        elif MultilayerType == 'DepthGraded':
            if len(LayerMaterial) != 2:
                raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial > -- should be a 2-element array, like ["'"W"'","'"Si"'"].')
            try:
                n = 2*Repetition
            except: raise Exception("%% DarpanX_Error: < Repetition > is required for DepthGraded Structure.")
            try:
                D_max=float(D_max)
                D_min=float(D_min)
                C=float(C)
            except: raise Exception("%% DarpanX_Error: < D_max >,< D_min >,< C >, < Gamma > all are required for DepthGraded Structure.")
            n=len(LayerMaterial)
            if n > 2:
                if len(SigmaValues) != 1 and len(SigmaValues) != 4 and len(SigmaValues) != n+1:raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions len(SigmaValues) should be- 1 or 4 or '+str(n+1))
            elif n==2:
                if len(SigmaValues) != 1 and len(SigmaValues) != 3:raise Exception("%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions. len(SigmaValues) should be- 1 or 3")
        
        elif MultilayerType == 'ClusterGraded':
            if len(LayerMaterial) != 2:
                raise Exception('%% DarpanX_Error: Incorrect < LayerMaterial > -- should be a 2-element array, like ["'"W"'","'"Si"'"].')
            try:
                NumStuck = int(NumStack)
                n = 2*Repetition
            except: raise Exception("%% DarpanX_Error: < Repetition >, < NumStack > are required for ClusterGraded Structure.")
            try:
                if len(Gamma) != NumStack: raise Exception("%% DarpanX_Error: Incorrect dimension of < Gamma >. Dimension should be- "+str(NumStuck))
                if len(Period) != NumStack: raise Exception("%% DarpanX_Error: Incorrect dimension of < Priod >. Dimension should be- "+str(NumStuck))
            except: raise Exception("%% DarpanX_Error: Undefined < Period >, < Gamma > or Incorrect dimension. Dimension should be- "+str(NumStuck))
            #n=len(LayerMaterial)
            if n > 2:
                if len(SigmaValues) != 1 and len(SigmaValues) != 4 and len(SigmaValues) != n+1 and len(SigmaValues) != (NumStuck*2+2):
                    raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions.')# len(SigmaValues) should be- 1 or 4 or '+str(n+1))
            elif n==2:
                if len(SigmaValues) != 1 and len(SigmaValues) != 3:raise Exception("%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions. len(SigmaValues) should be- 1 or 3")



        elif MultilayerType == 'UserDefined':
            try:
                LayerNum=int(LayerNum)
            except: raise Exception("%% DarpanX_Error: Undefined or incorrect  < LayerNum >.")
            if len(LayerMaterial) != LayerNum:
                raise Exception('%% DarpanX_Error: Incorrect dimension of  < LayerMaterial > -- should be a '+str(LayerNum)+'-element array, like ["'"W"'","'"Si"'",...].')
            #if DensityCorrection == 'yes' and len(OriginalDensity) != LayerNum:
            #    raise Exception('%% DarpanX_Error: Incorrect dimension of < OriginalDensity >..It should be='+str(len(LayerMaterial)))
            if len(SigmaValues) != len(LayerNum+1) and len(SigmaValues) !=1:raise Exception('%% DarpanX_Error: Undefined < SigmaValues > or incorrect dimensions. len(SigmaValues) should be- 1 or '+str(LayerNum+1))
        else:raise Exception('%% DarpanX_Error: Undefined < MultilayerType > --options are < ["'"UserDefinedML"'","'"SingleLayer"'", "'"DepthGraded"'", "'"ClusterGraded"'", "'"UserDefined"'"] >.')
    try:
        if DensityCorrection =='yes' and len(LayerDensity) != len(LayerMaterial):
            raise Exception('%% DarpanX_Error: Incorrect dimension of < LayerDensity >..It should be='+str(len(LayerMaterial)))
        else: pass
    except: raise Exception('%% DarpanX_Error: Incorrect or undefined dimension of < LayerDensity > or < LayerMaterial >.')

def err_check_read_infile(UserInFile):
    try:
        UserInFile=str(UserInFile)
    except:raise Exception("%% DarpanX_Error: Incorrect or undefined < UserInFile >. Define it like, UserInFile='inputs.drpnx' ")

    while path.exists(UserInFile) is False:
        print('%% DarpanX_Input: Input config file '+UserInFile+' not exist..')
        print('Give the input Config file name:')
        UserInFile = str(input())

    #print("%% DarpanX_message: Configuraton file (UserInFile) used:  "+UserInFile)
    return UserInFile

def read_infile(UserInFile):

    UserInFile=err_check_read_infile(UserInFile)
    with open(UserInFile) as userinputs:
        l1=[]
        l2=[]
        inp=[]
        for line in userinputs:
            try:
                if line.strip():
                    l=(line.replace(" ", "").strip().split('='))
                    inp=inp+[line]
                    if str(l[0]) and str(l[1]):
                        l1=l1+[l[0]]
                        l2=l2+[l[1]]
            except:pass
            if '# ======== Outputs ==========' in line:
                inp=inp[2:-1]
                break
    try:
        NK_dir=l2[l1.index('NK_dir')][1:-1]
    except:
        print("%% DarpanX_Settings: Undefined < NK_dir >. Set default direcrory (nk_data/nk)")
        dir1=get_dir()
        NK_dir=os.path.join(str(dir1),"..","nk_data/nist")
    Xscan=l2[l1.index('Xscan')][1:-1]
    if Xscan != 'Theta' and Xscan != 'Energy': raise Exception("%% DarpanX_Error: Undefined or Incorrect < Xscan > in Config file.-- options are [Theta] or [Energy] ")
    try:
        if Xscan == 'Theta':
            Energy=l2[l1.index('Energy')][1:-1]
            Energy=kev2a(float(Energy))
            Theta=None
        elif Xscan == 'Energy':
            Theta=l2[l1.index('Theta')][1:-1]
            Theta=deg2rad(float(Theta))
            Energy=None
    except:raise Exception("%% DarpanX_Error: Undefined or Incorrect incident beam energy or angle in Config file.")
    try:
        MultilayerType=l2[l1.index('MultilayerType')][1:-1]
    except:MultilayerType=None
    try:
        TopLayerNum=int(l2[l1.index('TopLayerNum')][1:-1])
    except:TopLayerNum=0
    try:
        EachRepetitionLayerNum=int(l2[l1.index('EachRepetitionLayerNum')][1:-1])
    except:EachRepetitionLayerNum=None
    try:
        LayerMaterial=l2[l1.index('LayerMaterial')]
        LayerMaterial=list(LayerMaterial[1:-1].split(','))
    except:LayerMaterial=None
    try:
        TopLayerThick=l2[l1.index('TopLayerThick')]
        TopLayerThick=list(TopLayerThick[1:-1].split(','))
        TopLayerThick=[float(i) for i in TopLayerThick]
    except:TopLayerThick=None
    try:
        EachRepetitionLayerThick=l2[l1.index('EachRepetitionLayerThick')]
        EachRepetitionLayerThick=list(EachRepetitionLayerThick[1:-1].split(','))
        EachRepetitionLayerThick=[float(i) for i in EachRepetitionLayerThick]
    except:EachRepetitionLayerThick=None
    try:
        Repetition=int(l2[l1.index('Repetition')][1:-1])
    except:Repetition=None
    try:
        SubstrateMaterial=l2[l1.index('SubstrateMaterial')][1:-1]
    except:SubstrateMaterial='Vacuum'
    try:
        AmbientMaterial=l2[l1.index('AmbientMaterial')][1:-1]
    except:AmbientMaterial='Vacuum'
    try:
        LayerNum=int(l2[l1.index('LayerNum')][1:-1])
    except:LayerNum=None
    try:
        NumStack=int(l2[l1.index('NumStack')][1:-1])
    except:NumStack=None
    try:
        Period=l2[l1.index('Period')]
        Period=list(Period[1:-1].split(','))
        Period=[float(i) for i in Period]
    except:Period=None
    try:
        Gamma=l2[l1.index('Gamma')]
        Gamma=list(Gamma[1:-1].split(','))
        Gamma=[float(i) for i in Gamma]
    except:Gamma=None
    try:
        D_max=float(l2[l1.index('D_max')][1:-1])
    except:D_max=None
    try:
        D_min=float(l2[l1.index('D_min')][1:-1])
    except:D_min=None
    try:
        C=float(l2[l1.index('C')][1:-1])
    except:C=None
    try:
        GammaTop=float(l2[l1.index('GammaTop')][1:-1])
    except:GammaTop=None
    try:
        Z_Array=l2[l1.index('Z_Array')]
        Z_Array=list(Z_Array[1:-1].split(','))
        Z_Array=[float(i) for i in Z_Array]
    except:Z_Array=None
    try:
        SigmaValues=l2[l1.index('SigmaValues')]
        SigmaValues=list(SigmaValues[1:-1].split(','))
        SigmaValues=[float(i) for i in SigmaValues]
    except:
        print("%% DarpanX_Settings: Undefined < SigmaValues >. Set SigmaValues=[0.0]")
        SigmaValues=[0.0]
    try:
        SigmaCorrMethod=l2[l1.index('SigmaCorrMethod')[1:-1]]
    except:SigmaCorrMethod='NevotCroce'
    try:
        InterfaceProf=l2[l1.index('InterfaceProf')][1:-1]
    except:InterfaceProf='ErrorFunction'
    try:
        ProjEffCorr=l2[l1.index('ProjEffCorr')][1:-1]
    except:ProjEffCorr='no'
    try:
        SampleSize=float(l2[l1.index('SampleSize')][1:-1])
    except:SampleSize=100
    try:
        IncidentBeamDia=float(l2[l1.index('IncidentBeamDia')][1:-1])
    except:IncidentBeamDia=0.1
    try:
        BeamPolarization=int(l2[l1.index('BeamPolarization')][1:-1])
    except:
        print("%% DarpanX_Settings: Incorrect or Undefined < BeamPolarization >. Set to 0")
        BeamPolarization=0
    try:
        DetectorPA=int(l2[l1.index('DetectorPA')][1:-1])
    except:
        print("%% DarpanX_Settings: Incorrect or Undefined < DetectorPA >. Set to 1")
        DetectorPA=1
    try:
        InsRes=float(l2[l1.index('InsRes')][1:-1])
    except:
        print("%% DarpanX_Settings: Incorrect or Undefined < InsRes >. Set to 0.0")
        InsRes=0
    try:
        DensityCorrection=l2[l1.index('DensityCorrection')][1:-1]
    except:DensityCorrection='no'
    #try:
    #    OriginalDensity=l2[l1.index('OriginalDensity')]
    #    OriginalDensity=list(OriginalDensity[1:-1].split(','))
    #    OriginalDensity=[float(i) for i in OriginalDensity]
    #except:OriginalDensity=[]
    try:
        LayerDensity=l2[l1.index('LayerDensity')]
        LayerDensity=list(LayerDensity[1:-1].split(','))
        LayerDensity=[float(i) for i in LayerDensity]
    except:LayerDensity=[]
    try:
        SubstrateDensity=l2[l1.index('SubstrateDensit')][1:-1]
        SubstrateDensity=float(SubstrateDensit)
    except:SubstrateDensity=None
    try:
        AmbientDensity=l2[l1.index('AmbientDensity')][1:-1]
        AmbientDensity=float(AmbientDensity)
    except:AmbientDensity=None
    try:
        NumCore=l2[l1.index('NumCore')][1:-1]
        NumCore=int(NumCore)
        if NumCore > 0:
            cpu=mlp.cpu_count()
            if NumCore > cpu:
                print("%% DarpanX_message: NumCore is exciding the maximum value. Set NumCore = "+str(cpu))
                NumCore=cpu
            print("%% DarpanX_message: Parallel processing is using with no.cores = "+str(NumCore))
        else: print("%% DarpanX_Error: < NumCore > should be an integer, like 1,2,3...")
    except:NumCore=None
    try:
        ShowPar=l2[l1.index('ShowPar')][1:-1]
        ShowPar=str(ShowPar)
    except:ShowPar = 'yes'

    #if MultilayerType == 'UserDefined':
    #    print("%% DarpanX_Input: Give the input Z_Array and Sigmavalues from ASCII files")
    #else: NK_Array=None
    read_infile.NK_dir=NK_dir
    read_infile.Xscan = Xscan
    read_infile.Energy = Energy
    read_infile.Theta = Theta
    read_infile.MultilayerType = MultilayerType
    read_infile.AmbientMaterial = AmbientMaterial
    read_infile.SubstrateMaterial = SubstrateMaterial
    read_infile.LayerMaterial = LayerMaterial
    read_infile.Repetition = Repetition
    read_infile.TopLayerNum = TopLayerNum
    read_infile.TopLayerThick = TopLayerThick
    read_infile.EachRepetitionLayerNum = EachRepetitionLayerNum
    read_infile.EachRepetitionLayerThick = EachRepetitionLayerThick
    read_infile.Period = Period
    read_infile.Gamma = Gamma
    read_infile.D_max = D_max
    read_infile.D_min = D_min
    read_infile.C = C
    read_infile.GammaTop = GammaTop
    read_infile.NumStack = NumStack
    read_infile.LayerNum = LayerNum
    read_infile.SigmaValues = SigmaValues
    read_infile.DensityCorrection = DensityCorrection
    read_infile.InterfaceProf = InterfaceProf
    read_infile.SigmaCorrMethod = SigmaCorrMethod
    read_infile.BeamPolarization= BeamPolarization
    read_infile.DetectorPA = DetectorPA
    read_infile.ProjEffCorr = ProjEffCorr
    read_infile.SampleSize = SampleSize
    read_infile.IncidentBeamDia = IncidentBeamDia
    #read_infile.OriginalDensity = OriginalDensity
    read_infile.LayerDensity = LayerDensity
    read_infile.InsRes = InsRes
    read_infile.NumCore = NumCore
    read_infile.Z_Array = Z_Array
    #read_infile.NK_Array = NK_Array
    read_infile.AmbientDensity = AmbientDensity
    read_infile.SubstrateDensity = SubstrateDensity
    read_infile.ShowPar = ShowPar
    read_infile.inp = inp

def find_duplicate(x):
    x=np.array(x)
    res = []
    for i in x:
        if i not in res:
            res.append(i)
    return res
