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
import numpy as np
from xspec import*
import sys
from os import path
import time
import multiprocessing as mulp
import matplotlib.pyplot as plt
#from tabulate import tabulate
import datetime
from darpanx.get_dir import*

UserInFile=sys.argv
err_check_read_infile(UserInFile)
print('  ')
print('=======================================================')
print("  ")
print("%% DarpanX_message: Configuraton file used:  "+UserInFile)

read_infile(UserInFile)

inp=read_infile.inp
NK_dir=read_infile.NK_dir
Xscan =read_infile. Xscan
Energy = read_infile.Energy
Theta = read_infile.Theta
MultilayerType = read_infile.MultilayerType
AmbientMaterial = read_infile.AmbientMaterial
SubstrateMaterial = read_infile.SubstrateMaterial
LayerMaterial = read_infile.LayerMaterial
Repetition = read_infile.Repetition
TopLayerNum = read_infile.TopLayerNum
TopLayerThick = read_infile.TopLayerThick
EachRepetitionLayerNum = read_infile.EachRepetitionLayerNum
EachRepetitionLayerThick = read_infile.EachRepetitionLayerThick
Period = read_infile.Period
Gamma = read_infile.Gamma
D_max = read_infile.D_max
D_min = read_infile.D_min
C = read_infile.C
GammaTop = read_infile.GammaTop
NumStack = read_infile.NumStack
LayerNum = read_infile.LayerNum
SigmaValues = read_infile.SigmaValues
DensityCorrection = read_infile.DensityCorrection
InterfaceProf = read_infile.InterfaceProf
SigmaCorrMethod = read_infile.SigmaCorrMethod
BeamPolarization= read_infile.BeamPolarization
DetectorPA = read_infile.DetectorPA
ProjEffCorr = read_infile.ProjEffCorr
SampleSize = read_infile.SampleSize
IncidentBeamDia = read_infile.IncidentBeamDia
#OriginalDensity = read_infile.OriginalDensity
Density = read_infile.LayerDensity
InsRes = read_infile.InsRes
NumCore = read_infile.NumCore
Z_Array = read_infile.Z_Array
#NK_Array = read_infile.NK_Array
AmbientDensity = read_infile.AmbientDensity
SubstrateDensity = read_infile.SubstrateDensity
ShowPar = read_infile.ShowPar

check_error(MultilayerType =MultilayerType ,LayerMaterial=LayerMaterial, Repetition=Repetition, NumStack=NumStack,
        EachRepetitionLayerNum=EachRepetitionLayerNum, LayerNum=LayerNum, Period=Period,Gamma=Gamma,D_max=D_max,
        D_min=D_min,C=C,GammaTop=GammaTop, SigmaValues=SigmaValues,TopLayerThick=TopLayerThick, EachRepetitionLayerThick=EachRepetitionLayerThick,
        Z_Array=Z_Array,NumCore=NumCore,AmbientDensity=AmbientDensity,SubstrateDensity=SubstrateDensity,
        SubstrateMaterial=SubstrateMaterial, TopLayerNum=TopLayerNum,AmbientMaterial=AmbientMaterial,DensityCorrection=DensityCorrection,
        InterfaceProf=InterfaceProf,SigmaCorrMethod=SigmaCorrMethod, BeamPolarization=BeamPolarization,DetectorPA=DetectorPA, ProjEffCorr=ProjEffCorr,
        SampleSize=SampleSize,IncidentBeamDia=IncidentBeamDia,LayerDensity=Density, InsRes=InsRes)

if ShowPar == 'yes':
    specific(MultilayerType =MultilayerType ,LayerMaterial=LayerMaterial, Repetition=Repetition, NumStack=NumStack,
        EachRepetitionLayerNum=EachRepetitionLayerNum, LayerNum=LayerNum, Period=Period,Gamma=Gamma,D_max=D_max,
        D_min=D_min,C=C,GammaTop=GammaTop, SigmaValues=SigmaValues,TopLayerThick=TopLayerThick, EachRepetitionLayerThick=EachRepetitionLayerThick,
        Z_Array=Z_Array,NumCore=NumCore,AmbientDensity=AmbientDensity,SubstrateDensity=SubstrateDensity,
        SubstrateMaterial=SubstrateMaterial, TopLayerNum=TopLayerNum,AmbientMaterial=AmbientMaterial,DensityCorrection=DensityCorrection,
        InterfaceProf=InterfaceProf,SigmaCorrMethod=SigmaCorrMethod, BeamPolarization=BeamPolarization,DetectorPA=DetectorPA, ProjEffCorr=ProjEffCorr,
        SampleSize=SampleSize,IncidentBeamDia=IncidentBeamDia,LayerDensity=Density, InsRes=InsRes)
    if Xscan == 'Theta':
        EneR=a2kev(Energy)
        TheR=None
    elif Xscan == 'Energy':
        TheR=rad2deg(Theta)
        EneR=None
    indepen_var(Xscan=Xscan,Theta=TheR,Energy=EneR)
print('  ')
print('=======================================================')
print("  ")
a=ML_Structure(NK_dir=NK_dir, LayerMaterial=LayerMaterial, Repetition=Repetition,NumStack=NumStack, MultilayerType=MultilayerType,
                    EachRepetitionLayerNum=EachRepetitionLayerNum,LayerNum=LayerNum, SubstrateMaterial=SubstrateMaterial, TopLayerNum
                    =TopLayerNum,AmbientMaterial=AmbientMaterial,DensityCorrection=DensityCorrection, InterfaceProf=InterfaceProf, 
                    SigmaCorrMethod=SigmaCorrMethod, BeamPolarization=BeamPolarization,DetectorPA=DetectorPA, ProjEffCorr=ProjEffCorr,
                    SampleSize=SampleSize, InsRes=InsRes)

#material_prop(LayerDensity=Density,AmbientDensity=a.AmbientDensity,SubstrateDensity=a.SubstrateDensity)

if Xscan == 'Theta': 
    ncsl=a.Materials_NK(Energy)
    if DensityCorrection != 'yes':NK_Array=a.refractiveindexarray(ncsl,AmbientDensity=AmbientDensity,SubstrateDensity=SubstrateDensity)

def theta_scan(i):
    Rs,Rp,Ra = a.get_reflectivity(z,NK_Array,sigma,Theta_in[i],Energy,IncidentBeamDia=IncidentBeamDia)
    return Rs,Rp,Ra

def energy_scan(i):
    global NK_Array
    ncsl=a.Materials_NK(Energy_in[i])
    NK_Array=a.refractiveindexarray(ncsl,LayerDensity=Density,AmbientDensity=AmbientDensity,SubstrateDensity=SubstrateDensity)
    Rs,Rp,Ra = a.get_reflectivity(z,NK_Array,sigma,Theta,Energy_in[i],IncidentBeamDia=IncidentBeamDia)
    return Rs,Rp,Ra

def darpanx_Th(engs, params, flux):
    global Theta_in,Period,Gamma,D_max,D_min,C,GammaTop,TopLayerThick,EachRepetitionLayerThick, Z_Array, NK_Array, z, sigma,IncidentBeamDia
    bin_num = len(engs)-1
    Theta_in=deg2rad(np.array(engs))
    
    if MultilayerType == 'UserDefinedML':
        TopLayerThick = params[0:TopLayerNum]
        indp=TopLayerNum+EachRepetitionLayerNum
        EachRepetitionLayerThick = params[TopLayerNum:indp]
        if DensityCorrection == 'yes':
            indp2=indp+len(LayerMaterial)
            Density = params[indp:indp2]
            NK_Array=a.refractiveindexarray(ncsl,LayerDensity=Density,AmbientDensity=AmbientDensity,SubstrateDensity=SubstrateDensity)
            if ProjEffCorr == 'yes':
                IncidentBeamDia=params[indp2]
                SigmaValues=params[indp2+1:-1]
            else: SigmaValues=params[indp2:-1]
        elif ProjEffCorr == 'yes':
            IncidentBeamDia=params[indp]
            SigmaValues=params[indp+1:-1]
        else: SigmaValues=params[indp:-1]

    elif MultilayerType == 'SingleLayer':
        Period = params[0]
        if DensityCorrection == 'yes': 
            Density = [params[1]]
            NK_Array=a.refractiveindexarray(ncsl,LayerDensity=Density,AmbientDensity=AmbientDensity,SubstrateDensity=SubstrateDensity)
            if ProjEffCorr == 'yes':
                IncidentBeamDia=params[2]
                SigmaValues=params[3:-1]
            else: SigmaValues=params[2:-1]
        elif ProjEffCorr == 'yes':
            IncidentBeamDia=params[1]
            SigmaValues=params[2:-1]
        else: SigmaValues=params[1:-1]

    elif MultilayerType == 'BiLayer':
        Period = params[0]
        Gamma = params[1]
        if DensityCorrection == 'yes': 
            Density = params[2:4]
            NK_Array=a.refractiveindexarray(ncsl,LayerDensity=Density,AmbientDensity=AmbientDensity,SubstrateDensity=SubstrateDensity)
            if ProjEffCorr == 'yes':
                IncidentBeamDia=params[4]
                SigmaValues=params[5:-1]
            else:SigmaValues=params[4:-1]
        elif ProjEffCorr == 'yes':
            IncidentBeamDia=params[2]
            SigmaValues=params[3:-1]
        else: SigmaValues=params[2:-1]

    elif MultilayerType == 'ClusterGraded':
        indp=NumStack*2
        Period = params[0:indp:2]
        Gamma = params[1:indp:2]
        if DensityCorrection == 'yes': 
            indp2=indp+2
            Density = params[indp:indp2]
            NK_Array=a.refractiveindexarray(ncsl,LayerDensity=Density,AmbientDensity=AmbientDensity,SubstrateDensity=SubstrateDensity)
            if ProjEffCorr == 'yes':
                IncidentBeamDia=params[indp2]
                SigmaValues=params[indp2+1:-1]
            else:SigmaValues=params[indp2:-1]
        elif ProjEffCorr == 'yes':
            IncidentBeamDia=params[indp]
            SigmaValues=params[indp+1:-1]
        else: SigmaValues=params[indp:-1]

    elif MultilayerType == 'DepthGraded':
        D_max = params[0]
        D_min = params[1]
        GammaTop = params[2]
        Gamma = params[3]
        C = params[4]
        if DensityCorrection == 'yes':
            indp=4+len(LayerMaterial)+1
            Density = params[5:indp]
            NK_Array=a.refractiveindexarray(ncsl,LayerDensity=Density,AmbientDensity=AmbientDensity,SubstrateDensity=SubstrateDensity)
            if ProjEffCorr == 'yes':
                IncidentBeamDia=params[indp]
                SigmaValues=params[indp+1:-1]
            else:SigmaValues=params[indp:-1]
        elif ProjEffCorr == 'yes':
            IncidentBeamDia=params[5]
            SigmaValues=params[6:-1]
        else: SigmaValues=params[5:-1]

    elif MultilayerType == 'UserDefined':
        Z_Array = params[0:LayerNum]
        indp=2*LayerNum
        if DensityCorrection == 'yes': 
            Density = params[LayerNum:indp]
            NK_Array=a.refractiveindexarray(ncsl,LayerDensity=Density,AmbientDensity=AmbientDensity,SubstrateDensity=SubstrateDensity)
            indp2=indp+LayerNum
            if ProjEffCorr == 'yes':
                IncidentBeamDia=params[indp2]
                SigmaValues=params[indp2+1:-1]
            else:SigmaValues=params[indp2:-1]
        elif ProjEffCorr == 'yes':
            IncidentBeamDia=params[indp]
            SigmaValues=params[indp+1:-1]
        else: SigmaValues=params[indp:-1]
    
    #if DensityCorrection is 'yes': NK_Array=a.refractiveindexarray(ncsl,LayerDensity=Density,AmbientDensity=AmbientDensity,SubstrateDensity=SubstrateDensity)

    z=a.LStructure(Period=Period,Gamma=Gamma,D_max=D_max,D_min=D_min,C=C,GammaTop=GammaTop,TopLayerThick=TopLayerThick, 
            EachRepetitionLayerThick=EachRepetitionLayerThick, Z_Array=Z_Array)
    sigma=a.roughnessarray(SigmaValues)
    theta_ind = range(bin_num)
    if NumCore != None:
        pool = mulp.Pool(NumCore)
        R = pool.map(theta_scan, theta_ind) # R-> will be a list[Rs,Rp,Ra]
        pool.close()
        flux[:]=np.array(R)[:,2] # Average component of Refl.
    else:
        for i in theta_ind:
            Rs,Rp,Ra = a.get_reflectivity(z,NK_Array,sigma,Theta_in[i],Energy,IncidentBeamDia=IncidentBeamDia)
            flux[i] = Ra


def darpanx_En(engs, params, flux):
    global Energy_in,Period,Gamma,D_max,D_min,C,GammaTop,TopLayerThick,EachRepetitionLayerThick, Z_Array, z, sigma
    bin_num = len(engs)-1
    Energy_in=kev2a(np.array(engs))
    if MultilayerType == 'UserDefinedML':
        TopLayerThick = params[0:TopLayerNum]
        indp=TopLayerNum+EachRepetitionLayerNum
        EachRepetitionLayerThick = params[TopLayerNum:indp]
        if DensityCorrection == 'yes':
            indp2=indp+len(LayerMaterial)
            Density = params[indp:indp2]
            #NK_Array=a.refractiveindexarray(ncsl,Density=Density)
            if ProjEffCorr == 'yes':
                IncidentBeamDia=params[indp2]
                SigmaValues=params[indp2+1:-1]
            else: SigmaValues=params[indp2:-1]
        elif ProjEffCorr == 'yes':
            IncidentBeamDia=params[indp]
            SigmaValues=params[indp+1:-1]
        else: SigmaValues=params[indp:-1]

    elif MultilayerType == 'SingleLayer':
        Period = params[0]
        if DensityCorrection == 'yes':
            Density = params[1]
            #NK_Array=a.refractiveindexarray(ncsl,Density=Density)
            if ProjEffCorr == 'yes':
                IncidentBeamDia=params[2]
                SigmaValues=params[3:-1]
            else: SigmaValues=params[2:-1]
        elif ProjEffCorr == 'yes':
            IncidentBeamDia=params[1]
            SigmaValues=params[2:-1]
        else: SigmaValues=params[1:-1]

    elif MultilayerType == 'BiLayer':
        Period = params[0]
        Gamma = params[1]
        if DensityCorrection == 'yes':
            Density = params[2:4]
            #NK_Array=a.refractiveindexarray(ncsl,LayerDensity=Density)
            if ProjEffCorr == 'yes':
                IncidentBeamDia=params[4]
                SigmaValues=params[5:-1]
            else:SigmaValues=params[4:-1]
        elif ProjEffCorr == 'yes':
            IncidentBeamDia=params[2]
            SigmaValues=params[3:-1]
        else: SigmaValues=params[2:-1]

    elif MultilayerType == 'ClusterGraded':
        indp=NumStack*2
        Period = params[0:indp:2]
        Gamma = params[1:indp:2]
        if DensityCorrection == 'yes':
            indp2=indp+2
            Density = params[indp:indp2]
            #NK_Array=a.refractiveindexarray(ncsl,Density=Density)
            if ProjEffCorr == 'yes':
                IncidentBeamDia=params[indp2]
                SigmaValues=params[indp2+1:-1]
            else:SigmaValues=params[indp2:-1]
        elif ProjEffCorr == 'yes':
            IncidentBeamDia=params[indp]
            SigmaValues=params[indp+1:-1]
        else: SigmaValues=params[indp:-1]

    elif MultilayerType == 'DepthGraded':
        D_max = params[0]
        D_min = params[1]
        GammaTop = params[2]
        Gamma = params[3]
        C = params[4]
        if DensityCorrection == 'yes':
            indp=4+len(LayerMaterial)+1
            Density = params[5:indp]
            #NK_Array=a.refractiveindexarray(ncsl,Density=Density)
            if ProjEffCorr == 'yes':
                IncidentBeamDia=params[indp]
                SigmaValues=params[indp+1:-1]
            else:SigmaValues=params[indp:-1]
        elif ProjEffCorr == 'yes':
            IncidentBeamDia=params[5]
            SigmaValues=params[6:-1]
        else: SigmaValues=params[5:-1]

    elif MultilayerType == 'UserDefined':
        Z_Array = params[0:LayerNum]
        indp=2*LayerNum
        if DensityCorrection == 'yes':
            Density = params[LayerNum:indp]
            #NK_Array=a.refractiveindexarray(ncsl,Density=Density)
            indp2=indp+LayerNum
            if ProjEffCorr == 'yes':
                IncidentBeamDia=params[indp2]
                SigmaValues=params[indp2+1:-1]
            else:SigmaValues=params[indp2:-1]
        elif ProjEffCorr == 'yes':
            IncidentBeamDia=params[indp]
            SigmaValues=params[indp+1:-1]
        else: SigmaValues=params[indp:-1]

    #if DensityCorrection is 'yes': NK_Array=a.refractiveindexarray(ncsl,Density=Density)

    z=a.LStructure(Period=Period,Gamma=Gamma,D_max=D_max,D_min=D_min,C=C,GammaTop=GammaTop,TopLayerThick=TopLayerThick,
            EachRepetitionLayerThick=EachRepetitionLayerThick, Z_Array=Z_Array)
    sigma=a.roughnessarray(SigmaValues)
    energy_ind = range(bin_num)
    if NumCore != None:
        pool = mulp.Pool(NumCore)
        R = pool.map(energy_scan, energy_ind) # R-> will be a list[Rs,Rp,Ra]
        pool.close()
        flux[:]=np.array(R)[:,2] # Average component of Refl.
    else:
        for i in energy_ind:
            Rs,Rp,Ra = a.get_reflectivity(z,NK_Array,sigma,Theta,Energy_in,IncidentBeamDia=IncidentBeamDia)
            flux[i] = Ra

def parinfo():
    sigma_parInfo=[]
    for i in range(len(SigmaValues)):
        sig="sigma"+str(i)+" Angs "+str(SigmaValues[i])+" 0.0 0.0 "+str(SigmaValues[i]+5)+" "+str(SigmaValues[i]+10)+" 1.0"
        sigma_parInfo = sigma_parInfo + [sig]
    sigma_parInfo=tuple(sigma_parInfo)

    density_parInfo=[]
    if DensityCorrection == 'yes': 
        for i in range(len(LayerMaterial)):
            den="rho_"+LayerMaterial[i]+" g/cm3 "+str(Density[i])+" 0.0 0.0 "+str(0.1)+" "+str(Density[i]+1.0)+" 1.0"
            density_parInfo=density_parInfo + [den]
    density_parInfo=tuple(density_parInfo)

    if ProjEffCorr == 'yes':
        beam_dia_parInfo=["BeamDia mm "+str(IncidentBeamDia)+" 0.0001 0.0001 "+
                str(IncidentBeamDia+0.2)+" "+str(IncidentBeamDia+0.3)+ " 0.01"]
    else:
        beam_dia_parInfo=[]
    beam_dia_parInfo=tuple(beam_dia_parInfo)

    if MultilayerType == 'UserDefinedML':
        thickness_parInfo=[]
        for i in range(TopLayerNum):
            d="d_"+LayerMaterial[i]+" Angs "+str(TopLayerThick[i])+" 0.0 0.0 "+\
            str(TopLayerThick[i]+50)+" "+str(TopLayerThick[i]+80)+" 1.0"
            thickness_parInfo=thickness_parInfo+[d]
        for i in range(TopLayerNum,TopLayerNum+EachRepetitionLayerNum):
            ind=i-TopLayerNum
            d="d_"+LayerMaterial[i]+" Angs "+str(EachRepetitionLayerThick[ind])+" 0.0 0.0 "+\
            str(EachRepetitionLayerThick[ind]+50)+" "+str(EachRepetitionLayerThick[ind]+80)+" 1.0"
            thickness_parInfo=thickness_parInfo+[d]
        thickness_parInfo=tuple(thickness_parInfo)

    elif MultilayerType == 'SingleLayer':
        thickness_parInfo=["d_"+LayerMaterial[0]+" Angs "+str(Period[0])+" 0.0 0.0 "+str(Period[0]+50)+" "+str(Period[0]+80)+" 1.0"]
        thickness_parInfo=tuple(thickness_parInfo)

    elif MultilayerType == 'BiLayer':
        thickness_parInfo=["Period Angs "+str(Period[0])+" 0.0 0.0 "+str(Period[0]+50.)+" "+str(Period[0]+80.)+" 1.0",
                           "Gamma \"\" "+str(Gamma[0])+" 0.0 0.1 0.4 1.0 1.0"]
        thickness_parInfo=tuple(thickness_parInfo)

    elif MultilayerType == 'DepthGraded':
        thickness_parInfo=["D_max Angs "+str(D_max)+" 0.0 0.0 "+str(D_max+50)+" "+str(D_max+80)+" 1.0",
                           "D_min Angs "+str(D_min)+" 0.0 0.0 "+str(D_min+10)+" "+str(D_min+15)+" 1.0",
                           "GammaTop \"\" "+str(GammaTop)+" 0.0 0.1 0.4 1.0 0.1",
                           "Gamma \"\" "+str(Gamma[0])+" 0.0 0.1 0.4 1.0 0.1",
                           "C(slope) \"\" "+str(C)+" 0.0 0.0 0.4 1.0 0.01"]
        thickness_parInfo=tuple(thickness_parInfo)
    
    elif MultilayerType == 'ClusterGraded':
        thickness_parInfo=[]
        #try:
        for i in range(NumStack):
            d=["Period_St"+str(i)+" Angs "+str(Period[i])+" 0.0 0.0 "+str(Period[i]+10)+" "+str(Period[i]+20)+" 1.0",
                "Gamma_St"+str(i)+" \"\" "+str(Gamma[i])+" 0.0 0.0 0.4 1.0 0.1"]
            thickness_parInfo=thickness_parInfo+d
        thickness_parInfo=tuple(thickness_parInfo)
        #except: raise Exception("%% DarpanX_Error: Incorrect dimension of < period > and < gamma >. Dimension should be- "+str(NumStack))
    elif MultilayerType == 'UserDefined':
        thickness_parInfo=[]
        for i in range(LayerNum):
            d="d_"+LayerMaterial[i]+" Angs "+str(Z_Array[i])+" 0.0 0.0 "+str(Z_Array[i]+50)+" "+str(Z_Array[i]+80)+" 1.0"
            thickness_parInfo=thickness_parInfo+[d]
        thickness_parInfo=tuple(thickness_parInfo)

    else:
        raise Exception("%% DarpanX_Error: Undefined or in correct < MultilayerType >")

    parInfo = thickness_parInfo + density_parInfo + beam_dia_parInfo + sigma_parInfo
    return parInfo

parInfo=parinfo()


def addmodel():
    if Xscan =='Theta':
        AllModels.addPyMod(darpanx_Th, parInfo, 'add')
    elif Xscan =='Energy': AllModels.addPyMod(darpanx_En, parInfo, 'add')

addmodel()

