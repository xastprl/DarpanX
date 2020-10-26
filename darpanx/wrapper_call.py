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

from darpanx.get_dir import*
import sys
def call_wrap(UserInFile): 
    wrap=get_wrap()
    sys.argv = UserInFile
    exec(open(wrap).read(),globals(), globals())


def iplot(device='/xw',xlog=False,ylog=True,xaxis='kev',comp=' ',OutFile=None,xlabel='Theta (degree)',ylabel='counts/s/degree'):
    if OutFile == None:
        Plot.device = device
    else: Plot.device = OutFile
    Plot.xAxis=xaxis
    Plot.xLog=xlog
    Plot.yLog=ylog
    #if Xscan == 'Theta':
    if comp != ' ':
        Plot.addCommand("window 2") 
    Plot.addCommand("label x "+xlabel)
    Plot.addCommand("window 1")
    Plot.addCommand("label y normalized "+ylabel)
    #Plot("data "+comp)
    try:
        Plot("data "+comp)
    except:
        print("%% DarpanX_message: No data loaded... Ploting model component.")
        Plot("mo")
    if OutFile != None: Plot.device = 'close'

def reiplot(comp=' ',xlabel='Theta (degree)',ylabel='counts/s/degree'):
    if comp != ' ':
        Plot.addCommand("window 2")
    Plot.addCommand("label x "+xlabel)
    Plot.addCommand("window 2")
    Plot.addCommand("label y normalized "+ylabel)
    Plot("data "+comp)

def plot(xlim=None,ylim=None,OutFile=None,OutFileFormat='pdf',Struc='no',Scale ='no',title='DarpanX Output',ylog='yes',xlog='no'):
    Z_Array=z
    Plot.xAxis='ch'
    Plot.xLog=False
    Plot.yLog=True
    Plot("data")
    #try:
    #    Plot("data")
    #except:
    #    print("%% DarpanX_message: No data loaded... Ploting model component.")
    #    Plot("mo")
    refl=np.array(Plot.y())
    yErrs = np.array(Plot.yErr())
    systErr=float(AllModels.systematic)
    yErrs=np.sqrt(yErrs**2 + (refl*systErr)**2)
    refl_model=Plot.model()
    Plot.xAxis='kev'
    Plot("data")
    xErrs = np.array(Plot.xErr())
    theta=Plot.x()
    fig, ax = plt.subplots(num=None, figsize=(10, 8), dpi=80, facecolor='w', edgecolor='k')
    plt.title(title,fontsize=18,color='k',loc='left')
    if xlog == 'yes':plt.xscale('log')
    if ylog == 'yes':plt.yscale('log')
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    if xlim != None:ax.set_xlim(xlim[0],xlim[1])
    if ylim != None:ax.set_ylim(ylim[0],ylim[1])
    plt.ylabel('Reflectivity', fontsize=24)
    plt.xlabel('Grazing incidence angle (degree)', fontsize=24)
    #plt.plot(theta, refl, 'r+',markersize=10,linewidth=2,label='Data')
    plt.errorbar(theta,refl , yerr=yErrs,xerr=xErrs, marker='p', ms=1, mew=1,barsabove=True,elinewidth=1,capsize=1,ecolor='k',color='red',label='Data')
    plt.plot(theta, refl_model,'b--',markersize=10,linewidth=1.2,label='Model')
    #plt.step(theta,refl_model,'b-',linewidth=1.,label='Model', where='mid')
    plt.legend(fontsize=22,loc='upper right')
    if OutFile != None:
        plt.savefig(OutFile+'.'+OutFileFormat)
    def substrate_fill():
        totalthick=(3.0*sum(Z_Array)/100.0)#1.8
        plt.fill_between([0,1.1], y1=-totalthick,y2=0.1, color='k',alpha=1.0)
        plt.annotate(SubstrateMaterial, (0.50,-totalthick),color='white',fontsize=15.0) 
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
        if MultilayerType == 'UserDefinedML':
            substrate_fill()
            y1=0.1
            for i in range(Repetition):
                y2=y1
                for j in range(EachRepetitionLayerNum-1,-1,-1):
                    y2=y2+EachRepetitionLayerThick[j]
                    if i == 0:plt.fill_between([0,1.1], y1=y1,y2=y2, color=cl[j],alpha=1.0,label=LayerMaterial[TopLayerNum+j])
                    else: plt.fill_between([0,1.1], y1=y1,y2=y2, color=cl[j],alpha=1.0)
                    y1=y2
            if TopLayerNum > 0:
                for k in range(TopLayerNum-1,-1,-1):
                    y2=y2+TopLayerThick[k]
                    plt.fill_between([0,1.1], y1=y1,y2=y2, color=cl1[k],alpha=1.0,label=LayerMaterial[k])
                    y1=y2
            ax1.legend(loc='lower left', bbox_to_anchor= (0.02, -0.02), ncol=len(LayerMaterial),borderaxespad=0, frameon=False,fontsize=18)
        elif MultilayerType == 'BiLayer' or MultilayerType == 'DepthGraded' or MultilayerType == 'ClusterGraded' or MultilayerType == 'SingleLayer':
            substrate_fill()          
            y1=0.1
            for i in range(len(Z_Array)-1,-1,-1):
                if (i % 2.0 == 0):
                    clbil='red'
                    mat=LayerMaterial[0]
                else:
                    clbil='blue'
                    mat=LayerMaterial[1]
                y2=y1+Z_Array[i]
                if i > 1:plt.fill_between([0,1.1], y1=y1,y2=y2, color=clbil,alpha=1.0)
                elif i == 1: plt.fill_between([0,1.1], y1=y1,y2=y2, color=clbil,alpha=1.0,label=mat)
                elif i ==0: plt.fill_between([0,1.1], y1=y1,y2=y2, color=clbil,alpha=1.0,label=mat)
                y1=y2
            ax1.legend(loc='lower left', bbox_to_anchor= (0.02, -0.02), ncol=len(LayerMaterial),borderaxespad=0, frameon=False,fontsize=18)
       
        elif MultilayerType == 'UserDefined':
            substrate_fill()
            y1=0.1
            lm=np.array(LayerMaterial)
            dup_mat=find_duplicate(LayerMaterial)
            cl_usd=lm.astype('U256')

            for i in range(len(dup_mat)):
                c=np.where(lm == dup_mat[i])[0]
                cl_usd[c]=cl1[i]
                c35=0
                for k in range(len(c)):
                    c35=c35+1
                    if c[k] == 0:
                        y1=0.1
                        y2=y1+Z_Array[c[k]]
                        #print(c[k],y1,y2)
                    else:
                        y1=0.1+sum(Z_Array[0:c[k]])
                        y2=y1+Z_Array[c[k]]
                        #print(c[k],y1,y2)
                    if c35 == 1:
                        plt.fill_between([0,1.1], y1=y1,y2=y2, color=cl_usd[c[k]],alpha=1.0,label=LayerMaterial[c[k]])
                    else:plt.fill_between([0,1.1], y1=y1,y2=y2, color=cl_usd[c[k]],alpha=1.0)

            ax1.legend(loc='lower left', bbox_to_anchor= (0.02, -0.02), ncol=len(dup_mat),borderaxespad=0, frameon=False,fontsize=18)
        if Scale == 'yes':
            TZ=sum(Z_Array)
            TZ1="%0.2f"%TZ
            plt.annotate("", xy=(1.0, 0.1), xytext=(1.0, 0.1+TZ),arrowprops=dict(arrowstyle="<->", connectionstyle="arc3",color='white',lw=2.0))
            plt.text(1.0, (0.1+TZ/2.0), str(TZ1)+'$\AA$',{'color': 'black', 'fontsize': 24, 'ha': 'center', 'va': 'center','bbox': dict(boxstyle="round", fc="white", ec="black", pad=0.2)},rotation=90)
            TZ=0
            TZ1=0

        '''
        elif MultilayerType == 'UserDefined':
            substrate_fill()
            y1=0.1
            for k in range(LayerNum):
                y2=y1+Z_Array[k]
                plt.fill_between([0,1.1], y1=y1,y2=y2, color=cl1[k],alpha=1.0,label=LayerMaterial[k])
                y1=y2
        '''
        if OutFile is not None:
            plt.savefig(OutFile+'_structure.'+OutFileFormat)
    plt.show(block=False)

def save(OutFile='test_darpanx',SaveAll='yes'):
    Xset.save(OutFile+'.xcm', info='a')
    current_time = datetime.datetime.now()
    with open(UserInFile) as f:
        with open(OutFile+".drpnx", "w") as f1:
            f1.write("# This file is generted by DarpanX on: "+str(os.uname()[1])+" at "+str(current_time.year)+"-"+str(current_time.month)+"-"+str(current_time.day)+". Don't edit it:")
            f1.write('\n')
            f1.write("# ======== Inputs ==========")
            f1.write('\n')
            #for line in f:
            #    #if "ROW" in line:
            #    if line.strip():
            #        f1.write(line)
            #    if '# ======== Outputs ==========' in line: break
            for i in inp:
                f1.write(i)
            #if SaveAll == 'yes':
            #s=AllData(1)
            #dat=s.values
            m=AllModels(1)
            mo_refl = m.folded(1)
            #Plot.xAxis='ch'
            #Plot.xLog=False
            #Plot.yLog=True
            #Plot("data")
            #refl=Plot.y()
            #mo_refl=Plot.model()
            Plot.xAxis='kev'
            Plot("data")
            theta=Plot.x()
            #parms_name=AllModels(1).parameterNames
            f1.write('\n')
            f1.write("# ======== Outputs ==========")
            f1.write('\n')
            n=AllModels(1).nParameters
            chi=(Fit.statistic)
            dof=Fit.dof
            reduce_chi=chi/dof
            f1.write("# Parameter name = [value  fit-delta  min  bot  top  max]")
            f1.write('\n')
            for i in range(n):
                f1.write(str(AllModels(1)(i+1).name)+" = "+str(AllModels(1)(i+1).values))
                f1.write('\n')
            try:
                f1.write('\n')
                f1.write("offset"," = "+str(s.response.gain.offset.values))
                f1.write("gain"," = "+str(s.response.gain.slope.values))
                f1.write('\n')
            except: pass
            f1.write('%s\t'%"CHI = "+str(chi))
            f1.write('\n')
            f1.write('%s\t'%"DOF = "+str(dof))
            f1.write('\n')
            f1.write('%s\t'%"Red-Chi = "+str(reduce_chi))
            f1.write('\n')
            f1.write('\n')
            if SaveAll == 'yes':
                f1.write("# ======== Data, Model, Out ==========")
                f1.write('\n')
                if Xscan == 'Theta':
                    #f1.write("# Theta (deg), data-reflectivity, model-reflectivity")
                    f1.write("# Theta (deg), model-reflectivity")
                    f1.write('\n')
                    #TH=(Theta_in[0:-1:2]+Theta_in[1::2])/2.0 #cosider mid point of the xspec binning
                    TH=theta
                    for i in range(len(TH)):
                        #f1.write('%0.10f\t' % rad2deg(TH[i]))
                        f1.write('%0.10f\t' % (TH[i]))
                        #f1.write('%0.10f\t' % dat[i])
                        f1.write('%0.10f\t' % mo_refl[i])
                        f1.write('\n')
                elif Xscan == 'Energy':
                    #f1.write("# Energy (keV), data-reflectivity, model-reflectivity")
                    f1.write("# Energy (keV), model-reflectivity")
                    f1.write('\n')
                    #EN=(Energy_in[0:-1:2]+Energy_in[1::2])/2.0
                    EN=theta
                    for i in range(len(EN)):
                        #f1.write('%0.10f\t' % a2kev(Energy_in[i]))
                        f1.write('%0.10f\t' % EN[i])
                        #f1.write('%0.10f\t' % dat[i])
                        f1.write('%0.10f\t' % mo_refl[i])
                        f1.write('\n')
    with open(OutFile+'.py',"w") as f0:
        f0.write("from xspec import*")
        f0.write('\n')
        f0.write("import darpanx as drp")
        f0.write('\n')
        f0.write("UserInFile="+"'"+OutFile+".drpnx"+"'")
        f0.write('\n')
        #f0.write("wrap=drp.get_wrap()")
        f0.write("drp.call_wrap(UserInFile)")
        f0.write('\n')
        #f0.write("exec(open(wrap).read(),globals())")
        #pv=int(sys.version[0])
        #if pv >= 3:f0.write("exec(open(wrap).read())")
        #elif pv < 3: f0.write("execfile(wrap)")
        #else: print("%% DarpanX_Error: Not getting python version")
        #f0.write('\n')
        f0.write("Xset.restore("+"'"+OutFile+".xcm"+"'"+")")
        f0.write('\n')
        f0.write("drp.iplot()")
        f0.close()
def make_model():
    if Xscan == 'Theta':
        if InsRes > 0:
            print("Initial Parameters")
            m=Model("gsmooth*darpanx_Th")
            AllModels(1)(1).values=[InsRes, -0.00001, 0.0, 0.0, InsRes, InsRes+0.01]
        elif InsRes == 0: m=Model("darpanx_Th")
        if ProjEffCorr == 'yes': m.darpanx_Th.BeamDia.frozen=True
        m.darpanx_Th.norm.values=[1.0, 0.01, 0.0, 0.5, 1.0, 2.0]
        AllModels.show()
    elif Xscan == 'Energy':
        if InsRes > 0:
            m=Model("gsmooth*darpanx_En")
            AllModels(1)(1).values=[InsRes, -0.00001, 0.0, 0.0, InsRes, InsRes+0.01]
        elif InsRes == 0: m=Model("darpanx_En")
        if ProjEffCorr == 'yes':m.darpanx_En.BeamDia.frozen=True
        m.darpanx_En.norm.values=[1.0, 0.01, 0.0, 0.5, 1.0, 2.0]
        AllModels.show()

def gainfit(status='on'):
    s=AllData(1)
    s.response.gain.slope = 1.0
    s.response.gain.offset.values=[0.0, 0.01, -0.1, -0.1, 0.1, 0.1]
    s.response.gain.slope.frozen=True
    if status == 'off':s.response.gain.off()


