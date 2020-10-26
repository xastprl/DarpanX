import matplotlib.pyplot as plt
import numpy as np

file='Example6_XRR_W_B4C_ClusterGraded_ModelFit'
data_dir='./'
num_block=5

file_name=file+'.xcm'
par,delta,hlo,lo,h,hh=np.loadtxt(data_dir+file_name, skiprows=18,unpack=True)

g=par[3:(2*num_block)+2:2]

d=par[2:(2*num_block)+2:2]

d_w=[1]*num_block
d_b4c=[1]*num_block

for i in range(num_block): d_w[i]=d[i]*g[i] ; d_b4c[i]=d[i]*(1.0-g[i])

import matplotlib.pyplot as plt

fig, ax = plt.subplots(num=None, figsize=(13, 8), dpi=80, facecolor='w', edgecolor='k')
ax2=ax.twinx() # twin object for two different y-axis on the sample plot
ax.set_xlim(0.5,num_block+0.5)
ax2.set_xlim(0.5,num_block+0.5)
ax.set_xlabel("Block number",fontsize=21)
ax.set_ylabel("Period ($\AA$)",fontsize=22)
ax2.set_ylabel("Gamma",fontsize=24)
ax.set_xticks(np.arange(1, num_block+1,1 ))
ax2.set_xticks(np.arange(1, num_block+1,1 ))
ax.tick_params(labelsize=20)
ax2.tick_params(labelsize=20)
lns1 = ax.plot(range(1,num_block+1),d, 'g^' ,markersize=11,linewidth=3,label='Period')
lns2 = ax2.plot(range(1,num_block+1), g,'m*',markersize=11,linewidth=3,label='Gamma')
lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc='lower left',fontsize='x-large',prop={'size':22},labelspacing=0.3)
#ax.grid()
#ax.axhline(sum(d)/num_block, color='g', lw=2, alpha=0.5)
#ax2.axhline(sum(g)/num_block, color='m', lw=2, alpha=0.5)
plt.savefig(data_dir+'Example6_XRR_Cluster_Graded_FitPar_plot.pdf')
plt.show()
