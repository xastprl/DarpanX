from xspec import*
import darpanx as drp
UserInFile='Example6_XRR_W_B4C_ClusterGraded_ModelFit.drpnx'
drp.call_wrap(UserInFile)
Xset.restore('Example6_XRR_W_B4C_ClusterGraded_ModelFit.xcm')
drp.iplot()
