from xspec import*
import darpanx as drp
UserInFile='Example6_XRR_W_B4C_BiLayer_ModelFit.drpnx'
drp.call_wrap(UserInFile)
Xset.restore('Example6_XRR_W_B4C_BiLayer_ModelFit.xcm')
drp.iplot()
