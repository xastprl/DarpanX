from xspec import*
import darpanx as drp
UserInFile='Example2_XRR_SingleSiModel.drpnx'
drp.call_wrap(UserInFile)
Xset.restore('Example2_XRR_SingleSiModel.xcm')
drp.iplot(comp="delc")
