import swifttools as st
import numpy as np

obsid='00010627113'
tsrcra = '275.09134'
tsrcdec = '7.1853'
obsmode='Pointing'
stage2chatter = 0
stage3chatter = 0
exprwtgrade = '0'
exprpcgrade = '0-12'
datadir = '/Volumes/HighEnergy1/archive/Swift/J1820/'
caldbenv = '/Users/tdincer/Softwares/CALDB'
caldbpath = caldbenv + '/data/swift/xrt/cpf/rmf'

wtsrcregionfile = ''
pileup = 1
pcsrcregionfile = ''
pcsnapshot = 0



srcinnerradius = '0"'
srcouterradius = '70.8"'

# Here is the steps
#st.stage2(obsid, tsrcra, tsrcdec, exprpcgrade, exprwtgrade, obsmode, stage2chatter, datadir)
#st.snapshots(obsid, tsrcra, tsrcdec, exprpcgrade, exprwtgrade, obsmode, datadir, caldbpath, stage3chatter, pcsnapshot, wtsrcregionfile, pcsrcregionfile, pileup, srcinnerradius, srcouterradius)
#grppha sw00010627113xwtw2posr.pha sw00010627113xwtw2posr.grp comm="group min 1 & chkey RESPFILE swxwt0s6psf1_20131212v001.rmf & chkey ANCRFILE sw00010627113xwtw2posr.arf & exit"