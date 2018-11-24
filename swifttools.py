import os
import sys
import glob
import numpy as np
from numpy import sqrt, indices, bincount, arange, pi
import subprocess
from astropy.io import fits
from astropy import units
from astropy import constants as const


def stage2(obsid, tsrcra, tsrcdec, exprpcgrade, exprwtgrade, obsmode, stage2chatter, datadir):
    """
    This function extracts the Stage 2 data products by applying default screening criteria to the raw data.
    The default criteria can be seen in XRT Swift Guide.
    :param obsid: observation ID.
    :param tsrcra: right ascension of the source
    :param tsrcdec: declination of the source
    :param exprpcgrade: pc grade
    :param exprwtgrade: wt grade
    :param obsmode: observering mode
    :param stage2chatter: chatter degree
    :param datadir: path to the input data directory
    How the values of these parameters should be entered is explained in an example below.

    Folder Structure:
    This function creates a folder with the naming convention obsid_res (e.g. 00010264028_res). The Stage 2 products are
    stored in the '00010264028_res/Stage2' folder.

    Usage:
    obsid = '00010264028'
    tsrcra = '233.832250'
    tsrcdec = '-57.229667'
    obsmode = 'Pointing'
    stage2chatter = 0
    exprwtgrade = '0'
    exprpcgrade = '0-12'
    datadir = '/Volumes/HighEnergy1/archive/Swift/MAXIJ1535/'
    stage2(obsid, tsrcra, tsrcdec, exprpcgrade, exprwtgrade, obsmode, stage2chatter, datadir)
    """

    print('Producing Stage2 for ObsID:'), obsid  # Analyzing 000701476
    resdir = obsid + '_res'

    os.mkdir(resdir)
    os.chdir(resdir)

    os.system('punlearn xrtpipeline')
    xrtpipeline_string = 'xrtpipeline indir=%s%s outdir=Stage2 steminputs=sw%s srcra=%s srcdec=%s exprpcgrade=%s ' \
                         'exprwtgrade=%s createexpomap=yes cleancols=yes chatter=%s cleanup=yes ' \
                         'exitstage=2 obsmode=%s' \
                         % (datadir, obsid, obsid, str(tsrcra), str(tsrcdec), exprpcgrade, exprwtgrade,
                            str(stage2chatter), obsmode)
    print(xrtpipeline_string)
    p = subprocess.Popen(xrtpipeline_string, shell=True)
    p.communicate()
    
    
def snapshots(obsid, tsrcra, tsrcdec, exprpcgrade, exprwtgrade, obsmode, datadir, caldbpath, stage3chatter, pcsnapshot,
              wtsrcregionfile, pcsrcregionfile, pileup, srcinnerradius, srcouterradius):
    """
    :param obsid: Observation ID (e.g. '00010264028').
    :param tsrcra: right ascension of the source (e.g. '233.832250').
    :param tsrcdec: declination of the source (e.g. '-57.229667').
    :param exprpcgrade: pc grade.
    :param exprwtgrade: wt grade.
    :param obsmode: observing mode (e.g. 'Pointing').
    :param datadir: path to the input data directory.
    :param caldbpath: path to the calibration folder (e.g. '/Users/tdincer/Softwares/CALDB').
    :param stage3chatter: chatter degree.
    :param pcsnapshot:
    :param wtsrcregionfile: source region file for wt mode data. If this parameter is empty (e.g ''), the code will use
           a default region.
    :param pcsrcregionfile:
    :param pileup: turn on/off pile-up thread (0 for on, 1 for off).
                   it helps to choose from one of predefined regions, namely circle and annulus.
    :param srcinnerradius: inner radius for the annulus source region.
    :param srcouterradius: outer radius for the annulus source region.
    """

    # Create the folders
    basedir = os.getcwd()
    resdir = obsid + '_res'
    os.chdir(resdir)
        
    try:
        os.mkdir('WT')
    except OSError:
        print "Can't create WT directory, it already exists!"
        pass
    try:
        os.mkdir('WT/GTI')
    except OSError:
        print "Can't create WT directory, it already exists!"
        pass
    try:
        os.mkdir('PC')
    except OSError:
        print "Can't create PC directory, it already exists!"
        pass
    try:
        os.mkdir('PC/GTI')
    except OSError:
        print "Can't create WT directory, it already exists!"
        pass
   
    # Find the level2 event files in pointing mode.
    wtevtfile = glob.glob('./Stage2/*wt*po_cl.evt')
    pcevtfile = glob.glob('./Stage2/*pc*po_cl.evt')

    if (wtevtfile == []) and (pcevtfile == []):
        print '@@@ There is neither WT nor PC mode pointing event file!'
        sys.exit()
    elif (pcevtfile == []) and (wtevtfile != []):
        print '@@@ There is only WT mode pointing event file!'
        print 'WT event file: ', wtevtfile
    elif (pcevtfile != []) and (wtevtfile == []):
        print '@@@ There is only PC mode pointing event file!'
        print 'PC event file: ', pcevtfile
    elif (pcevtfile != []) and (wtevtfile != []):
        print "@@@ There is both WT and PC mode pointing data!"
        print 'WT event file: ', wtevtfile
        print 'PC event file: ', pcevtfile
        
    if wtevtfile != []:
        wtgtifile = glob.glob('./Stage2/*wt*po_clgti.fits')[0]
        hdufile = fits.open(wtgtifile)
        hdr = hdufile[1].header
        nrowswt = hdr['NAXIS2']
        print 'Number of rows in GTI: ', nrowswt

        for ii in np.arange(1, nrowswt+1):
            if hdufile[1].data[ii-1][1]-hdufile[1].data[ii-1][0] > 400:
                fselect_expr = 'fselect %s+1 %swt_gti%s.fits "#row==%s"' % (wtgtifile, obsid, ii, ii)
                print fselect_expr
                p = subprocess.Popen(fselect_expr, shell=True)
                p.communicate()
                mv_expr = 'mv %swt_gti%s.fits WT/GTI' % (obsid, ii)
                os.system(mv_expr)
            else:
                print "GTI is < 400s"

        nrowsforloop = nrowswt + 1
        for ii in np.arange(1, nrowsforloop):
            if hdufile[1].data[ii - 1][1] - hdufile[1].data[ii - 1][0] > 400:
                try:
                    os.mkdir('WT/Snapshot' +str(ii))
                except OSError:
                    print "Can't create WT/Snapshot%s directory, it already exists!" %(str(ii))
                    pass
            
     
    #if pcevtfile != []:
    #    pcgtifile = glob.glob('./Stage2/*pc*po_clgti.fits')[0]
    #    hdufile = fits.open(pcgtifile)
    #    hdr = hdufile[1].header
    #    nrowspc = hdr['NAXIS2']
    #    print 'Number of rows in GTI: ', nrowspc

    #    for ii in np.arange(1, nrowspc+1):
    #        fselect_expr = 'fselect %s+1 %spc_gti%s.fits "#row==%s"' % (pcgtifile, obsid, ii, ii)
    #        print fselect_expr
    #        p = subprocess.Popen(fselect_expr, shell=True)
    #        p.communicate()
    #        mv_expr = 'mv %spc_gti%s.fits PC/GTI' % (obsid, ii)
    #        os.system(mv_expr)
            
    #    nrowsforloop = nrowspc + 1
    #    for ii in np.arange(1, nrowsforloop):
    #        os.mkdir('PC/Snapshot' +str(ii))
            
    # Find the attitudefile
    attfile = glob.glob(datadir+obsid+'/auxil/*pat.fits*')
    # Find the HD file
    hdfile = glob.glob('./Stage2/*xhdtc.hk')
    print("-----------------         Input Files        ----------------")
    print("Attitude File Name : "), attfile
    print("HD File Name       : "), hdfile
    
    if wtevtfile != []:
        os.chdir('WT')
        for snapshotfolder in glob.glob('Snapshot*'):
            print snapshotfolder
            sn = snapshotfolder[-1]  # This will be an integer.
            os.chdir(str(snapshotfolder))

            if exprwtgrade == '0-2':
                os.system('cp ' + caldbpath + '/swxwt0to2s6psf1_20131212v001.rmf .')
                os.system('cp ' + caldbpath + '/swxwt0to2s6psf2_20131212v001.rmf .')
                os.system('cp ' + caldbpath + '/swxwt0to2s6psf3_20131212v001.rmf .')
                defwtrmf = 'swxwt0to2s6psf1_20131212v001.rmf'
            elif exprwtgrade == '0':
                os.system('cp ' + caldbpath + '/swxwt0s6psf1_20131212v001.rmf .')
                os.system('cp ' + caldbpath + '/swxwt0s6psf2_20131212v001.rmf .')
                os.system('cp ' + caldbpath + '/swxwt0s6psf3_20131212v001.rmf .')
                defwtrmf = 'swxwt0s6psf1_20131212v001.rmf'

            if wtsrcregionfile == '':
                if pileup == 0:
                    srcregexp = "fk5; circle(%s , %s, %s)" % (tsrcra, tsrcdec, srcouterradius)
                elif pileup == 1:
                    srcregexp = "fk5; annulus(%s , %s, %s, %s)" % (tsrcra, tsrcdec, srcinnerradius, srcouterradius)

                srcregion = open('src.reg', 'w')
                srcregion.write(srcregexp)
                srcregion.close()
            else:
                pass

            bkgregexp = 'fk5; annulus(' + tsrcra + ',' + tsrcdec + ', 165", 200")'

            bkgregion = open('bkg.reg', 'w')
            bkgregion.write(bkgregexp)
            bkgregion.close()

            usrgtifile = '../GTI/%swt_gti%s.fits' % (obsid, sn)
            os.system('punlearn xrtpipeline')
            xrtpip_str = 'xrtpipeline indir=%s%s outdir=./ steminputs=sw%s srcra=%s ' \
                         'srcdec=%s exprwtgrade=%s createexpomap=yes extractproducts=yes ' \
                         'cleancols=yes chatter=%s exitstage=3 obsmode=%s ' \
                         'usrgtifile=%s wtregionfile=src.reg ' \
                         'wtrmffile=%s psfflag=yes datamode=WT clobber=yes wtbinsize=0.00390625' % (datadir, \
                         obsid, obsid, str(tsrcra), str(tsrcdec), exprwtgrade, \
                         str(stage3chatter), obsmode, usrgtifile, defwtrmf)
            print os.getcwd()
            print xrtpip_str              
            p = subprocess.Popen(xrtpip_str, shell=True)
            p.communicate()

            # Create the background spectrum
            bkgspec = open('bkg_spec_creator.xco', 'w')
            bkgspec.write('test\n')
            bkgspec.write('read event\n')
            bkgspec.write('.\n')
            cleanevent = glob.glob('*po_cl.evt')[0]
            bkgspec.write(cleanevent+'\n')
            bkgspec.write('yes\n')
            bkgspec.write('filter region bkg.reg\n')
            bkgspec.write('extract spectrum\n')
            bkgspec.write('yes\n')
            bkgspec.write('save spectrum bkg.pha\n')
            bkgspec.write('quit\n')
            bkgspec.write('yes\n')
            bkgspec.close()

            os.system('xselect @bkg_spec_creator.xco')

            os.chdir('../')
            
    if pcevtfile != []:
        pcdir = basedir+'/'+resdir+'/PC'
        os.chdir(pcdir)

        if exprpcgrade == '0-12':
            os.system('cp ' + caldbpath + '/swxpc0to12s6_20130101v014.rmf .')
            defpcrmf = 'swxpc0to12s6_20130101v014.rmf'
        elif exprwtgrade == '0-4':
            os.system('cp ' + caldbpath + '/swxpc0to4s6_20130101v014.rmf .')
            defpcrmf = 'swxpc0to4s6_20130101v014.rmf'
        elif exprwtgrade == '0':
            os.system('cp ' + caldbpath + '/swxpc0s6_20130101v014.rmf .')
            defpcrmf = 'swxpc0s6_20130101v014.rmf'

        if pcsrcregionfile == '':
            if pileup == 0:
                srcregexp = "fk5; circle(%s , %s, %s)" % (tsrcra, tsrcdec, srcouterradius)
            elif pileup == 1:
                srcregexp = "fk5; annulus(%s , %s, %s, %s)" % (tsrcra, tsrcdec, srcinnerradius, srcouterradius)

        bkgregexp = 'fk5; annulus(' + tsrcra + ',' + tsrcdec + ', 220", 250")'

        srcregion = open('src.reg', 'w')
        srcregion.write(srcregexp)
        srcregion.close()

        bkgregion = open('bkg.reg', 'w')
        bkgregion.write(bkgregexp)
        bkgregion.close()

        os.system('punlearn xrtpipeline')
        xrtpip_str = 'xrtpipeline indir=%s%s outdir=./ steminputs=sw%s srcra=%s ' \
                     'srcdec=%s exprpcgrade=%s createexpomap=yes extractproducts=yes ' \
                     'cleancols=yes chatter=%s exitstage=3 obsmode=%s ' \
                     'pcregionfile=src.reg ' \
                     'pcrmffile=%s psfflag=yes datamode=PC clobber=yes' % (datadir, \
                     obsid, obsid, str(tsrcra), str(tsrcdec), exprpcgrade, \
                     str(stage3chatter), obsmode, defpcrmf)
        print os.getcwd()
        print xrtpip_str
        p = subprocess.Popen(xrtpip_str, shell=True)
        p.communicate()

        # Create the background spectrum
        bkgspec = open('bkg_spec_creator.xco', 'w')
        bkgspec.write('test\n')
        bkgspec.write('read event\n')
        bkgspec.write('.\n')
        cleanevent = glob.glob('*po_cl.evt')[0]
        bkgspec.write(cleanevent+'\n')
        bkgspec.write('yes\n')
        bkgspec.write('filter region bkg.reg\n')
        bkgspec.write('extract spectrum\n')
        bkgspec.write('yes\n')
        bkgspec.write('save spectrum bkg.pha\n')
        bkgspec.write('quit\n')
        bkgspec.write('yes\n')
        bkgspec.close()

        os.system('xselect @bkg_spec_creator.xco')

        os.chdir('../')


def singleshot(obsid, tsrcra, tsrcdec, exprpcgrade, exprwtgrade, obsmode, datadir, caldbpath, stage3chatter, pcsnapshot,
               wtsrcregionfile, pcsrcregionfile, pileup, srcinnerradius, srcouterradius):
    """
    :param obsid: Observation ID (e.g. '00010264028').
    :param tsrcra: right ascension of the source (e.g. '233.832250').
    :param tsrcdec: declination of the source (e.g. '-57.229667').
    :param exprpcgrade: pc grade.
    :param exprwtgrade: wt grade.
    :param obsmode: observing mode (e.g. 'Pointing').
    :param datadir: path to the input data directory.
    :param caldbpath: path to the calibration folder (e.g. '/Users/tdincer/Softwares/CALDB').
    :param stage3chatter: chatter degree.
    :param pcsnapshot:
    :param wtsrcregionfile: source region file for wt mode data. If this parameter is empty (e.g ''), the code will use
           a default region.
    :param pcsrcregionfile:
    :param pileup: turn on/off pile-up thread (0 for on, 1 for off).
                   it helps to choose from one of predefined regions, namely circle and annulus.
    :param srcinnerradius: inner radius for the annulus source region.
    :param srcouterradius: outer radius for the annulus source region.
    """

    # Create the folders
    basedir = os.getcwd()
    resdir = obsid + '_res'
    os.chdir(resdir)

    try:
        os.mkdir('WT')
    except OSError:
        print "Can't create WT directory, it already exists!"
        pass
    try:
        os.mkdir('WT/GTI')
    except OSError:
        print "Can't create WT directory, it already exists!"
        pass
    try:
        os.mkdir('PC')
    except OSError:
        print "Can't create PC directory, it already exists!"
        pass
    try:
        os.mkdir('PC/GTI')
    except OSError:
        print "Can't create WT directory, it already exists!"
        pass

    # Find the level2 event files in pointing mode.
    wtevtfile = glob.glob('./Stage2/*wt*po_cl.evt')
    pcevtfile = glob.glob('./Stage2/*pc*po_cl.evt')

    if (wtevtfile == []) and (pcevtfile == []):
        print '@@@ There is neither WT nor PC mode pointing event file!'
        sys.exit()
    elif (pcevtfile == []) and (wtevtfile != []):
        print '@@@ There is only WT mode pointing event file!'
        print 'WT event file: ', wtevtfile
    elif (pcevtfile != []) and (wtevtfile == []):
        print '@@@ There is only PC mode pointing event file!'
        print 'PC event file: ', pcevtfile
    elif (pcevtfile != []) and (wtevtfile != []):
        print "@@@ There is both WT and PC mode pointing data!"
        print 'WT event file: ', wtevtfile
        print 'PC event file: ', pcevtfile

    if wtevtfile != []:
        wtgtifile = glob.glob('./Stage2/*wt*po_clgti.fits')[0]
        hdufile = fits.open(wtgtifile)
        hdr = hdufile[1].header
        nrowswt = hdr['NAXIS2']
        print 'Number of rows in GTI: ', nrowswt

        for ii in np.arange(1, nrowswt + 1):
            if hdufile[1].data[ii - 1][1] - hdufile[1].data[ii - 1][0] > 400:
                fselect_expr = 'fselect %s+1 %swt_gti%s.fits "#row==%s"' % (wtgtifile, obsid, ii, ii)
                print fselect_expr
                p = subprocess.Popen(fselect_expr, shell=True)
                p.communicate()
                mv_expr = 'mv %swt_gti%s.fits WT/GTI' % (obsid, ii)
                os.system(mv_expr)
            else:
                print "GTI is < 400s"

        nrowsforloop = nrowswt + 1
        for ii in np.arange(1, nrowsforloop):
            if hdufile[1].data[ii - 1][1] - hdufile[1].data[ii - 1][0] > 400:
                try:
                    os.mkdir('WT/Snapshot' + str(ii))
                except OSError:
                    print "Can't create WT/Snapshot%s directory, it already exists!" % (str(ii))
                    pass

    # Find the attitudefile
    attfile = glob.glob(datadir + obsid + '/auxil/*pat.fits*')
    # Find the HD file
    hdfile = glob.glob('./Stage2/*xhdtc.hk')
    print("-----------------         Input Files        ----------------")
    print("Attitude File Name : "), attfile
    print("HD File Name       : "), hdfile

    if wtevtfile != []:
        os.chdir('WT')
        for snapshotfolder in glob.glob('Snapshot*'):
            print snapshotfolder
            sn = snapshotfolder[-1]  # This will be an integer.
            os.chdir(str(snapshotfolder))

            if exprwtgrade == '0-2':
                os.system('cp ' + caldbpath + '/swxwt0to2s6psf1_20131212v001.rmf .')
                os.system('cp ' + caldbpath + '/swxwt0to2s6psf2_20131212v001.rmf .')
                os.system('cp ' + caldbpath + '/swxwt0to2s6psf3_20131212v001.rmf .')
                defwtrmf = 'swxwt0to2s6psf1_20131212v001.rmf'
            elif exprwtgrade == '0':
                os.system('cp ' + caldbpath + '/swxwt0s6psf1_20131212v001.rmf .')
                os.system('cp ' + caldbpath + '/swxwt0s6psf2_20131212v001.rmf .')
                os.system('cp ' + caldbpath + '/swxwt0s6psf3_20131212v001.rmf .')
                defwtrmf = 'swxwt0s6psf1_20131212v001.rmf'

            if wtsrcregionfile == '':
                if pileup == 0:
                    srcregexp = "fk5; circle(%s , %s, %s)" % (tsrcra, tsrcdec, srcouterradius)
                elif pileup == 1:
                    srcregexp = "fk5; annulus(%s , %s, %s, %s)" % (tsrcra, tsrcdec, srcinnerradius, srcouterradius)

                srcregion = open('src.reg', 'w')
                srcregion.write(srcregexp)
                srcregion.close()
            else:
                pass

            bkgregexp = 'fk5; annulus(' + tsrcra + ',' + tsrcdec + ', 165", 200")'

            bkgregion = open('bkg.reg', 'w')
            bkgregion.write(bkgregexp)
            bkgregion.close()

            #usrgtifile = '../GTI/%swt_gti%s.fits' % (obsid, sn)
            os.system('punlearn xrtpipeline')
            xrtpip_str = 'xrtpipeline indir=%s%s outdir=./ steminputs=sw%s srcra=%s ' \
                         'srcdec=%s exprwtgrade=%s createexpomap=yes extractproducts=yes ' \
                         'cleancols=yes chatter=%s exitstage=3 obsmode=%s ' \
                         'wtregionfile=src.reg ' \
                         'wtrmffile=%s psfflag=yes datamode=WT clobber=yes wtbinsize=0.00390625' % (datadir, \
                                                                                                    obsid, obsid,
                                                                                                    str(tsrcra),
                                                                                                    str(tsrcdec),
                                                                                                    exprwtgrade, \
                                                                                                    str(stage3chatter),
                                                                                                    obsmode,
                                                                                                    defwtrmf)
            print os.getcwd()
            print xrtpip_str
            p = subprocess.Popen(xrtpip_str, shell=True)
            p.communicate()

            # Create the background spectrum
            bkgspec = open('bkg_spec_creator.xco', 'w')
            bkgspec.write('test\n')
            bkgspec.write('read event\n')
            bkgspec.write('.\n')
            cleanevent = glob.glob('*po_cl.evt')[0]
            bkgspec.write(cleanevent + '\n')
            bkgspec.write('yes\n')
            bkgspec.write('filter region bkg.reg\n')
            bkgspec.write('extract spectrum\n')
            bkgspec.write('yes\n')
            bkgspec.write('save spectrum bkg.pha\n')
            bkgspec.write('quit\n')
            bkgspec.write('yes\n')
            bkgspec.close()

            os.system('xselect @bkg_spec_creator.xco')

            os.chdir('../')

    if pcevtfile != []:
        pcdir = basedir + '/' + resdir + '/PC'
        os.chdir(pcdir)

        if exprpcgrade == '0-12':
            os.system('cp ' + caldbpath + '/swxpc0to12s6_20130101v014.rmf .')
            defpcrmf = 'swxpc0to12s6_20130101v014.rmf'
        elif exprwtgrade == '0-4':
            os.system('cp ' + caldbpath + '/swxpc0to4s6_20130101v014.rmf .')
            defpcrmf = 'swxpc0to4s6_20130101v014.rmf'
        elif exprwtgrade == '0':
            os.system('cp ' + caldbpath + '/swxpc0s6_20130101v014.rmf .')
            defpcrmf = 'swxpc0s6_20130101v014.rmf'

        if pcsrcregionfile == '':
            if pileup == 0:
                srcregexp = "fk5; circle(%s , %s, %s)" % (tsrcra, tsrcdec, srcouterradius)
            elif pileup == 1:
                srcregexp = "fk5; annulus(%s , %s, %s, %s)" % (tsrcra, tsrcdec, srcinnerradius, srcouterradius)

        bkgregexp = 'fk5; annulus(' + tsrcra + ',' + tsrcdec + ', 220", 250")'

        srcregion = open('src.reg', 'w')
        srcregion.write(srcregexp)
        srcregion.close()

        bkgregion = open('bkg.reg', 'w')
        bkgregion.write(bkgregexp)
        bkgregion.close()

        os.system('punlearn xrtpipeline')
        xrtpip_str = 'xrtpipeline indir=%s%s outdir=./ steminputs=sw%s srcra=%s ' \
                     'srcdec=%s exprpcgrade=%s createexpomap=yes extractproducts=yes ' \
                     'cleancols=yes chatter=%s exitstage=3 obsmode=%s ' \
                     'pcregionfile=src.reg ' \
                     'pcrmffile=%s psfflag=yes datamode=PC clobber=yes' % (datadir, \
                                                                           obsid, obsid, str(tsrcra), str(tsrcdec),
                                                                           exprpcgrade, \
                                                                           str(stage3chatter), obsmode, defpcrmf)
        print os.getcwd()
        print xrtpip_str
        p = subprocess.Popen(xrtpip_str, shell=True)
        p.communicate()

        # Create the background spectrum
        bkgspec = open('bkg_spec_creator.xco', 'w')
        bkgspec.write('test\n')
        bkgspec.write('read event\n')
        bkgspec.write('.\n')
        cleanevent = glob.glob('*po_cl.evt')[0]
        bkgspec.write(cleanevent + '\n')
        bkgspec.write('yes\n')
        bkgspec.write('filter region bkg.reg\n')
        bkgspec.write('extract spectrum\n')
        bkgspec.write('yes\n')
        bkgspec.write('save spectrum bkg.pha\n')
        bkgspec.write('quit\n')
        bkgspec.write('yes\n')
        bkgspec.close()

        os.system('xselect @bkg_spec_creator.xco')

        os.chdir('../')


def radial_profile(data, center):
    y, x = indices(data.shape)
    r = sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(int)

    tbin = bincount(r.ravel(), data.ravel())
    cc = data.ravel()/data.ravel()
    nr = bincount(r.ravel(), cc)
    radialprofile = tbin / nr
    return radialprofile

def luminosity(flux, distance):#, mass):
    lum = flux * 4 * pi * (distance)**2.
    #eddlum = lum / 1.28e+38 * (mass / units.solMass)
    print 'Luminosity:', lum.to('erg/s')
    #print 'LEdd: ', eddlum.to()
