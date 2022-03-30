# pycharmm covalent rigid cdocker test case 

## Import module
import numpy as np
import pycharmm
import pycharmm.lib as lib
import pycharmm.read as read
import pycharmm.lingo as lingo
import pycharmm.settings as settings

from pycharmm.cdocker import Rigid_CDOCKER 


################################################################
##
##		Begin of pyCHARMM Covalent CDOCKER
##
################################################################

## Topology and parameter files 
settings.set_bomb_level(-1)
read.rtf('"../Toppar/top_all36_prot.rtf"')
read.rtf('"../Toppar/top_all36_cgenff.rtf"', append = True)
read.prm('"../Toppar/par_all36m_prot.prm"', flex = True)
read.prm('"../Toppar/par_all36_cgenff.prm"', append = True, flex = True)
settings.set_bomb_level(0)
lingo.charmm_script('stream "./ligandrtf"')

## Grid box information 
xcen = np.loadtxt('xcen', dtype = float)
ycen = np.loadtxt('ycen', dtype = float)
zcen = np.loadtxt('zcen', dtype = float)

## Covalent CDOCKER standard docking protocol 
sortedResult, dockResult = Rigid_CDOCKER(xcen = xcen, ycen = ycen, zcen = zcen, 
                           softGridFile = 'grid-emax-15-mine--120-maxe--2.bin',
                           hardGridFile = 'grid-emax-3-mine--20-maxe-40.bin',
                           nativeGridFile = 'grid-emax-100-mine--100-maxe-100.bin',
                           flag_save_top = False, flag_save_all = False, 
                           flag_use_hbond = True, flag_grid = True,
                           flag_fast_grid = True, flag_delete_grid = False, 
                           flag_fast_placement = False, confDir = './optimized/')

