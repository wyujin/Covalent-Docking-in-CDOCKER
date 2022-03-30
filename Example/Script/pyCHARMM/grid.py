# pycharmm covalent rigid cdocker test case 

## Import module
import numpy as np
import pycharmm
import pycharmm.lib as lib
import pycharmm.read as read
import pycharmm.grid as grid

from pycharmm.cdocker import Rigid_CDOCKER 


################################################################
##
##		Begin of pyCHARMM Covalent CDOCKER
##
################################################################

## Topology and parameter files 
read.rtf('"../Toppar/top_all36_prot.rtf"')
read.prm('"../Toppar/par_all36m_prot.prm"', flex = True)

## Grid box information 
xcen = np.loadtxt('xcen', dtype = float)
ycen = np.loadtxt('ycen', dtype = float)
zcen = np.loadtxt('zcen', dtype = float)
maxlen = np.loadtxt('maxlen', dtype = float)
cov_parameter = np.loadtxt('cov_parameter', dtype = float)
rcta = cov_parameter[0]
rctb = cov_parameter[1]
hmax = cov_parameter[2]

## Build receptor structure
read.psf_card("protein.psf", append = True)
read.pdb("protein.pdb", resid = True)

## Generate soft grid
cdocker = grid.CDOCKER()
settings = {'xCen' : xcen, 'yCen' : ycen, 'zCen' : zcen,
            'rcta' : rcta, 'rctb' : rctb, 'hMax' : hmax,
            'xMax' : maxlen, 'yMax' : maxlen, 'zMax' : maxlen,
            'emax' : 15, 'mine' : -120, 'maxe' : -2, 'flag_grhb' : True,
            'probeFile' : '../Toppar/fftdock_c36prot_cgenff_probes.txt'}
cdocker.setVar(settings)
cdocker.generate() 

## Generate hard grid
settings = {'emax' : 3, 'mine' : -20, 'maxe' : 40}
cdocker.setVar(settings)
cdocker.generate()

## Generate native grid
settings = {'emax' : 100, 'mine' : -100, 'maxe' : 100}
cdocker.setVar(settings)
cdocker.generate()

