#
# NuGridpy - Tools for accessing and visualising NuGrid data.
#
# Copyright 2007 - 2014 by the NuGrid Team.
# All rights reserved. See LICENSE.
#

"""
MESA output data loading and plotting

v0.2, 15OCT2012: NuGrid collaboration
(Sam Jones, Michael Bennett, Daniel Conti, William Hillary,
Falk Herwig, Christian Ritter)

v0.1, 23JUN2010: Falk Herwig

    mesa.py provides tools to get MESA stellar evolution data output
    into your favourite python session.  In the LOGS directory MESA
    outputs two types of files: history.data or star.log is a time
    evolution output, printing one line per so many cycles (e.g. each
    cycle) of all sorts of things.  profilennn.data or lognnn.data
    files are profile data files.  nnn is the number of profile.data
    or log.data files that is translated into model cycles in the
    profiles.index file.

    MESA allows users to freely define what should go into these two
    types of outputs, which means that column numbers can and do
    change.  mesa.py reads in both types of files and present them (as
    well as any header attributes) as arrays that can be referenced by
    the actual column name as defined in the header section of the
    files.  mesa.py then defines a (hopefully growing) set of standard
    plots that make use of the data just obtained.

    mesa.py is organised as a module that can be imported into any
    python or ipython session. It is related to nugridse.py which is a
    similar module to deal with 'se' output, used by the NuGrid
    collaboration.  mesa.py does not need se libraries.  The 'se' output
    files that can be written with MESA can be read and processed with
    the nugridse.py tool.

    mesa.py is providing two class objects, mesa_profile and
    history_data.  The first makes profile data available, the second
    reads and plots the history.data or star.log file.  Note that
    several instances of these can be initiated within one session and
    data from different instances (i.e. models, tracks etc) can be
    overplotted.

    Here is how a simple session could look like that is plotting an
    HRD (We prefer to load ipython with matplotlib and numpy support
    via the alias:

    alias mpython='ipython --pylab')

    >>> import mesa as ms
    >>> help ms
    ------> help(ms)

    >>> s=ms.history_data('.')
    >>> s.hrd()

    In order to find out what header attributes and columns are
    available in history.data or star.log use:

    >>> s.header_attr
    {'burn_min1': 50.0,
     'burn_min2': 1000.0,
     'c12_boundary_limit': 0.0001,
     'h1_boundary_limit': 0.0001,
     'he4_boundary_limit': 0.0001,
     'initial_mass': 2.0,
     'initial_z': 0.01}

    >>> s.cols
    {'center_c12': 38,
     'center_h1': 36,
     'center_he4': 37,
     ...

    In order to read the profile data from the first profile.data file
    in profiles.index, and then get the mass and temperature out and
    finally plot them try.  Typically you will have already a
    Kippenhahn diagram as a function of model number in front of you,
    and you want to access profile information for a given cycle
    number.  Typically you do not have profiles for all cycle
    numbers.  The best way to start a profile instance is with
    num_type='nearest_model' (check the docstring for other ways to
    select profiles for a profile instance):

    >>> a1=ms.mesa_profile('LOGS',59070)
    2001 in profiles.index file ...
    Found and load nearest profile for cycle 59000
    reading LOGS/profile1801.data ...
    Closing profile tool ...

    >>> T=a1.get('temperature')
    >>> mass=a1.get('mmid')
    >>> plot(mass,T)
    [<matplotlib.lines.Line2D object at 0x8456ed0>]

    Or, you could have had it easier in the following way:

    >>> a1.plot('mass','c12',logy=True,shape='-',legend='$^{12}\mathrm{C}$')

    where the superclass plot method interprets data column headers
    correctly and does all the work for you.

    Of course, a1.cols etc are available here as well and many other
    things. E.g. a.model contains an array with all the models for
    which profile.data or log.data  are available. You may initiate a profile object
    with a model number:

    >>> a2=ms.mesa_profile('.',55000,num_type='model')
    100 in profiles.index file ...
    reading ./profile87.data ...

    a1.log_ind (for any profile instance) provides a map of model
    number to profile file number.
    a1.cols and a1.header_attr gives the column names and header attributes.

"""
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from builtins import zip
from builtins import str
from builtins import range
from past.utils import old_div
import numpy as np
import numpy as np
import matplotlib
import matplotlib.pylab as pyl
import matplotlib.pyplot as pl
from matplotlib.collections import Collection
from matplotlib.artist import allow_rasterization
from matplotlib.patches import PathPatch
import os
import sys

from . import ascii_table
from . import utils as u
from .data_plot import *


def set_nugrid_path(path):
    """
        This function sets the path to the NuGrid VOSpace directory as a
        global variable, so that it need only be set once during an inter-
        active session.
        """
    global nugrid_path
    nugrid_path=path

def set_nice_params():
    fsize=18

    params = {'axes.labelsize':  fsize,
    #    'font.family':       'serif',
    'font.family':        'Times New Roman',
    'figure.facecolor':  'white',
    'text.fontsize':     fsize,
    'legend.fontsize':   fsize,
    'xtick.labelsize':   fsize*0.8,
    'ytick.labelsize':   fsize*0.8,
    'ytick.minor.pad': 8,
    'ytick.major.pad': 8,
    'xtick.minor.pad': 8,
    'xtick.major.pad': 8,
    'text.usetex':       False,
    'lines.markeredgewidth': 0}
    pl.rcParams.update(params)

class mesa_profile(DataPlot):
    """
    read profiles.index and prepare reading MESA profile files

    starts with reading profiles.index and creates hash array
    profile.data can then be accessed via prof_plot

    Parameters
    ----------
    sldir : string
        Directory path of LOGS.
    num : integer
        by default this is the i.  profile file (profile.data or
        log.data) available (e.g. num=1 is the 1. available profile
        file), however if you give
    num_type : string, optional
        If 'model' (exact) or 'nearest_model': get the profile
        profile.data file for model (or cycle number) used by the
        stellar evolution code

        If 'profile_num': num will be interpreted as the
        profile.data or log.data number profile_num (profile_num is
        the number that appears in the file names of type
        profile23.data or log23.data)

        If 'profiles_i': the ith file in profiles.index file

        If 'explicit': the exact file path of the profile file
        should be given in the variable give_filename.

        The default is "nearest_model".
    prof_ind_name : string, optional
        Use this optional argument if the profiles.index file has
        an alternative name, for example, do
        superpro=ms.profile('LOGS',1,prof_ind_name='super.prof')
    profile_prefix : string, optional
        Prifix in the profile name.  The default is "profile".
    data_suffix : string, optional
        Optional arguments that allow you to change the defaults
        for the profile.data profile files.
    mass : integer or float, optional
        The user may select a mass and metallicity instead of providing
        the sldir explicitly, if they are using the VOSpace data. If mass
        is provided then Z should also be provided.
        The default is None (i.e. user gives sldir explicitly)
    Z : float, optional
        See 'mass' above.
        The default is None (i.e. user gives sldir explicitly)
    data_set : string, optional
        Coose your data  of 'set1' or 'set1ext'.  The default is 'set1ext'.


    Examples
    --------
    initialise a mesa_profile instance for cycle 2000 like this:

    >>>my_profile1=ms.mesa_profile('LOGS',2000)

    or like this:

    >>>my_profile2=ms.mesa_profile(mass=2,Z=0.01,num=2000)
    """

    sldir = ''

    def __init__(self, sldir='./LOGS', num=1, num_type='nearest_model',
                 prof_ind_name='profiles.index',
                 profile_prefix='profile', data_suffix='.data', mass=None,
                 Z=None,give_filename=None,data_set='set1ext'):
        """
        read a profile.data profile file

        Parameters
        ----------
        sldir : string
            Directory path of LOGS.
        num : integer
            by default this is the i.  profile file (profile.data or
            log.data) available (e.g. num=1 is the 1. available profile
            file), however if you give
        num_type : string, optional
            If 'model' (exact) or 'nearest_model': get the profile
            profile.data file for model (or cycle number) used by the
            stellar evolution code

            If 'profile_num': num will be interpreted as the
            profile.data or log.data number profile_num (profile_num is
            the number that appears in the file names of type
            profile23.data or log23.data)

            If 'profiles_i': the ith file in profiles.index file

            If 'explicit': the exact file path of the profile file
            should be given in the variable give_filename.

            The default is "nearest_model".
        prof_ind_name : string, optional
            Use this optional argument if the profiles.index file has
            an alternative name, for example, do
            superpro=ms.profile('LOGS',1,prof_ind_name='super.prof')
        profile_prefix : string, optional
            Prifix in the profile name.  The default is "profile".
        data_suffix : string, optional
            Optional arguments that allow you to change the defaults
            for the profile.data profile files.
        mass : integer or float, optional
            The user may select a mass and metallicity instead of providing
            the sldir explicitly, if they are using the VOSpace data. If mass
            is provided then Z should also be provided.
            The default is None (i.e. user gives sldir explicitly)
        Z : float, optional
            See 'mass' above.
            The default is None (i.e. user gives sldir explicitly)

        Examples
        --------
        initialise a mesa_profile instance for cycle 2000 like this:

        >>>my_profile1=ms.mesa_profile('LOGS',2000)

        or like this:

        >>>my_profile2=ms.mesa_profile(mass=2,Z=0.01,num=2000)

        """

        self.prof_ind_name = prof_ind_name
        self.sldir         = sldir

        # seeker to find the data requested on VOspace:
        if mass is not None and Z is not None:
            try:
                print('nugrid_path = '+nugrid_path)
            except:
                raise IOError("nugrid_path has not been set. This is the path to the NuGrid VOSpace, e.g. /tmp/NuGrid. Set this using mesa.set_nugrid_path('path')")

            # which set? [find nearest]

            if (data_set=='set1ext'):
               setsZs=[0.02,0.01,6.e-3,1.e-3,1.e-4]
               setsnames=['set1.2','set1.1','set1.3a','set1.4a','set1.5a']
            elif (data_set=='set1'):
               setsZs=[0.02,0.01]
               setsnames=['set1.2','set1.1']
            else:
               raise IOError("Sorry. Requested data_set not available. Choose between set1ext and set1.")                  
       
            idx=np.abs(np.array(setsZs)-Z).argmin()
            setname=setsnames[idx]
            realZ=setsZs[idx]

            print('closest set is '+setname+' (Z = '+str(realZ)+')')

            mod_dir = nugrid_path+'/data/'+data_set+'/'+setname+'/see_wind/'
            if not os.path.exists(mod_dir):
                print('mod_dir = ', mod_dir)
                raise IOError("The data does not seem to be here. Please check that the NuGrid VOSpace is mounted and nugrid_path has been set correctly using mesa.set_nugrid_path('path')'.")

            # which mass? [find nearest]
            list=[el for el in os.listdir(mod_dir) if el[0]=='M']
            if len(list) == 0:
                raise IOError("Sorry. There is no data available for this set at present: "+mod_dir)
            setmasses=[el[1:el.index('Z')] for el in list]
            for i in range(len(setmasses)):
                if setmasses[i][-1]=='.': setmasses[i]=setmasses[i][:-1]
                setmasses[i] = float(setmasses[i])
            idx2=np.abs(np.array(setmasses)-mass).argmin()
            modname=list[idx2]
            realmass=setmasses[idx2]

            print('closest mass is '+str(realmass))

            mod_dir+=modname
            if 'LOGS' not in os.listdir(mod_dir):
                raise IOError("No 'LOGS' directory for this model. It may have been computed with the Geneva code. Try nugridse.py to explore the see_wind data for this model.")
            else:
                self.sldir=mod_dir+'/LOGS'
                sldir = mod_dir+'/LOGS'

        if give_filename is not None and num_type is not 'explicit':
            raise KeyError("Exact filename given but num_type is not explicit.")

        if num_type is 'nearest_model' or num_type is 'model':
            self._profiles_index()
        if num_type is 'nearest_model':
            amods=array(self.model)
            if amods[0]>num:
                 num = amods[0]  
            elif amods[-1]<num:
                 num = amods[-1]
            else:
                 nearmods=[where(amods<=num)[0][-1],where(amods>=num)[0][0]]
                 sometable={}
                 for thing in nearmods:
                     sometable[abs(self.model[thing]-num)]=thing
                 nearest = min(abs(self.model[nearmods[0]]-num),\
                   abs(self.model[nearmods[1]]-num))
                 num = self.model[sometable[nearest]]
            print('Found and load nearest profile for cycle '+str(num))
            num_type = 'model'
        if num_type is 'model':
            try:
                log_num=self.log_ind[num]
            except KeyError:
                print('There is no profile file for this model')
                print("You may retry with num_type='nearest_model'")
                return
        elif num_type is 'profiles_i':
            log_num=self._log_file_ind(num)
            if log_num == -1:
                print("Could not find a profile file with that number")
                return
        elif num_type is 'profile_num':
            log_num = num
        elif num_type is 'explicit':
            pass
        else:
            print('unknown num_type')
            return

        if num_type is 'explicit':
            filename = sldir+'/'+give_filename
            if not os.path.exists(filename):
                print('error: file '+give_filename+' not found in '+sldir)
        else:
            filename=self.sldir+'/'+profile_prefix+str(log_num)+data_suffix
            if not os.path.exists(filename):
                profile_prefix='log'
                filename=self.sldir+'/'+profile_prefix+str(log_num)+data_suffix
                if not os.path.exists(filename):
                    print('error: no profile.data file found in '+sldir)
                    print('error: no log.data file found in '+sldir)


        print('reading '+filename+' ...')
        header_attr = _read_mesafile(filename,only='header_attr')
        num_zones=int(header_attr['num_zones'])
        header_attr,cols,data = _read_mesafile(filename,data_rows=num_zones,only='all')

        self.cols        = cols
        self.header_attr = header_attr
        self.data        = data


    def __del__(self):
        print('Closing profile tool ...')

    def _profiles_index(self):
        """
        read profiles.index and make hash array

        Notes
        -----
        sets the attributes.

        log_ind : hash array that returns profile.data or log.data
        file number from model number.

        model : the models for which profile.data or log.data is
        available

        """

        prof_ind_name = self.prof_ind_name

        f = open(self.sldir+'/'+prof_ind_name,'r')
        line = f.readline()
        numlines=int(line.split()[0])
        print(str(numlines)+' in profiles.index file ...')

        model=[]
        log_file_num=[]
        for line in f:
            model.append(int(line.split()[0]))
            log_file_num.append(int(line.split()[2]))

        log_ind={}    # profile.data number from model
        for a,b in zip(model,log_file_num):
            log_ind[a] = b

        self.log_ind=log_ind
        self.model=model

# let's start with functions that aquire data

    def _log_file_ind(self,inum):
        """
        Information about available profile.data or log.data files.

        Parameters
        ----------
        inum : integer
            Attempt to get number of inum's profile.data file.
            inum_max: max number of profile.data or log.data files
            available

        """

        self._profiles_index()
        if inum <= 0:
            print("Smallest argument is 1")
            return

        inum_max = len(self.log_ind)
        inum -= 1

        if inum > inum_max:
            print('There are only '+str(inum_max)+' profile file available.')
            log_data_number = -1
            return log_data_number
        else:
            log_data_number=self.log_ind[self.model[inum]]
            print('The '+str(inum+1)+'. profile.data file is '+ \
                  str(log_data_number))
            return log_data_number

    def get(self,str_name):
        """
        return a column of data with the name str_name.

        Parameters
        ----------
        str_name : string
            Is the name of the column as printed in the
            profilennn.data or lognnn.data file; get the available
            columns from self.cols (where you replace self with the
            name of your instance)

        """

        column_array = self.data[:,self.cols[str_name]-1].astype('float')
        return column_array

    def write_PROM_HOTB_progenitor(self,name,description):
        """
        Write a progenitor file for the PROMETHEUS/HBOT supernova code.

        Parameters
        ----------
        name : string
            File name for the progenitor file
        description : string
            Information to be written into the file header.
        """
        try:
            from ProgenitorHotb_new import ProgenitorHotb_new
        except ImportError:
            print('Module ProgenitorHotb_new not found.')
            return
        import astronomy as ast
        nz=len(self.get('mass'))
        prog=ProgenitorHotb_new(nz)

        prog.header = '#'+description+'\n'

        prog.xzn    = self.get('rmid')[::-1]*ast.rsun_cm
        prog.massb  = self.get('mass')[::-1]
        prog.r_ob   = max(self.get('radius'))*ast.rsun_cm
        prog.temp   = 10.**self.get('logT')[::-1]*8.620689655172413e-11 # in MeV
        prog.stot   = self.get('entropy')[::-1]
        prog.ye     = self.get('ye')[::-1]
        prog.densty = 10.**self.get('logRho')[::-1]
        prog.press  = 10.**self.get('logP')[::-1]
        prog.eint   = self.get('energy')[::-1]
        prog.velx   = self.get('velocity')[::-1]

        nuclei=['neut','h1','he4','c12','o16','ne20','mg24','si28','s32',
                'ar36','ca40','ti44','cr48','fe52','fake']

        for i in range(len(nuclei)):
            if nuclei[i] == 'fake':
                ni56 = self.get('fe56')+self.get('cr56')
                prog.xnuc[:,i] = ni56[::-1]
            else:
                prog.xnuc[:,i] = self.get(nuclei[i])[::-1]

        prog.write(name)

    def write_STELLA_model(self,name):
        """
        Write an initial model in a format that may easily be read by the
        radiation hydrodynamics code STELLA.

        Parameters
        ----------
        name : string
            an identifier for the model. There are two output files from
            this method, which will be <name>.hyd and <name>.abn, which
            contain the profiles for the hydro and abundance variables,
            respectively.

        """

        # Hydro variables:
        zn = np.array(self.get('zone'),np.int64)
        Mr = self.get('mass')[::-1]
        dM = 10. ** self.get('logdq')[::-1] * self.header_attr['star_mass']
        R = self.get('radius')[::-1] * ast.rsun_cm
        dR = np.insert( np.diff(R), 0, R[0] )
        Rho = 10. ** self.get('logRho')[::-1]
        PRE = 10. ** self.get('logP')[::-1]
        T = 10. ** self.get('logT')[::-1]
        V = self.get('velocity')[::-1]

        # Abundances:
        def make_list(element,lowA,highA):
            l = []
            for i in range(lowA,highA+1):
                l.append(element+str(i))
            return l

        abun_avail = list(self.cols.keys())

        def elemental_abund(ilist,abun_avail):
            X = np.zeros(len(self.get('mass')))
            for a in ilist:
                if a in abun_avail:
                    X += self.get(a)[::-1]

            return X

        iH = ['h1','h2','prot']
        XH  = elemental_abund(iH, abun_avail)
        XHe = elemental_abund(make_list('he',1,5), abun_avail)
        XC  = elemental_abund(make_list('c',11,15), abun_avail)
        XN  = elemental_abund(make_list('n',12,16), abun_avail)
        XO  = elemental_abund(make_list('o',13,20), abun_avail)
        XNe = elemental_abund(make_list('ne',17,25), abun_avail)
        XNa = elemental_abund(make_list('na',20,25), abun_avail)
        XMg = elemental_abund(make_list('mg',21,28), abun_avail)
        XAl = elemental_abund(make_list('al',21,30), abun_avail)
        XSi = elemental_abund(make_list('si',25,34), abun_avail)
        XS  = elemental_abund(make_list('s',28,38), abun_avail)
        XAr = elemental_abund(make_list('ar',32,46), abun_avail)
        XCa = elemental_abund(make_list('ca',36,53), abun_avail)
        XFe = elemental_abund(make_list('fe',50,65), abun_avail)
        XCo = elemental_abund(make_list('co',52,66), abun_avail)
        XNi = elemental_abund(make_list('ni',54,71), abun_avail)

        XNi56 = self.get('ni56')

        # Write the output files:
        file_hyd = name+'.hyd'
        file_abn = name+'.abn'

        f = open(file_hyd,'w')
        # write header:
        f.write('  0.000E+00\n')
        f.write('# No.')
        f.write('Mr'.rjust(28)+
                'dM'.rjust(28)+
                'R'.rjust(28)+
                'dR'.rjust(28)+
                'Rho'.rjust(28)+
                'PRE'.rjust(28)+
                'T'.rjust(28)+
                'V'.rjust(28)+
                '\n')
        # write data:
        for i in range(len(zn)):
            f.write( str(zn[i]).rjust(5) +
                    '%.16E'.rjust(11) %Mr[i] +
                    '%.16E'.rjust(11) %dM[i] +
                    '%.16E'.rjust(11) %R[i] +
                    '%.16E'.rjust(11) %dR[i] +
                    '%.16E'.rjust(11) %Rho[i] +
                    '%.16E'.rjust(11) %PRE[i] +
                    '%.16E'.rjust(11) %T[i] +
                    '%.16E'.rjust(11) %V[i] +
                    '\n')

        f.close()

        f = open(file_abn,'w')
        # write header:
        f.write('# No.')
        f.write('Mr'.rjust(28)+
                'H'.rjust(28)+
                'He'.rjust(28)+
                'C'.rjust(28)+
                'N'.rjust(28)+
                'O'.rjust(28)+
                'Ne'.rjust(28)+
                'Na'.rjust(28)+
                'Mg'.rjust(28)+
                'Al'.rjust(28)+
                'Si'.rjust(28)+
                'S'.rjust(28)+
                'Ar'.rjust(28)+
                'Ca'.rjust(28)+
                'Fe'.rjust(28)+
                'Co'.rjust(28)+
                'Ni'.rjust(28)+
                'X(56Ni)'.rjust(28)+
                '\n')
        # write data:
        for i in range(len(zn)):
            f.write( str(zn[i]).rjust(5) +
                    '%.16E'.rjust(11) %Mr[i] +
                    '%.16E'.rjust(11) %XH[i] +
                    '%.16E'.rjust(11) %XHe[i] +
                    '%.16E'.rjust(11) %XC[i] +
                    '%.16E'.rjust(11) %XN[i] +
                    '%.16E'.rjust(11) %XO[i] +
                    '%.16E'.rjust(11) %XNe[i] +
                    '%.16E'.rjust(11) %XNa[i] +
                    '%.16E'.rjust(11) %XMg[i] +
                    '%.16E'.rjust(11) %XAl[i] +
                    '%.16E'.rjust(11) %XSi[i] +
                    '%.16E'.rjust(11) %XS[i] +
                    '%.16E'.rjust(11) %XAr[i] +
                    '%.16E'.rjust(11) %XCa[i] +
                    '%.16E'.rjust(11) %XFe[i] +
                    '%.16E'.rjust(11) %XCo[i] +
                    '%.16E'.rjust(11) %XNi[i] +
                    '%.16E'.rjust(11) %XNi56[i] +
                    '\n')

    def write_LEAFS_model(self,nzn=30000000,dr=5.e4,
                          rhostrip=5.e-4):
        """
        write an ascii file that will be read by Sam's version of
        inimod.F90 in order to make an initial model for LEAFS
        """

        from scipy import interpolate

        ye = self.get('ye')
        newye=[]
        rho = 10.**self.get('logRho')[::-1] # centre to surface
        # get index to strip all but the core:
        idx = np.abs(rho - rhostrip).argmin() + 1
        rho = rho[:idx]
        rhoc = rho[0]
        rad  = 10.**self.get('logR') * ast.rsun_cm
        rad = rad[::-1][:idx]
        ye = ye[::-1][:idx]

        print('there will be about ',old_div(rad[-1], dr), 'mass cells...')

        # add r = 0 point to all arrays
        rad = np.insert(rad,0,0)
        ye  = np.insert(ye,0,ye[0])
        rho = np.insert(rho,0,rho[0])

        print(rad)

        # interpolate
        fye  = interpolate.interp1d(rad,ye)
        frho = interpolate.interp1d(rad,rho)
        newye = []
        newrho = []
        newrad = []

        Tc   = 10.**self.get('logT')[-1]

        for i in range(nzn):
            if i * dr > rad[-1]: break
            newye.append(fye( i * dr ))
            newrho.append(frho( i * dr ))
            newrad.append( i * dr )

        f = open('M875.inimod','w')

        f.write(str(Tc)+' \n')
        f.write(str(rhoc)+' \n')
        for i in range(len(newye)):
            f.write(str(i+1)+'  '+str(newrad[i])+'  '+\
                    str(newrho[i])+'  '+str(newye[i])+' \n')

        f.close()

    def energy_profile(self,ixaxis):
        """
            Plot radial profile of key energy generations eps_nuc,
            eps_neu etc.

            Parameters
            ----------
            ixaxis : 'mass' or 'radius'
        """

        mass = self.get('mass')
        radius = self.get('radius') * ast.rsun_cm
        eps_nuc = self.get('eps_nuc')
        eps_neu = self.get('non_nuc_neu')

        if ixaxis == 'mass':
            xaxis = mass
            xlab = 'Mass / M$_\odot$'
        else:
            xaxis = old_div(radius, 1.e8) # Mm
            xlab = 'radius / Mm'

        pl.plot(xaxis, np.log10(eps_nuc),
                'k-',
                label='$\epsilon_\mathrm{nuc}>0$')
        pl.plot(xaxis, np.log10(-eps_nuc),
                'k--',
                label='$\epsilon_\mathrm{nuc}<0$')
        pl.plot(xaxis, np.log10(eps_neu),
                'r-',
                label='$\epsilon_\\nu$')

        pl.xlabel(xlab)
        pl.ylabel('$\log(\epsilon_\mathrm{nuc},\epsilon_\\nu)$')
        pl.legend(loc='best').draw_frame(False)

class history_data(DataPlot):
    """
    read history.data or star.log MESA output and plot various things,
    including HRD, Kippenhahn etc

    Parameters
    ----------
    sldir : string
        which LOGS directory.
    slname : string, optional
        If star.log is available instead, star.log file is read, this
        is an optional argument if history.data or star.log file has
        an alternative name.  The default is "history.data".
    clean_starlog : boolean, optional
        Request new cleaning of history.data or star.log, makes
        history.datasa or star.logsa which is the file that is actually
        read and plotted.  The default is False.
    mass : integer or float, optional
        The user may select a mass and metallicity instead of providing
        the sldir explicitly, if they are using the VOSpace data. If mass
        is provided then Z should also be provided.
        The default is None (i.e. user gives sldir explicitly)
    Z : float, optional
        See 'mass' above.
        The default is None (i.e. user gives sldir explicitly)
    data_set : string, optional
        Coose your data  of 'set1' or 'set1ext'.  The default is 'set1ext'.

    Examples
    --------
    use like this:

    >>> another=ms.history_data('LOGS',slname='anothername')

    or this:

    >>> ms.set_nugrid_path('/tmp/NuGrid')
    >>> anotherone=ms.history_data(mass=2,Z=0.01)
    """

    sldir  = ''
    slname = ''
    header_attr = []
    cols = []

    def __init__(self, sldir='./LOGS', slname='history.data',
                 clean_starlog=False, mass=None, Z=None,data_set='set1ext'):
        self.sldir  = sldir
        self.slname = slname
        self.clean_starlog  = clean_starlog

        # seeker to find the data requested on VOspace:
        if mass is not None and Z is not None:
            try:
                print('nugrid_path = '+nugrid_path)
            except:
                raise IOError("nugrid_path has not been set. This is the path to the NuGrid VOSpace, e.g. /tmp/NuGrid. Set this using mesa.set_nugrid_path('path')")

            # which set? [find nearest]

            if (data_set=='set1ext'):
               setsZs=[0.02,0.01,6.e-3,1.e-3,1.e-4]
               setsnames=['set1.2','set1.1','set1.3a','set1.4a','set1.5a']
            elif (data_set=='set1'):
               setsZs=[0.02,0.01]
               setsnames=['set1.2','set1.1']
            else:
               raise IOError("Sorry. Requested data_set not available. Choose between set1ext and set1.")                  

            idx=np.abs(np.array(setsZs)-Z).argmin()
            setname=setsnames[idx]
            realZ=setsZs[idx]

            print('closest set is '+setname+' (Z = '+str(realZ)+')')

            mod_dir = nugrid_path+'/data/'+data_set+'/'+setname+'/see_wind/'
            if not os.path.exists(mod_dir):
                print('mod_dir = ', mod_dir)
                raise IOError("The data does not seem to be here. Please check that the NuGrid VOSpace is mounted and nugrid_path has been set correctly using mesa.set_nugrid_path('path')'.")

            # which mass? [find nearest]
            list=[el for el in os.listdir(mod_dir) if el[0]=='M']
            if len(list) == 0:
                raise IOError("Sorry. There is no data available for this set at present: "+mod_dir)

            setmasses=[el[1:el.index('Z')] for el in list]
            for i in range(len(setmasses)):
                if setmasses[i][-1]=='.': setmasses[i]=setmasses[i][:-1]
                setmasses[i] = float(setmasses[i])
            idx2=np.abs(np.array(setmasses)-mass).argmin()
            modname=list[idx2]
            realmass=setmasses[idx2]

            print('closest mass is '+str(realmass))

            mod_dir+=modname
            if 'LOGS' not in os.listdir(mod_dir):
                raise IOError("No 'LOGS' directory for this model. It may have been computed with the Geneva code. Try nugridse.py to explore the see_wind data for this model.")
            else:
                self.sldir=mod_dir+'/LOGS'
                sldir=mod_dir+'/LOGS'

        if not os.path.exists(self.sldir+'/'+self.slname):
            if not os.path.exists(self.sldir+'/'+'star.log'):
                print('error: no history.data file found in '+sldir)
                print('error: no star.log file found in '+sldir)
            else:
                self.slname='star.log'
                self._read_starlog()
        else:
            self._read_starlog()

    def __del__(self):
        print('Closing', self.slname,' tool ...')

# let's start with functions that aquire data
    def _read_starlog(self):
        """ read history.data or star.log file again"""

        sldir   = self.sldir
        slname  = self.slname
        slaname = slname+'sa'

        if not os.path.exists(sldir+'/'+slaname):
            print('No '+self.slname+'sa file found, create new one from '+self.slname)
            _cleanstarlog(sldir+'/'+slname)
        else:
            if self.clean_starlog:
                print('Requested new '+self.slname+'sa; create new from '+self.slname)
                _cleanstarlog(sldir+'/'+slname)
            else:
                print('Using old '+self.slname+'sa file ...')

        cmd=os.popen('wc '+sldir+'/'+slaname)
        cmd_out=cmd.readline()
        cnum_cycles=cmd_out.split()[0]
        num_cycles=int(cnum_cycles) - 6

        filename=sldir+'/'+slaname

        header_attr,cols,data = _read_mesafile(filename,data_rows=num_cycles)

        self.cols        = cols
        self.header_attr = header_attr
        self.data        = data

    def get(self, str_name):
        """
        return a column of data with the name str_name.

        Parameters
        ----------
        str_name : string
            The name of the column as printed in history.data or
            star.log get the available columns from self.cols (where
            you replace self with the name of your instance

        """

        column_array = self.data[:,self.cols[str_name]-1].astype('float')
        return column_array

    def CO_ratio(self,ifig,ixaxis):
        """
        plot surface C/O ratio in Figure ifig with x-axis quantity ixaxis

        Parameters
        ----------
        ifig : integer
            Figure number in which to plot
        ixaxis : string
            what quantity is to be on the x-axis, either 'time' or 'model'
            The default is 'model'
        """

        def C_O(model):
            surface_c12=model.get('surface_c12')
            surface_o16=model.get('surface_o16')
            CORatio=old_div((surface_c12*4.),(surface_o16*3.))
            return CORatio

        if ixaxis=='time':
            xax=self.get('star_age')
        elif ixaxis=='model':
            xax=self.get('model_number')
        else:
            raise IOError("ixaxis not recognised")

        pl.figure(ifig)
        pl.plot(xax,C_O(self))

    def hrd(self,ifig=None,label=None,colour=None,s2ms=False,
            dashes=None,**kwargs):
        """
        Plot an HR diagram

        Parameters
        ----------
        ifig : integer or string
            Figure label, if None the current figure is used
            The default value is None.
        lims : list [x_lower, x_upper, y_lower, y_upper]
        label : string
            Label for the model
            The default value is None
        colour : string
            The colour of the line
            The default value is None
        s2ms : boolean, optional
            "Skip to Main Sequence"?
            The default is False.
        dashes : list, optional
            Custom dashing style. If None, ignore.
            The default is None.

        """

#        fsize=18
#
#        params = {'axes.labelsize':  fsize,
#        #    'font.family':       'serif',
#        'font.family':        'Times New Roman',
#        'figure.facecolor':  'white',
#        'text.fontsize':     fsize,
#        'legend.fontsize':   fsize,
#        'xtick.labelsize':   fsize*0.8,
#        'ytick.labelsize':   fsize*0.8,
#        'text.usetex':       False}
#
#        try:
#            pl.rcParams.update(params)
#        except:
#            pass

        if ifig is not None:
            pl.figure(ifig)

        if s2ms:
            h1=self.get('center_h1')
            idx=np.where(h1[0]-h1>=3.e-3)[0][0]
            skip=idx
        else:
            skip=0

        x = self.get('log_Teff')[skip:]
        y = self.get('log_L')[skip:]

        if label is not None:
            if colour is not None:
                line,=pl.plot(x,y,label=label,color=colour,**kwargs)
            else:
                line,=pl.plot(x,y,label=label,**kwargs)
        else:
            if colour is not None:
                line,=pl.plot(x,y,color=colour,**kwargs)
            else:
                line,=pl.plot(x,y,**kwargs)

        if dashes is not None:
            line.set_dashes(dashes)
        if label is not None:
            pl.legend(loc='best').draw_frame(False)


#        pyl.plot(self.data[:,self.cols['log_Teff']-1],\
#                 self.data[:,self.cols['log_L']-1],\
#                 label = "M="+str(self.header_attr['initial_mass'])+", Z="\
#                 +str(self.header_attr['initial_z']))

        pyl.xlabel('$\log T_{\\rm eff}$')
        pyl.ylabel('$\log L$')
        x1,x2=pl.xlim()
        if x2 > x1:
            ax=pl.gca()
            ax.invert_xaxis()
#            self._xlimrev()

    def hrd_key(self, key_str):
        """
        plot an HR diagram

        Parameters
        ----------
        key_str : string
            A label string

        """

        pyl.plot(self.data[:,self.cols['log_Teff']-1],\
                 self.data[:,self.cols['log_L']-1],label = key_str)
        pyl.legend()
        pyl.xlabel('log Teff')
        pyl.ylabel('log L')
        x1,x2=pl.xlim()
        if x2 > x1:
            self._xlimrev()

    def hrd_new(self, input_label="", skip=0):
        """
        plot an HR diagram with options to skip the first N lines and
        add a label string

        Parameters
        ----------
        input_label : string, optional
            Diagram label.  The default is "".
        skip : integer, optional
            Skip the first n lines.  The default is 0.

        """
        xl_old=pyl.gca().get_xlim()
        if input_label == "":
            my_label="M="+str(self.header_attr['initial_mass'])+", Z="+str(self.header_attr['initial_z'])
        else:
            my_label="M="+str(self.header_attr['initial_mass'])+", Z="+str(self.header_attr['initial_z'])+"; "+str(input_label)

        pyl.plot(self.data[skip:,self.cols['log_Teff']-1],self.data[skip:,self.cols['log_L']-1],label = my_label)
        pyl.legend(loc=0)
        xl_new=pyl.gca().get_xlim()
        pyl.xlabel('log Teff')
        pyl.ylabel('log L')
        if any(array(xl_old)==0):
            pyl.gca().set_xlim(max(xl_new),min(xl_new))
        elif any(array(xl_new)==0):
            pyl.gca().set_xlim(max(xl_old),min(xl_old))
        else:
            pyl.gca().set_xlim([max(xl_old+xl_new),min(xl_old+xl_new)])

    def xche4_teff(self,ifig=None,lims=[1.,0.,3.4,4.7],label=None,colour=None,
                   s2ms=True,dashes=None):
        """
        Plot effective temperature against central helium abundance.

        Parameters
        ----------
        ifig : integer or string
            Figure label, if None the current figure is used
            The default value is None.
        lims : list [x_lower, x_upper, y_lower, y_upper]
        label : string
        Label for the model
            The default value is None
        colour : string
            The colour of the line
            The default value is None
        s2ms : boolean, optional
            "Skip to Main Sequence"
            The default is True
        dashes : list, optional
            Custom dashing style. If None, ignore.
            The default is None.

        """
        fsize=18
        params = {'axes.labelsize':  fsize,
        #    'font.family':       'serif',
        'font.family':        'Times New Roman',
        'figure.facecolor':  'white',
        'text.fontsize':     fsize,
        'legend.fontsize':   fsize,
        'xtick.labelsize':   fsize*0.8,
        'ytick.labelsize':   fsize*0.8,
        'text.usetex':       False}

        try:
            pl.rcParams.update(params)
        except:
            pass


        if s2ms:
            h1=self.get('center_h1')
            idx=np.where(h1[0]-h1>=1.e-3)[0][0]
            skip=idx
        else:
            skip=0

        x = self.get('center_he4')[skip:]
        y = self.get('log_Teff')[skip:]
        if ifig is not None:
            pl.figure(ifig)
        if label is not None:
            if colour is not None:
                line,=pl.plot(x,y,label=label,color=colour)
            else:
                line,=pl.plot(x,y,label=label)
            pl.legend(loc='best').draw_frame(False)
        else:
            if colour is not None:
                line,=pl.plot(x,y,color=colour)
            else:
                line,=pl.plot(x,y)

        if dashes is not None:
            line.set_dashes(dashes)

        if label is not None:
            pl.legend(loc='best').draw_frame(False)

        pl.xlim(lims[:2])
        pl.ylim(lims[2:])
        pl.xlabel('$X_{\\rm c}(\,^4{\\rm He}\,)$')
        pl.ylabel('$\log\,T_{\\rm eff}$')


    def tcrhoc(self,ifig=None,lims=[3.,10.,8.,10.],label=None,colour=None,
               dashes=None):
        """
        Central temperature again central density plot

        Parameters
        ----------
        ifig : integer or string
            Figure label, if None the current figure is used
            The default value is None.
        lims : list [x_lower, x_upper, y_lower, y_upper]
        label : string
            Label for the model
            The default value is None
        colour : string
            The colour of the line
            The default value is None
        dashes : list, optional
            Custom dashing style. If None, ignore.
            The default is None.

        """

#        fsize=18
#
#        params = {'axes.labelsize':  fsize,
#            #    'font.family':       'serif',
#        'font.family':        'Times New Roman',
#        'figure.facecolor':  'white',
#        'text.fontsize':     fsize,
#        'legend.fontsize':   fsize,
#        'xtick.labelsize':   fsize*0.8,
#        'ytick.labelsize':   fsize*0.8,
#        'text.usetex':       False}
#
#        try:
#            pl.rcParams.update(params)
#        except:
#            pass

        if ifig is not None:
            pl.figure(ifig)

        if label is not None:
            if colour is not None:
                line,=pl.plot(self.get('log_center_Rho'),self.get('log_center_T'),label=label,
                        color=colour)
            else:
                line,=pl.plot(self.get('log_center_Rho'),self.get('log_center_T'),label=label)

        else:
            if colour is not None:
                line,=pl.plot(self.get('log_center_Rho'),self.get('log_center_T'),
                        color=colour)
            else:
                line,=pl.plot(self.get('log_center_Rho'),self.get('log_center_T'))
        if dashes is not None:
            line.set_dashes(dashes)
        if label is not None:
            pl.legend(loc='best').draw_frame(False)

        pl.xlim(lims[:2])
        pl.ylim(lims[2:])
        pl.xlabel('log $\\rho_{\\rm c}$')
        pl.ylabel('log $T_{\\rm c}$')

    def mdot_t(self,ifig=None,lims=[7.4,2.6,-8.5,-4.5],label=None,colour=None,s2ms=False,
               dashes=None):
        """
        Plot mass loss history as a function of log-time-left

        Parameters
        ----------
        ifig : integer or string
            Figure label, if None the current figure is used
            The default value is None.
        lims : list [x_lower, x_upper, y_lower, y_upper]
        label : string
            Label for the model
            The default value is None
        colour : string
            The colour of the line
            The default value is None
        s2ms : boolean, optional
            "skip to main sequence"
        dashes : list, optional
            Custom dashing style. If None, ignore.
            The default is None.

        """

        fsize=18

        params = {'axes.labelsize':  fsize,
        #    'font.family':       'serif',
        'font.family':        'Times New Roman',
        'figure.facecolor':  'white',
        'text.fontsize':     fsize,
        'legend.fontsize':   fsize,
        'xtick.labelsize':   fsize*0.8,
        'ytick.labelsize':   fsize*0.8,
        'text.usetex':       False}

        try:
            pl.rcParams.update(params)
        except:
            pass

        if ifig is not None:
            pl.figure(ifig)

        if s2ms:
            h1=self.get('center_h1')
            idx=np.where(h1[0]-h1>=3.e-3)[0][0]
            skip=idx
        else:
            skip=0

        gage= self.get('star_age')
        lage=np.zeros(len(gage))
        agemin = max(old_div(abs(gage[-1]-gage[-2]),5.),1.e-10)
        for i in np.arange(len(gage)):
            if gage[-1]-gage[i]>agemin:
                lage[i]=np.log10(gage[-1]-gage[i]+agemin)
            else :
                lage[i]=np.log10(agemin)
        x = lage[skip:]
        y = self.get('log_abs_mdot')[skip:]

        if ifig is not None:
            pl.figure(ifig)
        if label is not None:
            if colour is not None:
                line,=pl.plot(x,y,label=label,color=colour)
            else:
                line,=pl.plot(x,y,label=label)
        else:
            if colour is not None:
                line,=pl.plot(x,y,color=colour)
            else:
                line,=pl.plot(x,y)

        if dashes is not None:
            line.set_dashes(dashes)
        if label is not None:
            pl.legend(loc='best').draw_frame(False)

        pl.xlim(lims[:2])
        pl.ylim(lims[2:])
        pl.ylabel('$\mathrm{log}_{10}(\|\dot{M}\|/M_\odot\,\mathrm{yr}^{-1})$')
        pl.xlabel('$\mathrm{log}_{10}(t^*/\mathrm{yr})$')

    def mcc_t(self,ifig=None,lims=[0,15,0,25],label=None,colour=None,
              mask=False,s2ms=False,dashes=None):
        """
        Plot mass of convective core as a function of time.

        Parameters
        ----------
        ifig : integer or string
            Figure label, if None the current figure is used
            The default value is None.
        lims : list [x_lower, x_upper, y_lower, y_upper]
        label : string
            Label for the model
            The default value is None
        colour : string
            The colour of the line
            The default value is None
        mask : boolean, optional
            Do you want to try to hide numerical spikes in the
            plot?
            The default is False
        s2ms : boolean, optional
            skip to main squence?
        dashes : list, optional
            Custom dashing style. If None, ignore.
            The default is None.

        """

        fsize=18

        params = {'axes.labelsize':  fsize,
        #    'font.family':       'serif',
        'font.family':        'Times New Roman',
        'figure.facecolor':  'white',
        'text.fontsize':     fsize,
        'legend.fontsize':   fsize,
        'xtick.labelsize':   fsize*0.8,
        'ytick.labelsize':   fsize*0.8,
        'text.usetex':       False}

        try:
            pl.rcParams.update(params)
        except:
            pass

        if ifig is not None:
            pl.figure(ifig)

        if s2ms:
            h1=self.get('center_h1')
            idx=np.where(h1[0]-h1>=3.e-3)[0][0]
            skip=idx
        else:
            skip=0

        age= self.get('star_age')
        x1 = old_div(age, 1.e6)
        x2 = old_div(age, 1.e6)
        y1 = self.get('mix_qtop_1')*self.get('star_mass')
        y2 = self.get('mix_qtop_2')*self.get('star_mass')
        mt1 = self.get('mix_type_1')
        mt2 = self.get('mix_type_2')

        x1 = x1[skip:]
        x2 = x2[skip:]
        y1 = y1[skip:]
        y2 = y2[skip:]
        mt1 = mt1[skip:]
        mt2 = mt2[skip:]

        # Mask spikes...
        if mask:
            x1 = np.ma.masked_where(mt1 != 1, x1)
            x2 = np.ma.masked_where(mt2 != 1, x2)
            y1 = np.ma.masked_where(mt1 != 1, y1)
            y2 = np.ma.masked_where(mt2 != 1, y2)

        if ifig is not None:
            pl.figure(ifig)
        if label is not None:
            if colour is not None:
                line,=pl.plot(x1,y1,label=label,color=colour)
                line,=pl.plot(x2,y2,color=colour)
            else:
                line,=pl.plot(x1,y1,label=label)
                line,=pl.plot(x2,y2)
        else:
            if colour is not None:
                line,=pl.plot(x1,y1,color=colour)
                line,=pl.plot(x2,y2,color=colour)
            else:
                line,=pl.plot(x1,y1)
                line,=pl.plot(x2,y2)

        if dashes is not None:
            line.set_dashes(dashes)

        if label is not None:
            pl.legend(loc='best').draw_frame(False)

        pl.xlim(lims[:2])
        pl.ylim(lims[2:])
        pl.ylabel('$M/M_\odot}$')
        pl.xlabel('$t/{\\rm Myr}$')


    def kippenhahn_CO(self, num_frame, xax, t0_model=0,
                      title='Kippenhahn diagram', tp_agb=0.,
                      ylim_CO=[0,0]):
        """
        Kippenhahn plot as a function of time or model with CO ratio

        Parameters
        ----------
        num_frame : integer
            Number of frame to plot this plot into.
        xax : string
            Either model or time to indicate what is to be used on the
            x-axis.
        t0_model : integer, optional
            Model for the zero point in time, for AGB plots this would
            be usually the model of the 1st TP, which can be found with
            the Kippenhahn plot.  The default is 0.
        title : string, optional
            Figure title.  The defalut is "Kippenhahn diagram".
        tp_agb : float, optional
            If >= 0 then,
            ylim=[h1_min*1.-tp_agb/100 : h1_max*1.+tp_agb/100] with
            h1_min, h1_max the min and max H-free core mass coordinate.
            The defalut is 0.
        ylim_CO : list
            if ylim_CO is [0,0], then it is automaticly set.  The
            default is [0,0].

        """

        pyl.figure(num_frame)

        if xax == 'time':
            xaxisarray = self.get('star_age')
        elif xax == 'model':
            xaxisarray = self.get('model_number')
        else:
            print('kippenhahn_error: invalid string for x-axis selction.'+\
                  ' needs to be "time" or "model"')

        t0_mod=xaxisarray[t0_model]

        plot_bounds=True
        try:
            h1_boundary_mass  = self.get('h1_boundary_mass')
            he4_boundary_mass = self.get('he4_boundary_mass')
        except:
            try:
                h1_boundary_mass  = self.get('he_core_mass')
                he4_boundary_mass = self.get('c_core_mass')
            except:
                plot_bounds=False

        star_mass         = self.get('star_mass')
        mx1_bot           = self.get('mx1_bot')*star_mass
        mx1_top           = self.get('mx1_top')*star_mass
        mx2_bot           = self.get('mx2_bot')*star_mass
        mx2_top           = self.get('mx2_top')*star_mass
        surface_c12       = self.get('surface_c12')
        surface_o16       = self.get('surface_o16')

        COratio=old_div((surface_c12*4.),(surface_o16*3.))

        pyl.plot(xaxisarray[t0_model:]-t0_mod,COratio[t0_model:],'-k',label='CO ratio')
        pyl.ylabel('C/O ratio')
        pyl.legend(loc=4)
        if ylim_CO[0] is not 0 and  ylim_CO[1] is not 0:
            pyl.ylim(ylim_CO)
        if xax == 'time':
            pyl.xlabel('t / yrs')
        elif xax == 'model':
            pyl.xlabel('model number')

        pyl.twinx()
        if plot_bounds:
            pyl.plot(xaxisarray[t0_model:]-t0_mod,h1_boundary_mass[t0_model:],label='h1_boundary_mass')
            pyl.plot(xaxisarray[t0_model:]-t0_mod,he4_boundary_mass[t0_model:],label='he4_boundary_mass')
        pyl.plot(xaxisarray[t0_model:]-t0_mod,mx1_bot[t0_model:],',r',label='conv bound')
        pyl.plot(xaxisarray[t0_model:]-t0_mod,mx1_top[t0_model:],',r')
        pyl.plot(xaxisarray[t0_model:]-t0_mod,mx2_bot[t0_model:],',r')
        pyl.plot(xaxisarray[t0_model:]-t0_mod,mx2_top[t0_model:],',r')
        pyl.plot(xaxisarray[t0_model:]-t0_mod,star_mass[t0_model:],label='star_mass')
        pyl.ylabel('mass coordinate')
        pyl.legend(loc=2)
        if tp_agb > 0.:
            h1_min = min(h1_boundary_mass[t0_model:])
            h1_max = max(h1_boundary_mass[t0_model:])
            h1_min = h1_min*(1.-old_div(tp_agb,100.))
            h1_max = h1_max*(1.+old_div(tp_agb,100.))
            print('setting ylim to zoom in on H-burning:',h1_min,h1_max)
            pyl.ylim(h1_min,h1_max)

    def kippenhahn(self, num_frame, xax, t0_model=0,
                   title='Kippenhahn diagram', tp_agb=0., t_eps=5.e2,
                   plot_star_mass=True, symbol_size=8, c12_bm=False,
                   print_legend=True):
        """Kippenhahn plot as a function of time or model.

        Parameters
        ----------
        num_frame : integer
            Number of frame to plot this plot into, if <0 open no new
            figure.
        xax : string
            Either 'model', 'time' or 'logtimerev' to indicate what is
            to be used on the x-axis.
        t0_model : integer, optional
            If xax = 'time' then model for the zero point in time, for
            AGB plots this would be usually the model of the 1st TP,
            which can be found with the Kippenhahn plot.  The default
            is 0.
        title : string, optional
            The figure title.  The default is "Kippenhahn diagram".
        tp_agb : float, optional
            If > 0. then,
            ylim=[h1_min*1.-tp_agb/100 : h1_max*1.+tp_agb/100] with
            h1_min, h1_max the min and max H-free core mass coordinate.
            The default is 0. .
        t_eps : float, optional
            Final time for logtimerev.  The default is '5.e2'.
        plot_star_mass : boolean, optional
            If True, then plot the stellar mass as a line as well.  The
            default is True.
        symbol_size : integer, optional
            Size of convection boundary marker.  The default is 8.
        c12_bm : boolean, optional
            If we plot c12_boundary_mass or not.  The default is False.
        print_legend : boolean, optionla
            Show or do not show legend.  The defalut is True.

        """

        if num_frame >= 0:
            pyl.figure(num_frame)

        t0_mod=[]

        if xax == 'time':
            xaxisarray = self.get('star_age')
            if t0_model > 0:
                ind=self.get('model_number')
                t0_model=where(ind>t0_model)[0][0]
                t0_mod=xaxisarray[t0_model]
            else:
                t0_mod = 0.
            print('zero time is '+str(t0_mod))
        elif xax == 'model':
            xaxisarray = self.get('model_number')
            #t0_mod=xaxisarray[t0_model]
            t0_mod = 0.
        elif xax == 'logtimerev':
            xaxi    = self.get('star_age')
            xaxisarray = np.log10(np.max(xaxi)+t_eps-xaxi)
            t0_mod = 0.
        else:
            print('kippenhahn_error: invalid string for x-axis selction.'+\
                  ' needs to be "time" or "model"')


        plot_bounds=True
        try:
            h1_boundary_mass  = self.get('h1_boundary_mass')
            he4_boundary_mass = self.get('he4_boundary_mass')
            if c12_bm:
                c12_boundary_mass = self.get('c12_boundary_mass')
        except:
            try:
                h1_boundary_mass  = self.get('he_core_mass')
                he4_boundary_mass = self.get('c_core_mass')
                if c12_bm:
                    c12_boundary_mass = self.get('o_core_mass')
            except:
                plot_bounds=False

        star_mass         = self.get('star_mass')
        mx1_bot           = self.get('mx1_bot')*star_mass
        mx1_top           = self.get('mx1_top')*star_mass
        mx2_bot           = self.get('mx2_bot')*star_mass
        mx2_top           = self.get('mx2_top')*star_mass


        if xax == 'time':
            if t0_model>0:
                pyl.xlabel('$t - t_0$ $\mathrm{[yr]}$')
            else:
                pyl.xlabel('t / yrs')
        elif xax == 'model':
            pyl.xlabel('model number')
        elif xax == 'logtimerev':
            pyl.xlabel('$\log(t_{final} - t)$  $\mathrm{[yr]}$')

        pyl.plot(xaxisarray[t0_model:]-t0_mod,mx1_bot[t0_model:],linestyle='None',color='blue',alpha=0.3,marker='o',markersize=symbol_size,label='convection zones')
        pyl.plot(xaxisarray[t0_model:]-t0_mod,mx1_top[t0_model:],linestyle='None',color='blue',alpha=0.3,marker='o',markersize=symbol_size)
        pyl.plot(xaxisarray[t0_model:]-t0_mod,mx2_bot[t0_model:],linestyle='None',color='blue',alpha=0.3,marker='o',markersize=symbol_size)
        pyl.plot(xaxisarray[t0_model:]-t0_mod,mx2_top[t0_model:],linestyle='None',color='blue',alpha=0.3,marker='o',markersize=symbol_size)

        if plot_bounds:
            pyl.plot(xaxisarray[t0_model:]-t0_mod,h1_boundary_mass[t0_model:],color='red',linewidth=2,label='H-free core')
            pyl.plot(xaxisarray[t0_model:]-t0_mod,he4_boundary_mass[t0_model:],color='green',linewidth=2,linestyle='dashed',label='He-free core')
            if c12_bm:
                pyl.plot(xaxisarray[t0_model:]-t0_mod,c12_boundary_mass[t0_model:],color='purple',linewidth=2,linestyle='dotted',label='C-free core')
        if plot_star_mass is True:
            pyl.plot(xaxisarray[t0_model:]-t0_mod,star_mass[t0_model:],label='$M_\star$')
        pyl.ylabel('$m_\mathrm{r}/\mathrm{M}_\odot$')
        if print_legend:
            pyl.legend(loc=2)
        if tp_agb > 0.:
            h1_min = min(h1_boundary_mass[t0_model:])
            h1_max = max(h1_boundary_mass[t0_model:])
            h1_min = h1_min*(1.-old_div(tp_agb,100.))
            h1_max = h1_max*(1.+old_div(tp_agb,100.))
            print('setting ylim to zoom in on H-burning:',h1_min,h1_max)
            pyl.ylim(h1_min,h1_max)

    def t_surfabu(self, num_frame, xax, t0_model=0,
                  title='surface abundance', t_eps=1.e-3,
                  plot_CO_ratio=False):
        """
        t_surfabu plots surface abundance evolution as a function of
        time.

        Parameters
        ----------
        num_frame : integer
            Number of frame to plot this plot into, if <0 don't open
            figure.
        xax : string
            Either model, time or logrevtime to indicate what is to be
            used on the x-axis.
        t0_model : integer, optional
            Model for the zero point in time, for AGB plots this would
            be usually the model of the 1st TP, which can be found with
            the Kippenhahn plot.  The default is 0.
        title : string, optional
            Figure title.  The default is "surface abundance".
        t_eps : float, optional
            Time eps at end for logrevtime.  The default is 1.e-3.
        plot_CO_ratio : boolean, optional
            On second axis True/False.  The default is False.

        """
        if num_frame >= 0:
            pyl.figure(num_frame)

        if xax == 'time':
            xaxisarray = self.get('star_age')[t0_model:]
        elif xax == 'model':
            xaxisarray = self.get('model_number')[t0_model:]
        elif xax == 'logrevtime':
            xaxisarray = self.get('star_age')
            xaxisarray=np.log10(max(xaxisarray[t0_model:])+t_eps-xaxisarray[t0_model:])
        else:
            print('t-surfabu error: invalid string for x-axis selction.'+ \
                  ' needs to be "time" or "model"')

        star_mass         = self.get('star_mass')
        surface_c12       = self.get('surface_c12')
        surface_c13       = self.get('surface_c13')
        surface_n14       = self.get('surface_n14')
        surface_o16       = self.get('surface_o16')

        target_n14 = -3.5


        COratio=old_div((surface_c12*4.),(surface_o16*3.))
        t0_mod=xaxisarray[t0_model]
        log10_c12=np.log10(surface_c12[t0_model:])

        symbs=['k:','-','--','-.','b:','-','--','k-.',':','-','--','-.']

        pyl.plot(xaxisarray,log10_c12,\
                     symbs[0],label='$^{12}\mathrm{C}$')
        pyl.plot(xaxisarray,np.log10(surface_c13[t0_model:]),\
                     symbs[1],label='$^{13}\mathrm{C}$')
        pyl.plot(xaxisarray,np.log10(surface_n14[t0_model:]),\
                     symbs[2],label='$^{14}\mathrm{N}$')
        pyl.plot(xaxisarray,np.log10(surface_o16[t0_model:]),\
                     symbs[3],label='$^{16}\mathrm{O}$')
#                pyl.plot([min(xaxisarray[t0_model:]-t0_mod),max(xaxisarray[t0_model:]-t0_mod)],[target_n14,target_n14])

        pyl.ylabel('mass fraction $\log X$')
        pyl.legend(loc=2)

        if xax == 'time':
            pyl.xlabel('t / yrs')
        elif xax == 'model':
            pyl.xlabel('model number')
        elif xax == 'logrevtime':
            pyl.xlabel('$\\log t-tfinal$')
        if plot_CO_ratio:
            pyl.twinx()
            pyl.plot(xaxisarray,COratio[t0_model:],'-k',label='CO ratio')
            pyl.ylabel('C/O ratio')
            pyl.legend(loc=4)
        pyl.title(title)
        if xax == 'logrevtime':
            self._xlimrev()


# ... end t_surfabu

    def t_lumi(self,num_frame,xax):
        """
        Luminosity evolution as a function of time or model.

        Parameters
        ----------
        num_frame : integer
            Number of frame to plot this plot into.
        xax : string
            Either model or time to indicate what is to be used on the
            x-axis

        """

        pyl.figure(num_frame)

        if xax == 'time':
            xaxisarray = self.get('star_age')
        elif xax == 'model':
            xaxisarray = self.get('model_number')
        else:
            print('kippenhahn_error: invalid string for x-axis selction. needs to be "time" or "model"')


        logLH   = self.get('log_LH')
        logLHe  = self.get('log_LHe')

        pyl.plot(xaxisarray,logLH,label='L_(H)')
        pyl.plot(xaxisarray,logLHe,label='L(He)')
        pyl.ylabel('log L')
        pyl.legend(loc=2)


        if xax == 'time':
            pyl.xlabel('t / yrs')
        elif xax == 'model':
            pyl.xlabel('model number')

    def t_surf_parameter(self, num_frame, xax):
        """
        Surface parameter evolution as a function of time or model.

        Parameters
        ----------
        num_frame : integer
            Number of frame to plot this plot into.
        xax : string
            Either model or time to indicate what is to be used on the
            x-axis

        """

        pyl.figure(num_frame)

        if xax == 'time':
            xaxisarray = self.get('star_age')
        elif xax == 'model':
            xaxisarray = self.get('model_number')
        else:
            print('kippenhahn_error: invalid string for x-axis selction. needs to be "time" or "model"')


        logL    = self.get('log_L')
        logTeff    = self.get('log_Teff')

        pyl.plot(xaxisarray,logL,'-k',label='log L')
        pyl.plot(xaxisarray,logTeff,'-k',label='log Teff')
        pyl.ylabel('log L, log Teff')
        pyl.legend(loc=2)


        if xax == 'time':
            pyl.xlabel('t / yrs')
        elif xax == 'model':
            pyl.xlabel('model number')

    def _kip_vline(self, modstart, modstop, sparse, outfile,
                  xlims=[0.,0.], ylims=[0.,0.], ixaxis='log_time_left',
                  mix_zones=5, burn_zones=50):
        """
        *** DEPRECIATED and hence UNSUPPORTED ***
        This function creates a Kippenhahn plot with energy flux using
        vertical lines, better thermal pulse resolution.

        For a more comprehensive plot, your history.data or star.log
        file should contain columns called "mix_type_n", "mix_qtop_n",
        "burn_type_n" and "burn_qtop_n".  The number of columns
        (i.e. the bbiggest value of n) is what goes in the arguments as
        mix_zones and burn_zones.

        DO NOT WORRY! if you do not have these columns, just leave the
        default values alone and the script should recognise that you
        do not have these columns and make the most detailed plot that
        is available to you.

        Parameters
        ----------
        modstart : integer
            Model from which you want to plot (be careful if your
            history.data or star.log output is sparse...).
        modstop : integer
            Model to which you wish to plot.
        sparse : integer
            x-axis sparsity.
        outfile : string
            'filename + extension' where you want to save the figure.
        xlims, ylims : list, optional
            plot limits, however these are somewhat obsolete now that
            we have modstart and modstop.  Leaving them as 0. is
            probably no slower, and you can always zoom in afterwards
            in mpl.  The default is [0., 0.,].
        ixaxis : string, optional
            Either 'log_time_left', 'age', or 'model_number'.  The
            default "log_time_left".
        mix_zones, burn_zones : integer
            As described above, if you have more detailed output about
            your convection and energy generation boundaries in columns
            mix_type_n, mix_qtop_n, burn_type_n and burn_qtop_n, you
            need to specify the total number of columns for mixing zones
            and burning zones that you have.  Can't work this out from
            your history.data or star.log file?  Check the
            history_columns.list that you used, it'll be the number
            after "mixing regions" and "burning regions".  Can't see
            these columns?  leave it and 2 conv zones and 2 burn zones
            will be drawn using other data that you certainly should
            have in your history.data or star.log file.  The default for
            mix_zones is 5, the defalut for burn_zones is 50.

        """


        xxyy=[self.get('star_age')[modstart:modstop],self.get('star_age')[modstart:modstop]]
        mup = max(float(self.get('star_mass')[0])*1.02,1.0)
        nmodels=len(self.get('model_number')[modstart:modstop])

        Msol=1.98892E+33

        engenstyle = 'full'

        dx = sparse
        x = np.arange(0, nmodels, dx)

        btypemax = 20
        btypemin = -20
        btypealpha=0.

        ########################################################################
        #----------------------------------plot--------------------------------#
        fig = pl.figure()
#       fig.set_size_inches(16,9)
        fsize=15
        ax=pl.axes()

        if ixaxis == 'log_time_left':
        # log of time left until core collapse
            gage= self.get('star_age')
            lage=np.zeros(len(gage))
            agemin = max(old_div(abs(gage[-1]-gage[-2]),5.),1.e-10)
            for i in np.arange(len(gage)):
                if gage[-1]-gage[i]>agemin:
                    lage[i]=np.log10(gage[-1]-gage[i]+agemin)
                else :
                    lage[i]=np.log10(agemin)
            xxx = lage[modstart:modstop]
            print('plot versus time left')
            ax.set_xlabel('$\mathrm{log}_{10}(t^*) \, \mathrm{(yr)}$',fontsize=fsize)
        elif ixaxis =='model_number':
            xxx= self.get('model_number')[modstart:modstop]
            print('plot versus model number')
            ax.set_xlabel('Model number',fontsize=fsize)
        elif ixaxis =='age':
            xxx= old_div(self.get('star_age')[modstart:modstop],1.e6)
            print('plot versus age')
            ax.set_xlabel('Age [Myr]',fontsize=fsize)
        else:
            print('ixaxis must be one of: log_time_left, age or model_number')
            sys.exit()

        if xlims == [0.,0.]:
            xlims[0] = xxx[0]
            xlims[1] = xxx[-1]
        if ylims == [0.,0.]:
            ylims[0] = 0.
            ylims[1] = mup


        print('plotting patches')
        ax.plot(xxx[::dx],self.get('star_mass')[modstart:modstop][::dx],'k-')

        print('plotting abund boundaries')
        ax.plot(xxx,self.get('h1_boundary_mass')[modstart:modstop],label='H boundary')
        ax.plot(xxx,self.get('he4_boundary_mass')[modstart:modstop],label='He boundary')
#       ax.plot(xxx,self.get('c12_boundary_mass')[modstart:modstop],label='C boundary')

        ax.axis([xlims[0],xlims[1],ylims[0],ylims[1]])

        ax.set_ylabel('Mass [M$_\odot$]')

        ########################################################################

        try:
            self.get('burn_qtop_1')
        except:
            engenstyle = 'twozone'
        old_percent = 0
        if engenstyle == 'full':
            for i in range(len(x)):
                # writing status
                percent = int(i*100/(len(x) - 1))
                if percent >= old_percent + 5:
                    sys.stdout.flush()
                    sys.stdout.write("\r creating color map1 " + "...%d%%" % percent)
                    old_percent = percent

                for j in range(1,burn_zones+1):
                    ulimit=self.get('burn_qtop_'+str(j))[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                    if j==1:
                        llimit=0.0
                    else:
                        llimit=self.get('burn_qtop_'+str(j-1))[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                    btype=float(self.get('burn_type_'+str(j))[modstart:modstop][i*dx])
                    if llimit!=ulimit:
                        if btype>0.:
                            #btypealpha = btype/btypemax
                            #ax.axvline(xxx[i*dx],ymin=(llimit-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimit-ylims[0])/(ylims[1]-ylims[0]),color='b',alpha=btypealpha)
                            pass
                        if btype<0.:
                            #btypealpha = (btype/btypemin)/5
                            #ax.axvline(xxx[i*dx],ymin=(llimit-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimit-ylims[0])/(ylims[1]-ylims[0]),color='r',alpha=btypealpha)
                            pass

            print(' \n')

        old_percent = 0

        if engenstyle == 'twozone':
            for i in range(len(x)):
                # writing status
                percent = int(i*100/(len(x) - 1))
                if percent >= old_percent + 5:
                    sys.stdout.flush()
                    sys.stdout.write("\r creating color map1 " + "...%d%%" % percent)
                    old_percent = percent

                llimitl1=old_div(self.get('epsnuc_M_1')[modstart:modstop][i*dx],Msol)
                ulimitl1=old_div(self.get('epsnuc_M_4')[modstart:modstop][i*dx],Msol)
                llimitl2=old_div(self.get('epsnuc_M_5')[modstart:modstop][i*dx],Msol)
                ulimitl2=old_div(self.get('epsnuc_M_8')[modstart:modstop][i*dx],Msol)
                llimith1=old_div(self.get('epsnuc_M_2')[modstart:modstop][i*dx],Msol)
                ulimith1=old_div(self.get('epsnuc_M_3')[modstart:modstop][i*dx],Msol)
                llimith2=old_div(self.get('epsnuc_M_6')[modstart:modstop][i*dx],Msol)
                ulimith2=old_div(self.get('epsnuc_M_7')[modstart:modstop][i*dx],Msol)
                # lower thresh first, then upper thresh:
                #if llimitl1!=ulimitl1:
                    #ax.axvline(xxx[i*dx],ymin=(llimitl1-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimitl1-ylims[0])/(ylims[1]-ylims[0]),color='b',alpha=1.)
                #if llimitl2!=ulimitl2:
                    #ax.axvline(xxx[i*dx],ymin=(llimitl2-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimitl2-ylims[0])/(ylims[1]-ylims[0]),color='b',alpha=1.)
                #if llimith1!=ulimith1:
                    #ax.axvline(xxx[i*dx],ymin=(llimith1-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimith1-ylims[0])/(ylims[1]-ylims[0]),color='b',alpha=4.)
                #if llimith2!=ulimith2:
                    #ax.axvline(xxx[i*dx],ymin=(llimith2-ylims[0])/(ylims[1]-ylims[0]),ymax=(ulimith2-ylims[0])/(ylims[1]-ylims[0]),color='b',alpha=4.)

            print(' \n')

        mixstyle = 'full'
        try:
            self.get('mix_qtop_1')
        except:
            mixstyle = 'twozone'

        old_percent = 0

        if mixstyle == 'full':
            for i in range(len(x)):
            # writing reading status
                percent = int(i*100/(len(x) - 1))
                if percent >= old_percent + 5:
                    sys.stdout.flush()
                    sys.stdout.write("\r creating color map2 " + "...%d%%" % percent)
                    old_percent = percent

                for j in range(1,mix_zones+1):
                    ulimit=self.get('mix_qtop_'+str(j))[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                    if j==1:
                        llimit=0.0
                    else:
                        llimit=self.get('mix_qtop_'+str(j-1))[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                    mtype=self.get('mix_type_'+str(j))[modstart:modstop][i*dx]
                    if llimit!=ulimit:
                        if mtype == 1:
                            ax.axvline(xxx[i*dx],ymin=old_div((llimit-ylims[0]),(ylims[1]-ylims[0])),ymax=old_div((ulimit-ylims[0]),(ylims[1]-ylims[0])),color='k',alpha=3., linewidth=.5)

            print(' \n')

        old_percent = 0

        if mixstyle == 'twozone':
            for i in range(len(x)):
            # writing reading status
                percent = int(i*100/(len(x) - 1))
                if percent >= old_percent + 5:
                    sys.stdout.flush()
                    sys.stdout.write("\r creating color map2 " + "...%d%%" % percent)
                    old_percent = percent

                ulimit=self.get('conv_mx1_top')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                llimit=self.get('conv_mx1_bot')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                if llimit!=ulimit:
                    ax.axvline(xxx[i*dx],ymin=old_div((llimit-ylims[0]),(ylims[1]-ylims[0])),ymax=old_div((ulimit-ylims[0]),(ylims[1]-ylims[0])),color='k',alpha=5.,linewidth=.5)
                ulimit=self.get('conv_mx2_top')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                llimit=self.get('conv_mx2_bot')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                if llimit!=ulimit:
                    ax.axvline(xxx[i*dx],ymin=old_div((llimit-ylims[0]),(ylims[1]-ylims[0])),ymax=old_div((ulimit-ylims[0]),(ylims[1]-ylims[0])),color='k',alpha=3.,linewidth=.5)

            print(' \n')

        print('engenstyle was ', engenstyle)
        print('mixstyle was ', mixstyle)
        print('\n finished preparing color map')

        #fig.savefig(outfile)
        pl.show()

    def kip_cont(self, ifig=110, modstart=0, modstop=-1,t0_model=0,
                 outfile='out.png', xlims=[0.,0.], ylims=[0.,0.],
                 xres=1000, yres=1000, ixaxis='model_number',
                 mix_zones=20, burn_zones=20, plot_radius=False,
                 engenPlus=True, engenMinus=False,
                 landscape_plot=False, rad_lines=False, profiles=[],
                 showfig=True, outlines=True, boundaries=True,
                 c12_boundary=False, rasterise=False, yscale='1.',
                 engenlevels=None,CBM=False):
        """
        This function creates a Kippenhahn plot with energy flux using
        contours.

        This plot uses mixing_regions and burning_regions written to
        your history.data or star.log. Set both variables in the
        log_columns.list file to 20 as a start.

        The output log file should then contain columns called
        "mix_type_n", "mix_qtop_n", "burn_type_n" and "burn_qtop_n".
        The number of columns (i.e. the biggest value of n) is what
        goes in the arguments as mix_zones and burn_zones.  DO NOT
        WORRY! if you do not have these columns, just leave the default
        values alone and the script should recognise that you do not
        have these columns and make the most detailed plot that is
        available to you.

        Defaults are set to get some plot, that may not look great if
        you zoom in interactively. Play with xres and yres as well as
        setting the xlims to ylims to the region you are interested in.

        Parameters
        ----------
        ifig : integer, optional
            Figure frame number.  The default is 110.
        modstart : integer, optional
            Model from which you want to plot (be careful if your
            history.data or star.log output is sparse...).  If it is
            0 then it starts from the beginning, works even if
            log_cnt > 1.
            The default is 0.
        modstop : integer, optional
            Model to which you wish to plot, -1 corresponds to end
            [if log_cnt>1, devide modstart and modstop by log_cnt,
            this needs to be improved! SJ: this should be ficed now].
            The defalut is -1.
        t0_model : integer, optional
            Model number from which to reset the time to 0. Typically,
            if modstart!=0, t0_model=modstart is a good choice, but we
            leave the choice to the user in case the time is wished to
            start from 0 at a different key point of the evolution.
            The default value is 0.
        outfile : sting, optional
            'filename + extension' where you want to save the figure.
            The defalut is "out.png".
        xlims, ylims : list, optional
            Plot limits, however these are somewhat obsolete now that
            we have modstart and modstop.  Leaving them as 0. is
            probably no slower, and you can always zoom in afterwards
            in mpl.  ylims is important for well resolved thermal
            pulse etc plots; it's best to get the upper and lower limits
            of he-intershell using s.kippenhahn_CO(1,'model') first.
            The default is [0., 0.].
        xres, yres : integer, optional
            plot resolution. Needless to say that increasing these
            values will yield a nicer plot with some slow-down in
            plotting time.  You will most commonly change xres.  For a
            prelim plot, try xres~200, then bump it up to anywhere from
            1000-10000 for real nicely resolved, publication quality
            plots.  The default is 1000.
        ixaxis : string, optional
            Either 'log_time_left', 'age', or 'model_number'.  The
            default is "model_number".
        mix_zones, burn_zones : integer, optional
            As described above, if you have more detailed output about
            your convection and energy generation boundaries in columns
            mix_type_n, mix_qtop_n, burn_type_n and burn_qtop_n, you
            need to specify the total number of columns for mixing zones
            and burning zones that you have.  Can't work this out from
            your history.data or star.log file?  Check the
            history_columns.list that you used, it'll be the number
            after "mixing regions" and "burning regions".  Can't see
            these columns?  leave it and 2 conv zones and 2 burn zones
            will be drawn using other data that you certainly should
            have in your history.data or star.log file.  The defalut for
            both is 20.
        plot_radius : boolean, optional
            Whether on a second y-axis you want to plot the radius of
            the surface and the he-free core.  The default is False.
        engenPlus : boolean
            Plot energy generation contours for eps_nuc>0.  The default
            is True.
        endgenMinus : boolean, optional
            Plot energy generation contours for eos_nuc<0.  The default
            is True.
        landscape_plot : boolean, optionla
            The default is False.
        rad_lines : boolean, optional
            The deafault is False.
        profiles : list, optional
            The default is [].
        showfig : boolean, optional
            The default is True.
        outlines : boolean, optional
            Whether or not to plot outlines of conv zones in darker
            colour.
        boundaries : boolean, optional
            Whether or not to plot H-, He- and C-free boundaries.
        c12_boundary : boolean, optional
            The default is False.
        rasterise : boolean, optional
            Whether or not to rasterise the contour regions to make
            smaller vector graphics figures.  The default is False.
        yscale : string, optional
            Re-scale the y-axis by this amount
        engenlevels : list
            Give cusstom levels to the engenPlus contour. If None,
            the levels are chosen automatically.
            The default is None.
        CBM : boolean, optional
            plot contours for where CBM is active?

        Notes
        -----
        The parameter xlims is depricated.

        """

        # Find correct modstart and modstop:
        mod=np.array([int(i) for i in self.get('model_number')])
        mod1=np.abs(mod-modstart).argmin()
        mod2=np.abs(mod-modstop).argmin()
        if modstart != 0 : modstart=mod1
        if modstop != -1 : modstop=mod2

        xxyy=[self.get('star_age')[modstart:modstop],self.get('star_age')[modstart:modstop]]
        mup = max(float(self.get('star_mass')[0])*1.02,1.0)
        nmodels=len(self.get('model_number')[modstart:modstop])

        if ylims == [0.,0.]:
            mup   = max(float(self.get('star_mass')[0])*1.02,1.0)
            mDOWN = 0.
        else:
            mup = ylims[1]
            mDOWN = ylims[0]

        # y-axis resolution
        ny=yres
        #dy=mup/float(ny)
        dy = old_div((mup-mDOWN),float(ny))

        # x-axis resolution
        maxpoints=xres
        dx=int(max(1,old_div(nmodels,maxpoints)))

        #y = np.arange(0., mup, dy)
        y = np.arange(mDOWN, mup, dy)
        x = np.arange(0, nmodels, dx)
        Msol=1.98892E+33

        engenstyle = 'full'

        B1=np.zeros([len(y),len(x)],float)
        B2=np.zeros([len(y),len(x)],float)
        try:
            self.get('burn_qtop_1')
        except:
            engenstyle = 'twozone'
        if engenstyle == 'full' and (engenPlus == True or engenMinus == True):
            ulimit_array = np.array([self.get('burn_qtop_'+str(j))[modstart:modstop:dx]*self.get('star_mass')[modstart:modstop:dx] for j in range(1,burn_zones+1)])
            #ulimit_array = np.around(ulimit_array,decimals=len(str(dy))-2)
            llimit_array = np.delete(ulimit_array,-1,0)
            llimit_array = np.insert(ulimit_array,0,0.,0)
            #llimit_array = np.around(llimit_array,decimals=len(str(dy))-2)
            btype_array = np.array([self.get('burn_type_'+str(j))[modstart:modstop:dx] for j in range(1,burn_zones+1)])
            old_percent = 0
            for i in range(len(x)):
                # writing status
                percent = int(i*100/(len(x) - 1))
                if percent >= old_percent + 5:
                    sys.stdout.flush()
                    sys.stdout.write("\r creating color map burn " + "...%d%%" % percent)
                    old_percent = percent

                for j in range(burn_zones):
                    if btype_array[j,i] > 0. and abs(btype_array[j,i]) < 99.:
                        B1[(np.abs(y-llimit_array[j][i])).argmin():(np.abs(y-ulimit_array[j][i])).argmin()+1,i] = 10.0**(btype_array[j,i])
                    elif btype_array[j,i] < 0. and abs(btype_array[j,i]) < 99.:
                        B2[(np.abs(y-llimit_array[j][i])).argmin():(np.abs(y-ulimit_array[j][i])).argmin()+1,i] = 10.0**(abs(btype_array[j,i]))
            print(' \n')

        if engenstyle == 'twozone' and (engenPlus == True or engenMinus == True):
            V=np.zeros([len(y),len(x)],float)
            old_percent = 0
            for i in range(len(x)):
                # writing status
                percent = int(i*100/(len(x) - 1))
                if percent >= old_percent + 5:
                    sys.stdout.flush()
                    sys.stdout.write("\r creating color map1 " + "...%d%%" % percent)
                    old_percent = percent
                llimitl1=old_div(self.get('epsnuc_M_1')[modstart:modstop][i*dx],Msol)
                ulimitl1=old_div(self.get('epsnuc_M_4')[modstart:modstop][i*dx],Msol)
                llimitl2=old_div(self.get('epsnuc_M_5')[modstart:modstop][i*dx],Msol)
                ulimitl2=old_div(self.get('epsnuc_M_8')[modstart:modstop][i*dx],Msol)
                llimith1=old_div(self.get('epsnuc_M_2')[modstart:modstop][i*dx],Msol)
                ulimith1=old_div(self.get('epsnuc_M_3')[modstart:modstop][i*dx],Msol)
                llimith2=old_div(self.get('epsnuc_M_6')[modstart:modstop][i*dx],Msol)
                ulimith2=old_div(self.get('epsnuc_M_7')[modstart:modstop][i*dx],Msol)
                # lower thresh first, then upper thresh:
                if llimitl1!=ulimitl1:
                    for k in range(ny):
                        if llimitl1<=y[k] and ulimitl1>y[k]:
                            V[k,i]=10.
                if llimitl2!=ulimitl2:
                    for k in range(ny):
                        if llimitl2<=y[k] and ulimitl2>y[k]:
                            V[k,i]=10.
                if llimith1!=ulimith1:
                    for k in range(ny):
                        if llimith1<=y[k] and ulimith1>y[k]:
                            V[k,i]=30.
                if llimith2!=ulimith2:
                    for k in range(ny):
                        if llimith2<=y[k] and ulimith2>y[k]:
                            V[k,i]=30.
            print(' \n')

        mixstyle = 'full'
        try:
            self.get('mix_qtop_1')
        except:
            mixstyle = 'twozone'
        if mixstyle == 'full':
            old_percent = 0
            Z=np.zeros([len(y),len(x)],float)
            if CBM:
                Zcbm=np.zeros([len(y),len(x)],float)
            ulimit_array = np.array([self.get('mix_qtop_'+str(j))[modstart:modstop:dx]*self.get('star_mass')[modstart:modstop:dx] for j in range(1,mix_zones+1)])
            llimit_array = np.delete(ulimit_array,-1,0)
            llimit_array = np.insert(ulimit_array,0,0.,0)
            mtype_array = np.array([self.get('mix_type_'+str(j))[modstart:modstop:dx] for j in range(1,mix_zones+1)])
            for i in range(len(x)):
                # writing status
                percent = int(i*100/(len(x) - 1))
                if percent >= old_percent + 5:
                    sys.stdout.flush()
                    sys.stdout.write("\r creating color map mix " + "...%d%%" % percent)
                    old_percent = percent
                for j in range(mix_zones):
                    if mtype_array[j,i] == 1.:
                        Z[(np.abs(y-llimit_array[j][i])).argmin():(np.abs(y-ulimit_array[j][i])).argmin()+1,i] = 1.
                    if CBM:
                        if mtype_array[j,i] == 2.:
                            Zcbm[(np.abs(y-llimit_array[j][i])).argmin():(np.abs(y-ulimit_array[j][i])).argmin()+1,i] = 1.
            print(' \n')

        if mixstyle == 'twozone':
            Z=np.zeros([len(y),len(x)],float)
            old_percent = 0
            for i in range(len(x)):
            # writing reading status
                # writing status
                percent = int(i*100/(len(x) - 1))
                if percent >= old_percent + 5:
                    sys.stdout.flush()
                    sys.stdout.write("\r creating color map mix " + "...%d%%" % percent)
                    old_percent = percent

                ulimit=self.get('conv_mx1_top')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                llimit=self.get('conv_mx1_bot')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                if llimit!=ulimit:
                    for k in range(ny):
                        if llimit<=y[k] and ulimit>y[k]:
                            Z[k,i]=1.
                ulimit=self.get('conv_mx2_top')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                llimit=self.get('conv_mx2_bot')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                if llimit!=ulimit:
                    for k in range(ny):
                        if llimit<=y[k] and ulimit>y[k]:
                            Z[k,i]=1.
            print(' \n')

        if rad_lines == True:
            masses = np.arange(0.1,1.5,0.1)
            rads=[[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
            modno=[]
            for i in range(len(profiles)):
                p=mesa_profile('./LOGS',profiles[i])
                modno.append(p.header_attr['model_number'])
                for j in range(len(masses)):
                    idx=np.abs(p.get('mass')-masses[j]).argmin()
                    rads[j].append(p.get('radius')[idx])

        print('engenstyle was ', engenstyle)
        print('mixstyle was ', mixstyle)
        print('\n finished preparing color map')

        ########################################################################
        #----------------------------------plot--------------------------------#
        fig = pyl.figure(ifig)
#        fsize=20
        if landscape_plot == True:
            fig.set_size_inches(9,4)
            pl.gcf().subplots_adjust(bottom=0.2)
            pl.gcf().subplots_adjust(right=0.85)

#        params = {'axes.labelsize':  fsize,
#          'text.fontsize':   fsize,
#          'legend.fontsize': fsize,
#          'xtick.labelsize': fsize*0.8,
#          'ytick.labelsize': fsize*0.8,
#          'text.usetex': False}
#        pyl.rcParams.update(params)

        #ax=pl.axes([0.1,0.1,0.9,0.8])

        #fig=pl.figure()
        ax=pl.axes()

        if ixaxis == 'log_time_left':
        # log of time left until core collapse
            gage= self.get('star_age')
            lage=np.zeros(len(gage))
            agemin = max(old_div(abs(gage[-1]-gage[-2]),5.),1.e-10)
            for i in np.arange(len(gage)):
                if gage[-1]-gage[i]>agemin:
                    lage[i]=np.log10(gage[-1]-gage[i]+agemin)
                else :
                    lage[i]=np.log10(agemin)
            xxx = lage[modstart:modstop]
            print('plot versus time left')
            ax.set_xlabel('$\mathrm{log}_{10}(t^*) \, \mathrm{(yr)}$') #,fontsize=fsize)
            if xlims[1] == 0.:
                xlims = [xxx[0],xxx[-1]]
        elif ixaxis =='model_number':
            xxx= self.get('model_number')[modstart:modstop]
            print('plot versus model number')
            ax.set_xlabel('Model number') # ,fontsize=fsize)
            if xlims[1] == 0.:
                xlims = [self.get('model_number')[modstart],self.get('model_number')[modstop]]
        elif ixaxis =='age':
            if t0_model != 0:
                t0_mod=np.abs(mod-t0_model).argmin()
                xxx= self.get('star_age')[modstart:modstop] - self.get('star_age')[t0_mod]
                print('plot versus age')
                ax.set_xlabel('Age [yr] - '+str(self.get('star_age')[modstart])) #,fontsize=fsize)
            else:
                xxx= old_div(self.get('star_age')[modstart:modstop],1.e6)
                ax.set_xlabel('Age [Myr]') #,fontsize=fsize)
            if xlims[1] == 0.:
                xlims = [xxx[0],xxx[-1]]

        ax.set_ylabel('$\mathrm{Mass }(M_\odot)$')

        # some stuff for rasterizing only the contour part of the plot, for nice, but light, eps:
        class ListCollection(Collection):
            def __init__(self, collections, **kwargs):
                Collection.__init__(self, **kwargs)
                self.set_collections(collections)
            def set_collections(self, collections):
                self._collections = collections
            def get_collections(self):
                return self._collections
            @allow_rasterization
            def draw(self, renderer):
                for _c in self._collections:
                    _c.draw(renderer)

        def insert_rasterized_contour_plot(c):
            collections = c.collections
            for _c in collections:
                _c.remove()
            cc = ListCollection(collections, rasterized=True)
            ax = pl.gca()
            ax.add_artist(cc)
            return cc

        cmapMIX = matplotlib.colors.ListedColormap(['w','#8B8386']) # rose grey
        if CBM:
            cmapCBM = matplotlib.colors.ListedColormap(['w','g']) # green
        cmapB1  = pyl.cm.get_cmap('Blues')
        cmapB2  = pl.cm.get_cmap('Reds')

        ylims1=[0.,0.]
        ylims1[0]=ylims[0]
        ylims1[1]=ylims[1]
        if ylims == [0.,0.]:
            ylims[0] = 0.
            ylims[1] = mup
        if ylims[0] != 0.:
            ylab='$(\mathrm{Mass }$ - '+str(ylims[0])
            if yscale!='1.':
                ylab+=') / '+yscale+' $M_\odot$'
            else:
                ylab+=') / $M_\odot$'
            ax.set_ylabel(ylab)
            y = y - ylims[0]
            y = y*float(yscale) # SJONES tweak
            ylims[0] = y[0]
            ylims[1] = y[-1]

        print('plotting contours')
        CMIX    = ax.contourf(xxx[::dx],y,Z, cmap=cmapMIX,alpha=0.6,levels=[0.5,1.5])
        #CMIX    = ax.pcolor(xxx[::dx],y,Z, cmap=cmapMIX,alpha=0.6,vmin=0.5,vmax=1.5)
        if rasterise==True:
            insert_rasterized_contour_plot(CMIX)
        if outlines == True:
            CMIX_outlines    = ax.contour(xxx[::dx],y,Z, cmap=cmapMIX)
            if rasterise==True:
                insert_rasterized_contour_plot(CMIX_outlines)

        if CBM:
            CCBM    = ax.contourf(xxx[::dx],y,Zcbm, cmap=cmapCBM,alpha=0.6,levels=[0.5,1.5])
            if rasterise==True:
                insert_rasterized_contour_plot(CCBM)
            if outlines == True:
                CCBM_outlines    = ax.contour(xxx[::dx],y,Zcbm, cmap=cmapCBM)
                if rasterise==True:
                    insert_rasterized_contour_plot(CCBM_outlines)

        if engenstyle == 'full' and engenPlus == True:
            if engenlevels!= None:
                CBURN1  = ax.contourf(xxx[::dx],y,B1, cmap=cmapB1, alpha=0.5, locator=matplotlib.ticker.LogLocator(),levels=engenlevels)
                if outlines:
                    CB1_outlines  = ax.contour(xxx[::dx],y,B1, cmap=cmapB1, alpha=0.7, locator=matplotlib.ticker.LogLocator(),levels=engenlevels)
            else:
                CBURN1  = ax.contourf(xxx[::dx],y,B1, cmap=cmapB1, alpha=0.5, locator=matplotlib.ticker.LogLocator())
                if outlines:
                    CB1_outlines  = ax.contour(xxx[::dx],y,B1, cmap=cmapB1, alpha=0.7, locator=matplotlib.ticker.LogLocator())
            CBARBURN1 = pyl.colorbar(CBURN1)
            CBARBURN1.set_label('$|\epsilon_\mathrm{nuc}-\epsilon_{\\nu}| \; (\mathrm{erg\,g}^{-1}\mathrm{\,s}^{-1})$') #,fontsize=fsize)
            if rasterise==True:
                insert_rasterized_contour_plot(CBURN1)
                if outlines:
                    insert_rasterized_contour_plot(CB1_outlines)

        if engenstyle == 'full' and engenMinus == True:
            CBURN2  = ax.contourf(xxx[::dx],y,B2, cmap=cmapB2, alpha=0.5, locator=matplotlib.ticker.LogLocator())
            if outlines:
                CBURN2_outlines  = ax.contour(xxx[::dx],y,B2, cmap=cmapB2, alpha=0.7, locator=matplotlib.ticker.LogLocator())
            CBARBURN2 = pl.colorbar(CBURN2)
            if engenPlus == False:
                CBARBURN2.set_label('$|\epsilon_\mathrm{nuc}-\epsilon_{\\nu}| \; (\mathrm{erg\,g}^{-1}\mathrm{\,s}^{-1})$') #,fontsize=fsize)
            if rasterise==True:
                insert_rasterized_contour_plot(CBURN2)
                if outlines:
                    insert_rasterized_contour_plot(CB2_outlines)

        if engenstyle == 'twozone' and (engenPlus == True or engenMinus == True):
            ax.contourf(xxx[::dx],y,V, cmap=cmapB1, alpha=0.5)

        print('plotting patches')
        mtot=self.get('star_mass')[modstart:modstop][::dx]
        mtot1=(mtot-ylims1[0])*float(yscale)
        ax.plot(xxx[::dx],mtot1,'k-')

        if boundaries == True:
            print('plotting abund boundaries')
            try:
                bound=self.get('h1_boundary_mass')[modstart:modstop]
                bound1=(bound-ylims1[0])*float(yscale)
                ax.plot(xxx,bound1,label='H boundary',linestyle='-')

                bound=self.get('he4_boundary_mass')[modstart:modstop]
                bound1=(bound-ylims1[0])*float(yscale)
                ax.plot(xxx,bound1,label='He boundary',linestyle='--')

                bound=self.get('c12_boundary_mass')[modstart:modstop]
                bound1=(bound-ylims1[0])*float(yscale)
                ax.plot(xxx,bound1,label='C boundary',linestyle='-.')

            except:
                try:
                    bound=self.get('he_core_mass')[modstart:modstop]
                    bound1=(bound-ylims1[0])*float(yscale)
                    ax.plot(xxx,bound1,label='H boundary',linestyle='-')

                    bound=self.get('c_core_mass')[modstart:modstop]-ylims[0]
                    bound1=(bound-ylims1[0])*float(yscale)
                    ax.plot(xxx,bound1,label='He boundary',linestyle='--')

                    bound=self.get('o_core_mass')[modstart:modstop]-ylims[0]
                    bound1=(bound-ylims1[0])*float(yscale)
                    ax.plot(xxx,bound1,label='C boundary',linestyle='-.')

                    bound=self.get('si_core_mass')[modstart:modstop]-ylims[0]
                    bound1=(bound-ylims1[0])*float(yscale)
                    ax.plot(xxx,bound1,label='C boundary',linestyle='-.')

                    bound=self.get('fe_core_mass')[modstart:modstop]-ylims[0]
                    bound1=(bound-ylims1[0])*float(yscale)
                    ax.plot(xxx,bound1,label='C boundary',linestyle='-.')

                except:
#                    print 'problem to plot boundaries for this plot'
                    pass

        ax.axis([xlims[0],xlims[1],ylims[0],ylims[1]])

        if plot_radius == True:
            ax2=pyl.twinx()
            ax2.plot(xxx,np.log10(self.get('he4_boundary_radius')[modstart:modstop]),label='He boundary radius',color='k',linewidth=1.,linestyle='-.')
            ax2.plot(xxx,self.get('log_R')[modstart:modstop],label='radius',color='k',linewidth=1.,linestyle='-.')
            ax2.set_ylabel('log(radius)')
        if rad_lines == True:
            ax2=pyl.twinx()
            for i in range(len(masses)):
                ax2.plot(modno,np.log10(rads[i]),color='k')

        if outfile[-3:]=='png':
            fig.savefig(outfile,dpi=300)
        elif outfile[-3:]=='eps':
            fig.savefig(outfile,format='eps')
        elif outfile[-3:]=='pdf':
            fig.savefig(outfile,format='pdf')
        if showfig == True:
            pyl.show()
# we may or may not need this below
#        fig.clear()

    def kip_cont_COmdot(self, ifig=110, modstart=0, modstop=-1,t0_model=0,
                 outfile='out.png', xlims=[0.,0.], ylims=[0.,0.],
                 xres=1000, yres=1000, ixaxis='model_number',
                 mix_zones=20, burn_zones=20, plot_radius=False,
                 engenPlus=True, engenMinus=False,
                 landscape_plot=False, rad_lines=False, profiles=[],
                 showfig=True, outlines=True, boundaries=True,
                 c12_boundary=False, rasterise=False, yscale='1.',
                 engenlevels=None,CO_ratio=True):
        """
        This function creates a Kippenhahn plot with energy flux using
        contours.

        This plot uses mixing_regions and burning_regions written to
        your history.data or star.log. Set both variables in the
        log_columns.list file to 20 as a start.

        The output log file should then contain columns called
        "mix_type_n", "mix_qtop_n", "burn_type_n" and "burn_qtop_n".
        The number of columns (i.e. the biggest value of n) is what
        goes in the arguments as mix_zones and burn_zones.  DO NOT
        WORRY! if you do not have these columns, just leave the default
        values alone and the script should recognise that you do not
        have these columns and make the most detailed plot that is
        available to you.

        Defaults are set to get some plot, that may not look great if
        you zoom in interactively. Play with xres and yres as well as
        setting the xlims to ylims to the region you are interested in.

        Parameters
        ----------
        ifig : integer, optional
            Figure frame number.  The default is 110.
        modstart : integer, optional
            Model from which you want to plot (be careful if your
            history.data or star.log output is sparse...).  If it is
            0 then it starts from the beginning, works even if
            log_cnt > 1.
            The default is 0.
        modstop : integer, optional
            Model to which you wish to plot, -1 corresponds to end
            [if log_cnt>1, devide modstart and modstop by log_cnt,
            this needs to be improved! SJ: this should be ficed now].
            The defalut is -1.
        t0_model : integer, optional
            Model number from which to reset the time to 0. Typically,
            if modstart!=0, t0_model=modstart is a good choice, but we
            leave the choice to the user in case the time is wished to
            start from 0 at a different key point of the evolution.
            The default value is 0.
        outfile : sting, optional
            'filename + extension' where you want to save the figure.
            The defalut is "out.png".
        xlims, ylims : list, optional
            Plot limits, however these are somewhat obsolete now that
            we have modstart and modstop.  Leaving them as 0. is
            probably no slower, and you can always zoom in afterwards
            in mpl.  ylims is important for well resolved thermal
            pulse etc plots; it's best to get the upper and lower limits
            of he-intershell using s.kippenhahn_CO(1,'model') first.
            The default is [0., 0.].
        xres, yres : integer, optional
            plot resolution. Needless to say that increasing these
            values will yield a nicer plot with some slow-down in
            plotting time.  You will most commonly change xres.  For a
            prelim plot, try xres~200, then bump it up to anywhere from
            1000-10000 for real nicely resolved, publication quality
            plots.  The default is 1000.
        ixaxis : string, optional
            Either 'log_time_left', 'age', or 'model_number'.  The
            default is "model_number".
        mix_zones, burn_zones : integer, optional
            As described above, if you have more detailed output about
            your convection and energy generation boundaries in columns
            mix_type_n, mix_qtop_n, burn_type_n and burn_qtop_n, you
            need to specify the total number of columns for mixing zones
            and burning zones that you have.  Can't work this out from
            your history.data or star.log file?  Check the
            history_columns.list that you used, it'll be the number
            after "mixing regions" and "burning regions".  Can't see
            these columns?  leave it and 2 conv zones and 2 burn zones
            will be drawn using other data that you certainly should
            have in your history.data or star.log file.  The defalut for
            both is 20.
        plot_radius : boolean, optional
            Whether on a second y-axis you want to plot the radius of
            the surface and the he-free core.  The default is False.
        engenPlus : boolean
            Plot energy generation contours for eps_nuc>0.  The default
            is True.
        endgenMinus : boolean, optional
            Plot energy generation contours for eos_nuc<0.  The default
            is True.
        landscape_plot : boolean, optionla
            The default is False.
        rad_lines : boolean, optional
            The deafault is False.
        profiles : list, optional
            The default is [].
        showfig : boolean, optional
            The default is True.
        outlines : boolean, optional
            Whether or not to plot outlines of conv zones in darker
            colour.
        boundaries : boolean, optional
            Whether or not to plot H-, He- and C-free boundaries.
        c12_boundary : boolean, optional
            The default is False.
        rasterise : boolean, optional
            Whether or not to rasterise the contour regions to make
            smaller vector graphics figures.  The default is False.
        yscale : string, optional
            Re-scale the y-axis by this amount
        engenlevels : list
            Give cusstom levels to the engenPlus contour. If None,
            the levels are chosen automatically.
            The default is None.

        Notes
        -----
        The parameter xlims is depricated.

        """

        # Find correct modstart and modstop:
        mod=np.array([int(i) for i in self.get('model_number')])
        mod1=np.abs(mod-modstart).argmin()
        mod2=np.abs(mod-modstop).argmin()
        if modstart != 0 : modstart=mod1
        if modstop != -1 : modstop=mod2

        xxyy=[self.get('star_age')[modstart:modstop],self.get('star_age')[modstart:modstop]]
        mup = max(float(self.get('star_mass')[0])*1.02,1.0)
        nmodels=len(self.get('model_number')[modstart:modstop])

        if ylims == [0.,0.]:
            mup   = max(float(self.get('star_mass')[0])*1.02,1.0)
            mDOWN = 0.
        else:
            mup = ylims[1]
            mDOWN = ylims[0]

        # y-axis resolution
        ny=yres
        #dy=mup/float(ny)
        dy = old_div((mup-mDOWN),float(ny))

        # x-axis resolution
        maxpoints=xres
        dx=int(max(1,old_div(nmodels,maxpoints)))

        #y = np.arange(0., mup, dy)
        y = np.arange(mDOWN, mup, dy)
        x = np.arange(0, nmodels, dx)
        Msol=1.98892E+33

        engenstyle = 'full'

        B1=np.zeros([len(y),len(x)],float)
        B2=np.zeros([len(y),len(x)],float)
        try:
            self.get('burn_qtop_1')
        except:
            engenstyle = 'twozone'
        if engenstyle == 'full' and (engenPlus == True or engenMinus == True):
            ulimit_array = np.array([self.get('burn_qtop_'+str(j))[modstart:modstop:dx]*self.get('star_mass')[modstart:modstop:dx] for j in range(1,burn_zones+1)])
            #ulimit_array = np.around(ulimit_array,decimals=len(str(dy))-2)
            llimit_array = np.delete(ulimit_array,-1,0)
            llimit_array = np.insert(ulimit_array,0,0.,0)
            #llimit_array = np.around(llimit_array,decimals=len(str(dy))-2)
            btype_array = np.array([self.get('burn_type_'+str(j))[modstart:modstop:dx] for j in range(1,burn_zones+1)])
            old_percent = 0
            for i in range(len(x)):
                # writing status
                percent = int(i*100/(len(x) - 1))
                if percent >= old_percent + 5:
                    sys.stdout.flush()
                    sys.stdout.write("\r creating color map burn " + "...%d%%" % percent)
                    old_percent = percent

                for j in range(burn_zones):
                    if btype_array[j,i] > 0. and abs(btype_array[j,i]) < 99.:
                        B1[(np.abs(y-llimit_array[j][i])).argmin():(np.abs(y-ulimit_array[j][i])).argmin()+1,i] = 10.0**(btype_array[j,i])
                    elif btype_array[j,i] < 0. and abs(btype_array[j,i]) < 99.:
                        B2[(np.abs(y-llimit_array[j][i])).argmin():(np.abs(y-ulimit_array[j][i])).argmin()+1,i] = 10.0**(abs(btype_array[j,i]))
            print(' \n')

        if engenstyle == 'twozone' and (engenPlus == True or engenMinus == True):
            V=np.zeros([len(y),len(x)],float)
            old_percent = 0
            for i in range(len(x)):
                # writing status
                percent = int(i*100/(len(x) - 1))
                if percent >= old_percent + 5:
                    sys.stdout.flush()
                    sys.stdout.write("\r creating color map1 " + "...%d%%" % percent)
                    old_percent = percent
                llimitl1=old_div(self.get('epsnuc_M_1')[modstart:modstop][i*dx],Msol)
                ulimitl1=old_div(self.get('epsnuc_M_4')[modstart:modstop][i*dx],Msol)
                llimitl2=old_div(self.get('epsnuc_M_5')[modstart:modstop][i*dx],Msol)
                ulimitl2=old_div(self.get('epsnuc_M_8')[modstart:modstop][i*dx],Msol)
                llimith1=old_div(self.get('epsnuc_M_2')[modstart:modstop][i*dx],Msol)
                ulimith1=old_div(self.get('epsnuc_M_3')[modstart:modstop][i*dx],Msol)
                llimith2=old_div(self.get('epsnuc_M_6')[modstart:modstop][i*dx],Msol)
                ulimith2=old_div(self.get('epsnuc_M_7')[modstart:modstop][i*dx],Msol)
                # lower thresh first, then upper thresh:
                if llimitl1!=ulimitl1:
                    for k in range(ny):
                        if llimitl1<=y[k] and ulimitl1>y[k]:
                            V[k,i]=10.
                if llimitl2!=ulimitl2:
                    for k in range(ny):
                        if llimitl2<=y[k] and ulimitl2>y[k]:
                            V[k,i]=10.
                if llimith1!=ulimith1:
                    for k in range(ny):
                        if llimith1<=y[k] and ulimith1>y[k]:
                            V[k,i]=30.
                if llimith2!=ulimith2:
                    for k in range(ny):
                        if llimith2<=y[k] and ulimith2>y[k]:
                            V[k,i]=30.
            print(' \n')

        mixstyle = 'full'
        try:
            self.get('mix_qtop_1')
        except:
            mixstyle = 'twozone'
        if mixstyle == 'full':
            old_percent = 0
            Z=np.zeros([len(y),len(x)],float)
            ulimit_array = np.array([self.get('mix_qtop_'+str(j))[modstart:modstop:dx]*self.get('star_mass')[modstart:modstop:dx] for j in range(1,mix_zones+1)])
            llimit_array = np.delete(ulimit_array,-1,0)
            llimit_array = np.insert(ulimit_array,0,0.,0)
            mtype_array = np.array([self.get('mix_type_'+str(j))[modstart:modstop:dx] for j in range(1,mix_zones+1)])
            for i in range(len(x)):
                # writing status
                percent = int(i*100/(len(x) - 1))
                if percent >= old_percent + 5:
                    sys.stdout.flush()
                    sys.stdout.write("\r creating color map mix " + "...%d%%" % percent)
                    old_percent = percent
                for j in range(mix_zones):
                    if mtype_array[j,i] == 1.:
                        Z[(np.abs(y-llimit_array[j][i])).argmin():(np.abs(y-ulimit_array[j][i])).argmin()+1,i] = 1.
            print(' \n')

        if mixstyle == 'twozone':
            Z=np.zeros([len(y),len(x)],float)
            old_percent = 0
            for i in range(len(x)):
            # writing reading status
                # writing status
                percent = int(i*100/(len(x) - 1))
                if percent >= old_percent + 5:
                    sys.stdout.flush()
                    sys.stdout.write("\r creating color map mix " + "...%d%%" % percent)
                    old_percent = percent

                ulimit=self.get('conv_mx1_top')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                llimit=self.get('conv_mx1_bot')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                if llimit!=ulimit:
                    for k in range(ny):
                        if llimit<=y[k] and ulimit>y[k]:
                            Z[k,i]=1.
                ulimit=self.get('conv_mx2_top')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                llimit=self.get('conv_mx2_bot')[modstart:modstop][i*dx]*self.get('star_mass')[modstart:modstop][i*dx]
                if llimit!=ulimit:
                    for k in range(ny):
                        if llimit<=y[k] and ulimit>y[k]:
                            Z[k,i]=1.
            print(' \n')

        if rad_lines == True:
            masses = np.arange(0.1,1.5,0.1)
            rads=[[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
            modno=[]
            for i in range(len(profiles)):
                p=mesa_profile('./LOGS',profiles[i])
                modno.append(p.header_attr['model_number'])
                for j in range(len(masses)):
                    idx=np.abs(p.get('mass')-masses[j]).argmin()
                    rads[j].append(p.get('radius')[idx])

        print('engenstyle was ', engenstyle)
        print('mixstyle was ', mixstyle)
        print('\n finished preparing color map')

        ########################################################################
        #----------------------------------plot--------------------------------#
        fig = pyl.figure(ifig)
        fsize=20
        if landscape_plot == True:
            fig.set_size_inches(9,4)
            pl.gcf().subplots_adjust(bottom=0.2)
            pl.gcf().subplots_adjust(right=0.85)

#        params = {'axes.labelsize':  fsize,
#          'text.fontsize':   fsize,
#          'legend.fontsize': fsize,
#          'xtick.labelsize': fsize*0.8,
#          'ytick.labelsize': fsize*0.8,
#          'text.usetex': False}
#        pyl.rcParams.update(params)

        #ax=pl.axes([0.1,0.1,0.9,0.8])

        #fig=pl.figure()

        ax1 = pyl.axes(frameon=False)
        ax1.axes.get_xaxis().set_visible(False)
        ax1.axes.get_yaxis().set_visible(False)
        pyl.subplots_adjust(hspace=0.7)

        #axk=fig.add_subplot(211)
        ax = pyl.subplot2grid((3,3), (0, 0), rowspan=2, colspan=3)
        #ax=pl.axes()

        if ixaxis == 'log_time_left':
        # log of time left until core collapse
            gage= self.get('star_age')
            lage=np.zeros(len(gage))
            agemin = max(old_div(abs(gage[-1]-gage[-2]),5.),1.e-10)
            for i in np.arange(len(gage)):
                if gage[-1]-gage[i]>agemin:
                    lage[i]=np.log10(gage[-1]-gage[i]+agemin)
                else :
                    lage[i]=np.log10(agemin)
            xxx = lage[modstart:modstop]
            print('plot versus time left')
            ax.set_xlabel('$\mathrm{log}_{10}(t^*) \, \mathrm{(yr)}$',fontsize=fsize)
            if xlims[1] == 0.:
                xlims = [xxx[0],xxx[-1]]
        elif ixaxis =='model_number':
            xxx= self.get('model_number')[modstart:modstop]
            print('plot versus model number')
            ax.set_xlabel('Model number')#,fontsize=fsize)
            if xlims[1] == 0.:
                xlims = [self.get('model_number')[modstart],self.get('model_number')[modstop]]
        elif ixaxis =='age':
            if t0_model != 0:
                t0_mod=np.abs(mod-t0_model).argmin()
                xxx= self.get('star_age')[modstart:modstop] - self.get('star_age')[t0_mod]
                print('plot versus age')
                ax.set_xlabel('Age [yr] - '+str(self.get('star_age')[modstart]))#,fontsize=fsize)
            else:
                xxx= old_div(self.get('star_age')[modstart:modstop],1.e6)
                ax.set_xlabel('Age [Myr]')#,fontsize=fsize)
            if xlims[1] == 0.:
                xlims = [xxx[0],xxx[-1]]

        ax.set_ylabel('$\mathrm{Mass }(M_\odot)$')

        # some stuff for rasterizing only the contour part of the plot, for nice, but light, eps:
        class ListCollection(Collection):
            def __init__(self, collections, **kwargs):
                Collection.__init__(self, **kwargs)
                self.set_collections(collections)
            def set_collections(self, collections):
                self._collections = collections
            def get_collections(self):
                return self._collections
            @allow_rasterization
            def draw(self, renderer):
                for _c in self._collections:
                    _c.draw(renderer)

        def insert_rasterized_contour_plot(c):
            collections = c.collections
            for _c in collections:
                _c.remove()
            cc = ListCollection(collections, rasterized=True)
            ax = pl.gca()
            ax.add_artist(cc)
            return cc

        cmapMIX = matplotlib.colors.ListedColormap(['w','#8B8386']) # rose grey
        cmapB1  = pyl.cm.get_cmap('Blues')
        cmapB2  = pl.cm.get_cmap('Reds')

        ylims1=[0.,0.]
        ylims1[0]=ylims[0]
        ylims1[1]=ylims[1]
        if ylims == [0.,0.]:
            ylims[0] = 0.
            ylims[1] = mup
        if ylims[0] != 0.:
            ylab='$(\mathrm{Mass }$ - '+str(ylims[0])
            if yscale!='1.':
                ylab+=') / '+yscale+' $M_\odot$'
            else:
                ylab+=') / $M_\odot$'
            ax.set_ylabel(ylab)
            y = y - ylims[0]
            y = y*float(yscale) # SJONES tweak
            ylims[0] = y[0]
            ylims[1] = y[-1]

        print('plotting contours')
        CMIX    = ax.contourf(xxx[::dx],y,Z, cmap=cmapMIX,alpha=0.6,levels=[0.5,1.5])
        if rasterise==True:
            insert_rasterized_contour_plot(CMIX)
        if outlines == True:
            CMIX_outlines    = ax.contour(xxx[::dx],y,Z, cmap=cmapMIX)
            if rasterise==True:
                insert_rasterized_contour_plot(CMIX_outlines)

        if engenstyle == 'full' and engenPlus == True:
            if engenlevels!= None:
                CBURN1  = ax.contourf(xxx[::dx],y,B1, cmap=cmapB1, alpha=0.5, locator=matplotlib.ticker.LogLocator(),levels=engenlevels)
                if outlines:
                    CB1_outlines  = ax.contour(xxx[::dx],y,B1, cmap=cmapB1, alpha=0.7, locator=matplotlib.ticker.LogLocator(),levels=engenlevels)
            else:
                CBURN1  = ax.contourf(xxx[::dx],y,B1, cmap=cmapB1, alpha=0.5, locator=matplotlib.ticker.LogLocator())
                if outlines:
                    CB1_outlines  = ax.contour(xxx[::dx],y,B1, cmap=cmapB1, alpha=0.7, locator=matplotlib.ticker.LogLocator())
            CBARBURN1 = pyl.colorbar(CBURN1)
            CBARBURN1.set_label('$|\epsilon_\mathrm{nuc}-\epsilon_{\\nu}| \; (\mathrm{erg\,g}^{-1}\mathrm{\,s}^{-1})$',fontsize=fsize)
            if rasterise==True:
                insert_rasterized_contour_plot(CBURN1)
                if outlines:
                    insert_rasterized_contour_plot(CB1_outlines)

        if engenstyle == 'full' and engenMinus == True:
            CBURN2  = ax.contourf(xxx[::dx],y,B2, cmap=cmapB2, alpha=0.5, locator=matplotlib.ticker.LogLocator())
            if outlines:
                CBURN2_outlines  = ax.contour(xxx[::dx],y,B2, cmap=cmapB2, alpha=0.7, locator=matplotlib.ticker.LogLocator())
            CBARBURN2 = pl.colorbar(CBURN2)
            if engenPlus == False:
                CBARBURN2.set_label('$|\epsilon_\mathrm{nuc}-\epsilon_{\\nu}| \; (\mathrm{erg\,g}^{-1}\mathrm{\,s}^{-1})$',fontsize=fsize)
            if rasterise==True:
                insert_rasterized_contour_plot(CBURN2)
                if outlines:
                    insert_rasterized_contour_plot(CB2_outlines)

        if engenstyle == 'twozone' and (engenPlus == True or engenMinus == True):
            ax.contourf(xxx[::dx],y,V, cmap=cmapB1, alpha=0.5)

        print('plotting patches')
        mtot=self.get('star_mass')[modstart:modstop][::dx]
        mtot1=(mtot-ylims1[0])*float(yscale)
        ax.plot(xxx[::dx],mtot1,'k-')

        if boundaries == True:
            print('plotting abund boundaries')
            try:
                bound=self.get('h1_boundary_mass')[modstart:modstop]
                bound1=(bound-ylims1[0])*float(yscale)
                ax.plot(xxx,bound1,label='H boundary',linestyle='-')

                bound=self.get('he4_boundary_mass')[modstart:modstop]
                bound1=(bound-ylims1[0])*float(yscale)
                ax.plot(xxx,bound1,label='He boundary',linestyle='--')

                bound=self.get('c12_boundary_mass')[modstart:modstop]
                bound1=(bound-ylims1[0])*float(yscale)
                ax.plot(xxx,bound1,label='C boundary',linestyle='-.')

            except:
                try:
                    bound=self.get('he_core_mass')[modstart:modstop]
                    bound1=(bound-ylims1[0])*float(yscale)
                    ax.plot(xxx,bound1,label='H boundary',linestyle='-')

                    bound=self.get('c_core_mass')[modstart:modstop]-ylims[0]
                    bound1=(bound-ylims1[0])*float(yscale)
                    ax.plot(xxx,bound1,label='He boundary',linestyle='--')

                    bound=self.get('o_core_mass')[modstart:modstop]-ylims[0]
                    bound1=(bound-ylims1[0])*float(yscale)
                    ax.plot(xxx,bound1,label='C boundary',linestyle='-.')

                    bound=self.get('si_core_mass')[modstart:modstop]-ylims[0]
                    bound1=(bound-ylims1[0])*float(yscale)
                    ax.plot(xxx,bound1,label='C boundary',linestyle='-.')

                    bound=self.get('fe_core_mass')[modstart:modstop]-ylims[0]
                    bound1=(bound-ylims1[0])*float(yscale)
                    ax.plot(xxx,bound1,label='C boundary',linestyle='-.')

                except:
#                    print 'problem to plot boundaries for this plot'
                    pass
        if CO_ratio == True:
            surface_c12       = self.get('surface_c12')
            surface_o16       = self.get('surface_o16')
            COratio=old_div((surface_c12*4.),(surface_o16*3.))
            ax2=pyl.twinx()
            ax2.plot(xxx,COratio[modstart:modstop]-ylims[0],'-k',label='C/O ratio')
            ax2.axis([xlims[0],xlims[1],0,max(COratio)*1.1])
            ax2.legend(loc=1, fontsize=0.5*fsize)
            ax2.set_ylabel('C/O ratio')

        if plot_radius == True:
            ax2=pyl.twinx()
            try:
               ax2.plot(xxx,np.log10(self.get('he4_boundary_radius')[modstart:modstop]),label='He boundary radius',color='k',linewidth=1.,linestyle='-.')
            except:
               try:
                   ax2.plot(xxx,np.log10(self.get('c_core_mass')[modstart:modstop]),label='He boundary radius',color='k',linewidth=1.,linestyle='-.')
               except:
                   pass
            ax2.plot(xxx,self.get('log_R')[modstart:modstop],label='radius',color='k',linewidth=1.,linestyle='-.')
            ax2.set_ylabel('log(radius)')
        if rad_lines == True:
            ax2=pyl.twinx()
            for i in range(len(masses)):
                ax2.plot(modno,np.log10(rads[i]),color='k')

        ax.axis([xlims[0],xlims[1],ylims[0],ylims[1]])

        axm = pyl.subplot2grid((3,3), (2, 0), rowspan=1, colspan=3)

        Mdot=self.get('log_abs_mdot')[modstart:modstop]-ylims[0]
        axm.plot(xxx,Mdot,label='Log Mdot',linestyle='-')
        minMdot=min(Mdot)
        if minMdot==-99.0:
           minMdot=-20.0
        else:
           minMdot=minMdot
        axm.axis([xlims[0],xlims[1],minMdot,max(Mdot)*0.9])
        axm.legend(loc=4, fontsize=0.5*fsize)
        axm.set_ylabel('Log$_{10}$ |$\dot{M}$|')
        #axm.set_xlabel('Model number')

        axm2=pyl.twinx()
        LogLum=self.get('log_L')[modstart:modstop]-ylims[0]
        axm2.plot(xxx,LogLum,label='Log L',linestyle='-.')
        axm2.axis([xlims[0],xlims[1],min(LogLum)*1.1,max(LogLum)*1.1])
        axm2.set_ylabel('Log$_{10}$ L/L$_{\odot}$')
        axm2.legend(loc=1, fontsize=0.5*fsize)

        if outfile[-3:]=='png':
            fig.savefig(outfile,dpi=300)
        elif outfile[-3:]=='eps':
            fig.savefig(outfile,format='eps')
        elif outfile[-3:]=='pdf':
            fig.savefig(outfile,format='pdf')
        if showfig == True:
            pyl.show()
# we may or may not need this below
#        fig.clear()


    def find_first_TP(self):

        """
        Find first TP of the TPAGB phase and returns the model
        number at its LHe maximum.

        Parameters
        ----------

        """

        star_mass         = self.get('star_mass')
        he_lumi           = self.get('log_LHe')
        h_lumi            = self.get('log_LH')
        mx2_bot           = self.get('mx2_bot')*star_mass
        try:
           h1_boundary_mass  = self.get('h1_boundary_mass')
           he4_boundary_mass = self.get('he4_boundary_mass')
        except:
           try:
                h1_boundary_mass  = self.get('he_core_mass')
                he4_boundary_mass = self.get('c_core_mass')
           except:
                pass


        TP_bot=np.array(self.get('conv_mx2_bot'))*np.array(self.get('star_mass'))
        TP_top=np.array(self.get('conv_mx2_top'))*np.array(self.get('star_mass'))
        lum_array=[]
        activate=False
        models=[]
        pdcz_size=[]
        for i in range(len(h1_boundary_mass)):
            if (h1_boundary_mass[i]-he4_boundary_mass[i] <0.2) and (he4_boundary_mass[i]>0.2):
                if (mx2_bot[i]>he4_boundary_mass[i]) and (he_lumi[i]>h_lumi[i]):
                        if TP_top[i]>he4_boundary_mass[i]:
                                pdcz_size.append(TP_top[i]-TP_bot[i])
                                activate=True
                                lum_array.append(he_lumi[i])
                                models.append(i)
                                #print(TP_bot[i],TP_top[i])
                if (activate == True) and (he_lumi[i]<h_lumi[i]):
                        #if fake tp
                        if max(pdcz_size)<1e-5:
                                active=False
                                lum_array=[]
                                models=[]
                                print('fake tp')
                        else:
                                break
        t0_model = models[np.argmax(lum_array)]
        return t0_model


    def find_TPs_and_DUPs(self, percent=5., makefig=False):
        """
        Function which finds TPs and uses the calc_DUP_parameter
        function.  To calculate DUP parameter evolution dependent of
        the star or core mass.

        Parameters
        ----------
        fig : integer
            Figure number to plot.
        t0_model : integer
            First he-shell lum peak.
        percent : float
            dredge-up is defined as when the mass dredged up is a certain
            percent of the total mass dredged up during that event, which
            is set by the user in this variable.
            The default is 5.
        makefig :
            do you want a figure to be made?

        Returns
        -------
        TPmods : array
            model numbers at the peak of each thermal pulse
        DUPmods : array
            model numbers at the dredge-up, where dredge-up is defined as
            when the mass dredged up is a certain percent of the total mass
            dredged up during that event, which is set by the user
        TPend : array
            model numbers at the end of the PDCZ for each TP
        lambda : array
            DUP efficiency for each pulse
        """

        t0_model=self.find_first_TP()
        t0_idx=(t0_model-self.get("model_number")[0])
        first_TP_he_lum=10**(self.get("log_LHe")[t0_idx])
        he_lum=10**(self.get("log_LHe")[t0_idx:])
        h_lum=10**(self.get("log_LH")[t0_idx:])
        model=self.get("model_number")[t0_idx:]

        try:
           h1_bndry=self.get("h1_boundary_mass")[t0_idx:]
        except:
           try:
                h1_bndry=self.get('he_core_mass')[t0_idx:]
           except:
                pass
        # SJ find TPs by finding local maxima in He-burning luminosity and
        # checking that the he_lum is greater than the h_lum:
        maxima=[0]
        for i in range(2,len(model)-1):
            if he_lum[i] > he_lum[i-1] and he_lum[i] > he_lum[i+1]:
                if he_lum[i-1] > he_lum[i-2] and he_lum[i+1] > he_lum[i+2]:
                    if he_lum[i] > h_lum[i]:
                        maxima.append(i)

        # find DUPs when h-boundary first decreases by more than XX% of the total DUP
        # depth:
        DUPs=[]
        TPend=[]
        maxDUPs=[]
        for i in range(len(maxima)):
            idx1=maxima[i]
            try:
                idx2=maxima[i+1]
            except IndexError:
                idx2=-1
            bound=h1_bndry[idx1:idx2]
            bound0=bound[0]
            if bound0==min(bound) or bound0 < min(bound): # then no DUP
                DUP=idx1
                DUPs.append(DUP)
                maxDUPs.append(DUP)
            else:
                maxDUPs.append(idx1+bound.argmin()) # model number of deepest extend of 3DUP
                maxDUP=bound0-min(bound) # total mass dredged up in DUP
                db=bound - bound[0]
                db_maxDUP = old_div(db, maxDUP)
                DUP=np.where(db_maxDUP <= old_div(-float(percent),100.))[0][0]
                DUPs.append(DUP+idx1)
#                # Alternative definition, where envelope reaches mass coordinate
#                # where top of PDCZ had resided during the TP:
#                top=self.get('mx2_top')[idx1]
#                DUP=np.abs(bound-top).argmin()
#                DUPs.append(DUP+idx1)

        # find end of PDCZ by seeking from TP peak and checking mx2_bot:
            mx2b=self.get('mx2_bot')[t0_idx:][idx1:idx2]
            for i in range(len(mx2b)):
                if mx2b[i]==0.:
                    endTP=i+idx1
                    TPend.append(endTP)
                    break

        # 3DUP efficiency:
        lambd=[0.]
        for i in range(1,len(maxima)):
            dmenv = h1_bndry[maxima[i]] - h1_bndry[maxDUPs[i-1]]
            dmdredge = h1_bndry[maxima[i]] - h1_bndry[maxDUPs[i]]
            lambd.append(old_div(dmdredge,dmenv))

        TPmods = maxima + t0_idx
        DUPmods = DUPs + t0_idx
        TPend = TPend + t0_idx

        return TPmods, DUPmods, TPend, lambd

    def TPAGB_properties(self):

        """
        Temporary, use for now same function in nugrid_set.py!
        Returns many TPAGB parameters which are
        TPstart,TPmods,TP_max_env,TPend,min_m_TP,max_m_TP,DUPmods,DUPm_min_h
        Same function in nugrid_set.py.

        Parameters
        ----------

        """

        peak_lum_model,h1_mass_min_DUP_model=self.find_TP_attributes( 3, t0_model=self.find_first_TP(), color='r', marker_type='o')

        print('first tp')
        print(self.find_first_TP())
        print('peak lum mmmodel')
        print(peak_lum_model)
        print(h1_mass_min_DUP_model)

        TPmods=peak_lum_model

        DUPmods=h1_mass_min_DUP_model
        DUPmods1=[]
        for k in range(len(DUPmods)):
                DUPmods1.append(int(float(DUPmods[k]))+100) #to exclude HBB? effects

        DUPmods=DUPmods1



        TPstart=[]
        #find beginning of TP, goes from TP peak backwards
        # find end of PDCZ by seeking from TP peak and checking mx2_bot:
        models=self.get('model_number')
        mx2b_array=self.get('conv_mx2_bot')
        mx2t_array=self.get('conv_mx2_top')
        massbot=mx2b_array#*self.header_attr['initial_mass']
        masstop=mx2t_array#*self.header_attr['initial_mass']
        massenv=np.array(self.get('conv_mx1_bot'))*np.array(self.get('star_mass'))   #*self.header_attr['initial_mass']

        #h1_bdy=self.get('h1_boundary_mass')

        for k in range(len(TPmods)):
                idx=list(models).index(TPmods[k])
                mx2b=mx2b_array[:idx]
                for i in range(len(mx2b)-1,0,-1):
                        if mx2b[i]==0.:
                            startTP=models[i]
                            TPstart.append(int(float(startTP)))
                            break
        #Find end of TP, goes from TP forwards:
        TPend=[]
        max_m_TP=[]
        min_m_TP=[]
        DUP_m=[]
        TP_max_env=[]
        DUPm_min_h=[]
        flagdecline=False
        for k in range(len(TPmods)):
            idx=list(models).index(TPmods[k])
            mx2b=mx2b_array[idx:]
            mx2t=mx2t_array[idx:]
            refsize=mx2t[0]-mx2b[0]
            for i in range(len(mx2b)):
                if i==0:
                    continue
                if ((mx2t[i]-mx2b[i])<(0.5*refsize)) and (flagdecline==False):
                    flagdecline=True
                    refmasscoord=mx2t[i]
                    print('flagdecline to true')
                    continue
                if flagdecline==True:
                    if (mx2t[i]-mx2b[i])<(0.1*refsize):
                        #for the massive and HDUP AGB's where PDCZ conv zone becomes the Hdup CONV ZONE
                        if refmasscoord<mx2t[i]:
                            endTP=models[idx+i-1]
                            TPend.append(int(float(endTP)))
                            print('HDUp, TP end',endTP)
                            break
                        if (mx2t[i]-mx2b[i])<1e-5:
                            endTP=models[idx+i-1]
                            TPend.append(int(float(endTP)))
                            print('normal TPend',endTP)
                            break

                # if max(mx2t[0:(i-1)])>mx2t[i]:
                #         (max(mx2t[0:(i-1)]) - min(mx2b[0:(i-1)]))
                #         flag=True
                #         continue
                # if flag==True:
                #     endidx=idx+i
                #     endTP=models[endidx]
                #     TPend.append(int(float(endTP)))

                # if (mx2t[i]-mx2b[i])<1e-5:                        #mx2b[i])==0.:
                #     endidx=idx+i
                #     endTP=models[endidx]
                #     TPend.append(int(float(endTP)))
                #     break

            print('found TP boundaries',TPstart[-1],TPend[-1])
        #find max and minimum mass coord of TP at max Lum
            mtot=self.get('star_mass')
            masstop_tot=np.array(masstop)*np.array(mtot)
            idx_tpext=list(masstop_tot).index(max(masstop_tot[TPstart[k]:(TPend[k]-10)]))
            print('TP',k+1,TPmods[k])
            print(TPstart[k],TPend[k])
            print('INDEX',idx_tpext,models[idx_tpext])
            print(max(masstop_tot[TPstart[k]:(TPend[k]-10)]))
            mtot=self.get('star_mass')[idx_tpext]
            max_m_TP.append(masstop[idx_tpext]*mtot)
            min_m_TP.append(massbot[idx_tpext]*mtot)

            TP_max_env.append(massenv[idx_tpext])#*mtot)
            if k> (len(DUPmods)-1):
                    continue
            idx=list(models).index(DUPmods[k])
            mtot=self.get('star_mass')[idx]
            #DUP_m.append(h1_bdy[idx])#*mtot)
        #######identify if it is really a TDUP, Def.
            try:
                h1_bndry=self.get("h1_boundary_mass")[t0_idx:]
            except:
                try:
                    h1_bndry=self.get('he_core_mass')[t0_idx:]
                except:
                    pass

            if h1_bndry[idx]>=max_m_TP[-1]:
                print('Pulse',k+1,'model',TPmods[k],'skip')
                print(h1_bndry[idx],max_m_TP[-1])
                DUPmods[k] = -1
                DUPm_min_h.append( -1)
                continue

            DUPm_min_h.append(h1_bdy[idx])
        for k in range(len(TPmods)):
            print('#############')
            print('TP ',k+1)
            print('Start: ',TPstart[k])
            print('Peak' , TPmods[k],TP_max_env[k])
            print('(conv) PDCZ size: ',min_m_TP[k],' till ',max_m_TP[k])
            print('End',TPend[k])
            if k <=(len(DUPmods)-1):
                print(len(DUPmods),k)
                print('DUP max',DUPmods[k])
                print(DUPm_min_h[k])
            else:
                print('no DUP')

            return TPstart,TPmods,TP_max_env,TPend,min_m_TP,max_m_TP,DUPmods,DUPm_min_h


    def find_TP_attributes(self, t0_model=0, fig=10, color='k', marker_type='*',
                           h_core_mass=False, no_fig=False):
        """

        Function which finds TPs and uses the calc_DUP_parameter
        function.  To calculate DUP parameter evolution dependent of
        the star or core mass.
        
        Parameters
        ----------
        fig : integer
        Figure number to plot.
        t0_model : integer
        First he-shell lum peak
        color : string
        Color of the plot.
        marker_type : string
        marker type.
        h_core_mass : boolean, optional
        If True: plot dependence from h free core , else star mass.
        The default is False.
        no_fig : boolean, optional
        The default is False.
        
        """
       
        #if len(t0_model)==0:
        
        t0_idx=(t0_model-self.get("model_number")[0])
	#for first TPi

        first_TP_he_lum=10**(self.get("log_LHe")[t0_idx])
        he_lum=10**(self.get("log_LHe")[t0_idx:])
        h_lum=10**(self.get("log_LH")[t0_idx:])
        model=self.get("model_number")[t0_idx:]
        try:
           h1_bndry=self.get("h1_boundary_mass")[t0_idx:]
           he4_bdy=self.get("he4_boundary_mass")[t0_idx:]
        except:
           try:
                h1_bndry=self.get('he_core_mass')[t0_idx:]	
                he4_bdy=self.get("c_core_mass")[t0_idx:]	
           except:
                pass
        TP_bot=np.array(self.get('conv_mx2_bot')[t0_idx:])*np.array(self.get('star_mass')[t0_idx:])
        TP_top=np.array(self.get('conv_mx2_top')[t0_idx:])*np.array(self.get('star_mass')[t0_idx:])

	

        #define label
        z=self.header_attr["initial_z"]
        mass=self.header_attr["initial_mass"]
        leg=str(mass)+"M$_{\odot}$ Z= "+str(z)
        peak_lum_model=[]
        peak_lum_save=[]
        h1_mass_tp=[]
        h1_mass_min_DUP_model=[]
        ##TP identification with he lum if within 1% of first TP luminosity
        perc=0.01
        min_TP_distance=300 #model
        lum_1=[]
        model_1=[]
        h1_mass_model=[]
        TP_counter=0
        new_TP=True
        TP_interpulse=False
        interpulse_counter=0

        TP_size=TP_top[0]-TP_bot[0]
        lastDUP=False
        for i in range(len(he_lum)):
                #interpulse_counter+=1
            #if (h_lum[i]<he_lum[i]):
                #interpulse_counter=0
                #new_TP=True
                #TP_interpulse=False
                #if i > 0:
                #    h1_mass_1.append(h1_bndry[i])
                #    h1_mass_model.append(model[i])
	    #in case when He-lum is dominating till the end of the calculation (He-burner?
            if (TP_interpulse ==False) and (h_lum[i]<he_lum[i]):
               if (len(he_lum)-1)==i:
                        TP_interpulse=True	
			#interpulse_counter=1000 #value higher than 200

	    #if simulation stops during a TP
            if (len(he_lum)-1)==i:
                if (h_lum[i]<he_lum[i]):
                        if  (he4_bdy[i]<TP_bot[i]):
                                lastDUP=True
                                break 

            #print i
            if i ==0:
                 lum_1.append(first_TP_he_lum)
                 model_1.append(t0_model)
		
            if True: #(TP_interpulse==False):
                 if (h_lum[i]<he_lum[i]):
                      if  (he4_bdy[i]<TP_bot[i]) and ( (TP_top[i]-TP_bot[i]) > TP_size*0.1) :
                                lum_1.append(he_lum[i])
                                model_1.append(model[i])
		 		#print 'peak at model',model[i]
                                TP_interpulse=True	  	
       
            if ( ((len(he_lum)-1)==i) and (h_lum[i]>he_lum[i])) and TP_interpulse==False:
                        #print 'test for last DUP'
                        if min(h1_bndry[peak_lum_model[-1]-t0_idx:i])<h1_bndry[peak_lum_model[-1]-t0_idx]:
                                #print 'last DUP after last TP'
                                lastDUP=True
       	                        break
				

            if ((h_lum[i]>he_lum[i]) and (TP_interpulse==True)) or ( ((len(he_lum)-1)==i) and (TP_interpulse==True)):
		#make sure that pulse is fully computed
               	if (len(he_lum)-1)<(i+2000):
			#print 'test for last DUP'
                        if min(h1_bndry[peak_lum_model[-1]-t0_idx:i])<h1_bndry[peak_lum_model[-1]-t0_idx]:
				#print 'last DUP after last TP'
                                lastDUP=True	
                        break
		#print 'model',model[i] 
		#print 'lum1',lum_1
                TP_interpulse=True
                max_value1=np.array(lum_1).max()
                if len(lum_1)<10:
                        continue
                if len(peak_lum_model)>2:
                        if ( 10**((np.log10(prev_he_lum_start)+np.log10(peak_lum_save[-1]))/2.)  )  >max_value1:
                             continue
                max_value=np.array(lum_1).max()
                max_index = lum_1.index(max_value)
                prev_he_lum_start=lum_1[0]

                peak_lum_save.append(max_value)
                #print max_index,i
                peak_lum_model.append(model_1[max_index])
                #for DUP calc
		#print 'peak model',(model_1[max_index])
		#print 'current peak lum',peak_lum_model[-1]
                #max_lum_idx=h1_mass_model.index(model_1[max_index])
                #min_value=np.array(h1_mass_1[max_lum_idx:]).min()
                #min_index = h1_mass_1.index(min_value)
		#if interpulse_counter<1000:
                #	h1_mass_min_DUP_model.append(h1_mass_model[min_index])
		#else:
		#	h1_mass_min_DUP_model.append(-1)
                #TP_counter+=1
                lum_1=[]
                model_1=[]
                TP_interpulse=False
                #h1_mass_1=[]i
                #h1_mass_model=[]
                #new_TP=False
        #TP_interpulse=False
	#here check if h1_mass_min_DUP_model is really at the lowest point
        for k in range(len(peak_lum_model)-1):
                idx1=list(model).index(peak_lum_model[k])
                idx2=list(model).index(peak_lum_model[k+1])
                h1_mass_min_DUP_model.append(model[list(h1_bndry).index(min(h1_bndry[idx1:idx2]))]	)
        if lastDUP==True:
                idx1=list(model).index(peak_lum_model[-1])	
                idx2=-1
                h1_mass_min_DUP_model.append(model[list(h1_bndry).index(min(h1_bndry[idx1:idx2]))]      )
	#print peak_lum_model
	#print h1_mass_min_DUP_model
        #print peak_lum_model
        #print h1_mass_min_DUP_model
        #print h1_mass_tp
        modeln=[]
        return peak_lum_model,h1_mass_min_DUP_model



    def calc_DUP_parameter(self, modeln, label, fig=10, color='r', marker_type='*',
                           h_core_mass=False):
        """
        Method to calculate the DUP parameter evolution for different
        TPs specified specified by their model number.

        Parameters
        ----------
        fig : integer
            Figure number to plot.
        modeln : list
            Array containing pairs of models each corresponding to a
            TP. First model where h boundary mass will be taken before
            DUP, second model where DUP reaches lowest mass.
        leg : string
            Plot label.
        color : string
            Color of the plot.
        marker_type : string
            marker type.
        h_core_mass : boolean, optional
            If True: plot dependence from h free core , else star mass.
            The default is False.

        """
        number_DUP=(old_div(len(modeln),2) -1) #START WITH SECOND
        try:
           h1_bnd_m=self.get('h1_boundary_mass')
        except:
           try:
                h1_bnd_m=self.get('he_core_mass')
           except:
                pass
        star_mass=self.get('star_mass')
        age=self.get("star_age")
        firstTP=h1_bnd_m[modeln[0]]
        first_m_dredge=h1_bnd_m[modeln[1]]
        DUP_parameter=np.zeros(number_DUP)
        DUP_xaxis=np.zeros(number_DUP)
        j=0
        for i in np.arange(2,len(modeln),2):
            TP=h1_bnd_m[modeln[i]]
            m_dredge=h1_bnd_m[modeln[i+1]]
            if i ==2:
                last_m_dredge=first_m_dredge
            #print "testest"
            #print modeln[i]
            if h_core_mass==True:
                DUP_xaxis[j]=h1_bnd_m[modeln[i]]                        #age[modeln[i]] - age[modeln[0]]
            else:
                DUP_xaxis[j]=star_mass[modeln[i]]
            #DUP_xaxis[j]=modeln[i]
            DUP_parameter[j]=old_div((TP-m_dredge),(TP-last_m_dredge))
            last_m_dredge=m_dredge
            j+=1

        pl.figure(fig)
        pl.rcParams.update({'font.size': 18})
        pl.rc('xtick', labelsize=18)
        pl.rc('ytick', labelsize=18)

        pl.plot(DUP_xaxis,DUP_parameter,marker=marker_type,markersize=12,mfc=color,color='k',linestyle='-',label=label)
        if h_core_mass==True:
            pl.xlabel("$M_H$",fontsize=20)
        else:
            pl.xlabel("M/M$_{\odot}$",fontsize=24)
        pl.ylabel("$\lambda_{DUP}$",fontsize=24)
        pl.minorticks_on()
        pl.legend()



class star_log(history_data):
    """
    Class derived from history_data class (copy). Existing just (for
    compatibility reasons) for older mesa python scripts.

    """



# below are some utilities that the user typically never calls directly


def _read_mesafile(filename,data_rows=0,only='all'):
    """ private routine that is not directly called by the user"""
    f=open(filename,'r')
    vv=[]
    v=[]
    lines = []
    line  = ''
    for i in range(0,6):
        line = f.readline()
        lines.extend([line])

    hval  = lines[2].split()
    hlist = lines[1].split()
    header_attr = {}
    for a,b in zip(hlist,hval):
        header_attr[a] = float(b)
    if only is 'header_attr':
        return header_attr

    cols    = {}
    colnum  = lines[4].split()
    colname = lines[5].split()
    for a,b in zip(colname,colnum):
        cols[a] = int(b)

    data = []

    old_percent = 0
    for i in range(data_rows):
        # writing reading status
        percent = int(i*100/np.max([1, data_rows-1]))
        if percent >= old_percent + 5:
            sys.stdout.flush()
            sys.stdout.write("\r reading " + "...%d%%" % percent)
            old_percent = percent
        line = f.readline()
        v=line.split()
        try:
            vv=np.array(v,dtype='float64')
        except ValueError:
            for item in v:
                if item.__contains__('.') and not item.__contains__('E'):
                    v[v.index(item)]='0'
        data.append(vv)

    print(' \n')
    f.close()
    a=np.array(data)
    data = []
    return header_attr, cols, a


def _cleanstarlog(file_in):
    """
    cleaning history.data or star.log file, e.g. to take care of
    repetitive restarts.

    private, should not be called by user directly

    Parameters
    ----------
    file_in : string
        Typically the filename of the mesa output history.data or
        star.log file, creates a clean file called history.datasa or
        star.logsa.

    (thanks to Raphael for providing this tool)

    """

    file_out=file_in+'sa'
    f = open(file_in)
    lignes = f.readlines()
    f.close()

    nb    = np.array([],dtype=int)   # model number
    nb    = np.concatenate((nb    ,[  int(lignes[len(lignes)-1].split()[ 0])]))
    nbremove = np.array([],dtype=int)   # model number
    i=-1

    for i in np.arange(len(lignes)-1,0,-1):
        line = lignes[i-1]

        if i > 6 and line != "" :
            if int(line.split()[ 0])>=nb[-1]:
                nbremove = np.concatenate((nbremove,[i-1]))
            else:
                nb = np.concatenate((nb    ,[  int(line.split()[ 0])]))
    i=-1
    for j in nbremove:
        lignes.remove(lignes[j])

    fout = open(file_out,'w')
    for j in np.arange(len(lignes)):
        fout.write(lignes[j])
    fout.close()
