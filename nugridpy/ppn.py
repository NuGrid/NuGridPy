
#
# NuGridpy - Tools for accessing and visualising NuGrid data.
#
# Copyright 2007 - 2014 by the NuGrid Team.
# All rights reserved. See LICENSE.
#

"""
Analyze and visualize the output of single-zone ppn simulations

Contains:

class xtime: for analyzing x-time.data time evolutions

class abu_vector: for analyzing iso_massfxxxx.DAT time evolutions

PPN Assumptions
===============
The first non white space character in a header line is a '#'.  Header
attributes are separated from their value by white space or by white
space surrounding an equals sign.

A Header attribute is separated by the previous Header attribute by
white space or a line break.

There are only 6 data columns.  The first being the number, second
being Z, third being A, fourth being isomere state, fifth being
abundance_yps, and sixth the element name.

The first five columns consist purely of numbers, no strings are
allowed. Of the values in the final column, the name column, the first
two are letters specifying the element name, and the rest are spaces or
numbers (in that strict order), except for the element names: Neut and
Prot.

All the profile files in the directory have the same cycle attributes.
The cycle numbers of the 'filename'+xxxxx start at 0.

PPN files always end in .DAT and are not allowed any '.'.

There can not be any blank lines in the data files.

No cycle numbers are skipped, ie if cycle 0 and 3 are present, 1 and 2
must be here as well.

"""
from __future__ import print_function
from __future__ import absolute_import

from builtins import zip
from builtins import str
from builtins import range

from numpy import *
import matplotlib
from matplotlib.pylab import *
import os

from .data_plot import *
from .utils import *
from . import utils

#import pdb

class xtime(DataPlot):
    '''
    read and plot x-time.dat output files.

    Parameters
    ----------
    sldir : string, optional
        The directory of the specified filename.  The default is './'.
    fname : string, optional
        Specify alternative filename of file of type x-time.dat.  The
        default is 'x-time.dat'.

    Returns
    -------
    xtime
        A xtime instance.

    Examples
    --------
    >>> import ppn
    >>> f=ppn.xtime()
    There are 1099 species found.
    There are 19 time steps found.

    >>> f.col
    f.col_num  f.col_tot  f.cols

    >>> f.cols[0:10]
    ['age',
    't9',
    'rho',
    'sum_yps',
    'NEUT',
    'PROT',
    'H   2',
    'HE  3',
    'HE  4',
    'BE  7']

    >>> f.plot('HE  4')
    age HE  4

    >>> f.plot('C  12')
    age C  12

    >>> f.plot('PROT')

    '''
    sldir  = '' #Standard Directory
    data = []
    cols=[]     #list of column attribute names
    col_num = []#Dict of column attribute names and their associated values
    xdat = []
    ydat = []

    def __init__(self, sldir='./', fname='x-time.dat'):
        '''
        A xtime instance

        Parameters
        ----------
        sldir : string, optional
            The directory of the specified filename.  The default is './'.
        fname : string, optional
            Specify alternative filename of file of type x-time.dat.  The
            default is 'x-time.dat'.

        Returns
        -------
        xtime
            A xtime instance.

        '''
        self.sldir= sldir

        if not os.path.exists(sldir):  # If the path does not exist
            print('error: Directory, '+sldir+ ' not found')
            print('Now returning None')
            return None
        else:
            self.data, self.col_num, self.cols, self.col_tot, self.ilines =\
            self._readFile(fname,sldir)


    def _readFile(self, fname, sldir):
        '''
        Private method that reads in the data file and organizes it
        within this object.

        '''
        if sldir.endswith('/'):
            fname = str(sldir)+str(fname)
        else:
            fname = str(sldir)+'/'+str(fname)
        f=open(fname,'r')

        # read header line
        line=f.readline()
        cols  = []
        ispec = 0
        for i in range(1,len(line.split('|'))):
            col = line.split('|')[i].strip()
            if '-' in col:
                ispec += 1
                col   = col.split('-')[1]
            cols.append(col)
            col_num={}
        col_tot = len(cols)

        print('number of species: ', str(ispec))
        print('number of cols: ', str(col_tot))


        col_num={}
        for a,b in zip(cols,list(range(col_tot))):
            col_num[a]=b

        # read remainder of the file
        lines=f.readlines()
        data=[]
        for i in range(len(lines)):
            v=lines[i].split()
            vv=array(v,dtype='float')
            data.append(vv)
        ilines=i
        print("There are "+str(ilines)+" time steps found.")
        return data,col_num,cols,col_tot,ilines

    def get(self, col_str):
        '''
        get one data column with the data

        Parameters
        ----------
        col_str : string
            One of the column strings in self.cols.

        '''
        data_column=zeros(self.ilines)
        for i in range(self.ilines):
            data_column[i]=self.data[i][self.col_num[col_str]]
        return data_column

    def plot_xtime(self, y, x='time', label='default', labelx=None,
                   labely=None ,title=None, shape='.', logx=False,
                   logy=True, base=10):
        '''
        make a simple plot of two columns against each other.

        An example would be instance.plot_xtime('PB206', label='PB206 vs t_y'
        Recomend using the plot function DataPlot.plot() it has more
        functionality.

        Parameters
        ----------
        Y : string
            Column on Y-axis.
        X : string, optional
            Column on X-axis.  The default is "time".
        label : string, optional
            Legend label.  The default is "default".
        labelX : string, optional
            The label on the X axis.  The default is None.
        labelY : string, optional
            The label on the Y axis.  The default is None.
        title : string, optional
            The Title of the Graph.  The default is None.
        shape : string, optional
            What shape and colour the user would like their plot in.
            The default is '.'.
        logX : boolean, optional
            A boolean of weather the user wants the x axis
            logarithmically.  The default is False.
        logY : boolean, optional
            A boolean of weather the user wants the Y axis
            logarithmically.  The default is True.
        base : integer, optional
            The base of the logarithm.  The default is 10.

        Notes
        -----
        For all possable choices visit,
        <http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot>

        '''
        if label is 'default':
            lab_str=y
        else:
            lab_str=label

        try:
            self.get(x)
        except KeyError:
            x='age'

        DataPlot.plot(self,x,y,legend=lab_str,labelx=labelx, labely=labely,
                      title=title, shape=shape,logx=logx, logy=logy, base=base)
        '''
        print X,Y
        xdat=self.get(X)
        ydat=self.get(Y)
        self.xdat = xdat
        self.ydat = ydat


        plot(xdat,log10(ydat),label=lab_str)
        legend()
        '''

class abu_vector(DataPlot, Utils):
    '''
    Class for reading selem00xxxx.DAT files

    Parameters
    ----------
    sldir : string, optional
        where fname exists
    filenames : string, optional
        The default file names of the abundance vectors.  The
        defaults is "iso_massf".
    USEEPP : string, optional
        The default is "auto".

    Returns
    -------
    abu_vector
        A PPN instance.

    Examples
    --------
    >>> import ppn
    >>> p=ppn.abu_vector('./run')
    39 cycle numbers found in ./run
    Ranging from 0 to 38

    To find the cycle attributes:
    >>> p.cattrs
    ['mod', 'dzeit', 'agej', 't9', 'rho', 'densn', 'densp', 'densa']

    To find the data column attributes

    >>> p.dcols
    ['NUM', 'Z', 'A', 'ISOM', 'ABUNDANCE_MF', 'ISOTP']
    >>> p.get('Z',0)
    array([1, 2, 2, 4, 5, 3, 6, 6, 7, 7, 6, 7, 8, 8, 8, 9, 9, 9])
    >>> p.get('ABUNDANCE_MF',0)
    array([  1.43722000e-10,   1.00000000e-99,   9.81499000e-01,
    4.08738000e-20,   1.00000000e-99,   2.06944000e-21,
    3.42800000e-04,   9.62307000e-05,   1.05081000e-12,
    1.25518000e-02,   1.90131000e-08,   1.42230000e-07,
    4.98449000e-05,   4.80246000e-07,   4.24345000e-12,
    9.85201000e-17,   6.30866000e-16,   9.12726000e-11])

    or if the user wants the data from the first 3 cycles:

    >>> p.get('ABUNDANCE_MF',[0,1,2])
    [array([  1.43722000e-10,   ...,   9.81499000e-01,],
    array([  1.43722000e-10,   ...,   9.81499000e-01,],
    array([  1.43722000e-10,   ...,   9.81499000e-01,]
    >>> p.getElement('C 14',0)
    array([  1.50000000e+01,   6.00000000e+00,   1.40000000e+01,
    1.90131000e-08])
    >>> p.plot('abundance_yps', 'Z',0)

    plots data

    >>> p.iso_abund(0)

    Plots an isotope abundance distribution

    >>> p.abu_chart(0)

    Plots an isotope abundance chart

    One note about the plot functions, if instead of a single cycle the
    user inputs a list of cycles, the method will then, instead of
    plotting them, will then save a .png for each cycle.  Also if you
    just want a singular plot saved, the user can input their cycle,
    in a list like [0]. And that will save their plot.

    '''
    sldir = ''  #Standard Directory
    inputdir = '' # A copy of Standard Directory which never changes
    cattrs={} # cycle attributes
    dcols=[]  # list of the column attributes
    index=0   # index of were column data begins in the file
    files=[]  # list of files
    isotopes=[]# list of isotopes
    def __init__(self, sldir='./', filenames='iso_massf',
                 USEEPP='auto',verbose=False):
        '''
        initial method of this class

        Parameters
        ----------
        sldir : string, optional
            where fname exists
        filenames : string, optional
            The default file names of the abundance vectors.  The
            defaults is "iso_massf".
        USEEPP : string, optional
            The default is "auto".

        Returns
        -------
        abu_vector
            A PPN instance.

        '''
        self.debug=False
        self._stable_names() # provides in addition to stable_el from
                             # utils also just the stable element names
        self.sldir = sldir
        self.inputdir = ''
        self.startdir = os.getcwd()
        self.cattrs=[]
        self.dcols=[]
        self.files=[]
        self.isotopes=[]
        if not os.path.exists(sldir):  # If the path does not exist
            print('error: Directory, '+sldir+ ' not found')
            print('Now returning None')
            return None
        f=os.listdir(sldir) # reads the directory

        def iniabmode(): #Where we go if dealing with iniab files in a USEEPP directory
            print("---> In iniab mode")
            for file in f:
            # Removes any files that are not ppn files
                if filenames in file and 'ppn' in file and '~' not in file \
                and '#' not in file and  file[-1] is 'n' \
                and file[-2] is 'p' and file[-3] is 'p'\
                and 'restart' not in file:
                    self.files.append(file)
            self.files.sort()

            if len(self.files)==0:
                # If there are no Files in the Directory
                print('Error: no '+filenames+ ' named files exist in Directory')
                print('Now returning None')
                return None
            fname=self.files[len(self.files)-1]
            self.dcols,self.index=self._readPPN(fname,sldir)
            indexp_cyc2filels={}  # created index pointer from mod (cycle
            i = 0                 # name) to index in files array
            for file in self.files:
                mod=utils.iniabu(file)
                indexp_cyc2filels[mod] = i
                i += 1
            self.indexp_cyc2filels = indexp_cyc2filels

            for i in range(len(self.files)):
                self.files[i]=self.sldir+self.files[i]
            if not verbose:
                print(str(len(self.files))+' cycle numbers found in '+sldir)
                print('Ranging from '+str(self.files[0])[-9:-4]+' to '+str(self.files[len(self.files)-1])[-9:-4])
                print('Range may not be continuous. To display all available cycles, print <abu_vector_instance>.files')
            self.isotopes=self.get('ISOTP',self.files[0],numtype='file')

        def iso_massfmode():
            for file in f:
            # Removes any files that are not ppn files
                filelength=len(filenames)+4
                if filenames in file and 'DAT' in file and '~' not in file \
                and '#' not in file and len(file)>filelength \
                and 'restart' not in file:
                    self.files.append(file)
            self.files.sort()

            if len(self.files)==0:
                # If there are no Files in the Directory
                print('Error: no '+filenames+ ' named files exist in Directory')
                print('Now returning None')
                return None
            fname=self.files[len(self.files)-1]
            self.cattrs,self.dcols,self.index=self._readFile(fname,sldir)
            indexp_cyc2filels={}  # created index pointer from mod (cycle
            i = 0                 # name) to index in files array
            for file in self.files:
                mod=self.get('mod',fname=file,numtype='file')
                indexp_cyc2filels[mod] = i
                i += 1
            self.indexp_cyc2filels = indexp_cyc2filels

            for i in range(len(self.files)):
                self.files[i]=self.sldir+self.files[i]
            if not verbose:
                print(str(len(self.files))+' cycle numbers found in '+sldir)
                print('Ranging from '+str(self.files[0])[-9:-4]+' to '+str(self.files[len(self.files)-1])[-9:-4])
                print('Range may not be continuous. To display all available cycles, print <abu_vector_instance>.files')
            self.isotopes=self.get('ISOTP',self.files[0],numtype='file')


        '''
        In this chunk of code we find whether the specified directory includes USEEPP
        and if it does, go to 'iniabmode' otherwise (or if USEEPP='off' is thrown) we
        deal with 'iso_massfmode'
        '''
        if ('USEEPP' in sldir or 'USEEPP' in self.startdir) and (USEEPP == 'on' or USEEPP == 'auto'):
            print('--> Going to iniab mode')
            filenames='iniab'
            iniabmode()
        elif ('USEEPP' in sldir or 'USEEPP' in self.startdir) and USEEPP == 'off':
            print("--> Skipping iniab mode (USEEPP == 'off')")
            iso_massfmode()
        elif ('USEEPP' not in sldir or 'USEEPP' not in self.startdir) and USEEPP == 'on':
            print("--> Error: Not in, or referencing, a USEEPP directory. Cannot enter iniab mode")
            print("--> Skipping iniab mode")
            iso_massfmode()
        else:
            iso_massfmode()

    def getCycleData(self, attri, fname, numtype='cycNum'):
        """
        In this method a column of data for the associated cycle
        attribute is returned.

        Parameters
        ----------
        attri : string
            The name of the attribute we are looking for.
        fname : string
            The name of the file we are getting the data from or the
            cycle number found in the filename.
        numtype : string, optional
            Determines whether fname is the name of a file or, the
            cycle number.  If it is 'file' it will then  interpret it as
            a file, if it is 'cycNum' it will then  interpret it as a
            cycle number.  The default is "cycNum".

        """

        fname=self.findFile(fname,numtype)

        if self.inputdir == '':
            self.inputdir = self.sldir      # This chunk of code changes into the directory where fname is,
        os.chdir(self.inputdir)                           # and appends a '/' to the directory title so it accesses the
        self.sldir=os.getcwd() + '/'              # file correctly

        f=open(fname,'r')
        lines=f.readlines()

        if self.inputdir != './':                               #This chunk of code changes back into the directory you started in.
            os.chdir(self.startdir)
            self.sldir = self.inputdir


        for i in range(len(lines)):
            lines[i]=lines[i].strip()

        for i in range(len(lines)):
            if lines[i].startswith('#'):
                lines[i]=lines[i].strip('#')
                tmp=lines[i].split()
                tmp1=[]
                for j in range(len(tmp)):
                    if tmp[j] != '=' or '':
                        tmp1.append(tmp[j])
                tmp=tmp1
                for j in range(len(tmp)):
                    if tmp[j]== attri:
                        try:
                            if '.' in tmp[j+1]:
                                return float(tmp[j+1])
                            else:
                                return int(tmp[j+1])
                        except ValueError:
                            return str(tmp[j+1])

            elif lines[i].startswith('H'):
                continue
            else:
                print('This cycle attribute does not exist')
                print('Returning None')
                return None




    def getColData(self, attri, fname, numtype='cycNum'):
        """
        In this method a column of data for the associated column
        attribute is returned.

        Parameters
        ----------
        attri : string
            The name of the attribute we are looking for.
        fname : string
            The name of the file we are getting the data from or the
            cycle number found in the filename.
        numtype : string, optional
            Determines whether fname is the name of a file or, the
            cycle number.  If it is 'file' it will then  interpret it as
            a file, if it is 'cycNum' it will then  interpret it as a
            cycle number.  The default is "cycNum".

        """
        fname=self.findFile(fname,numtype)
        f=open(fname,'r')
        for i in range(self.index+1):
            f.readline()
        lines=f.readlines()
        for i in range(len(lines)):
            lines[i]=lines[i].strip()
            lines[i]=lines[i].split()
        index=0
        data=[]

        while index < len (self.dcols):
            if attri== self.dcols[index]:
                break
            index+=1

        for i in range(len(lines)):

            if index==5 and len(lines[i])==7:
                data.append(str(lines[i][index].capitalize())+'-'\
                +str(lines[i][index+1]))
            elif index==5 and len(lines[i])!=7:
                tmp=str(lines[i][index])
                if tmp[len(tmp)-1].isdigit():
                    tmp1=tmp[0]+tmp[1]
                    tmp1=tmp1.capitalize()
                    tmp2=''
                    for j in range(len(tmp)):
                        if j == 0 or j == 1:
                            continue
                        tmp2+=tmp[j]
                    data.append(tmp1+'-'+tmp2)
                elif tmp=='PROT':
                    data.append('H-1')
                elif tmp==('NEUT'or'NEUTR'or'nn'or'N   1'or'N-1'):
                    data.append('N-1')
                else:
                    data.append(tmp)
            elif index==0:
                data.append(int(lines[i][index]))
            else:
                data.append(float(lines[i][index]))

        return array(data)

    def getElement(self, attri, fname, numtype='cycNum'):
        '''
        In this method instead of getting a particular column of data,
        the program gets a particular row of data for a particular
        element name.

        attri : string
            The name of the attribute we are looking for. A complete
            list of them can be obtained by calling

            >>> get('element_name')
        fname : string
            The name of the file we are getting the data from or the
            cycle number found in the filename.
        numtype : string, optional
            Determines whether fname is the name of a file or, the
            cycle number.  If it is 'file' it will then  interpret it as
            a file, if it is 'cycNum' it will then  interpret it as a
            cycle number.  The default is "cycNum".

        Returns
        -------
        array
            A numpy array of the four element attributes, number, Z, A
            and abundance, in that order.

        Notes
        -----
        Warning
        '''
        element=[] #Variable for holding the list of element names
        number=[]  #Variable for holding the array of numbers
        z=[]       #Variable for holding the array of z
        a=[]       #Variable for holding the array of a
        abd=[]     #Variable for holding the array of Abundance
        data=[]    #variable for the final list of data

        fname=self.findFile(fname,numtype)
        f=open(fname,'r')
        for i in range(self.index+1):
            f.readline()
        lines=f.readlines()
        for i in range(len(lines)):
            lines[i]=lines[i].strip()
            lines[i]=lines[i].split()
        index=0
        data=[]

        while index < len (self.dcols):
            if attri== self.dcols[index]:
                break
            index+=1

        element=self.get(self.dcols[5],fname,numtype)
        number=[]
        z=[]
        a=[]
        isom=[]
        abd=[]
        for i in range(len(lines)):
            number.append(int(lines[i][0]))
            z.append(float(lines[i][1]))
            isom.append(float(lines[i][2]))
            abd.append(float(lines[i][1]))
        index=0 #Variable for determing the index in the data columns


        while index < len(element):
            if attri == element[index]:
                break

            index+=1

        data.append(number[index])
        data.append(z[index])
        data.append(a[index])
        data.append(isom[index])
        data.append(abd[index])

        return array(data)



    def get(self, attri, fname=None, numtype='cycNum', decayed=False):
        '''
        In this method all data for an entire cycle (basically the
        content of an iso_massfnnnn.DAT file) or a column of data for
        the associated attribute is returned.

        Parameters
        ----------
        attri : string or integer
            If attri is a string, attri is the cycle or name of the
            attribute we are looking for.

            If attri is an integer, attri is the cycle number (cycle arrays
            are not supported).
        fname : string, optional
            If attri is a string, fname is the name of the file we are
            getting the data from or the cycle number found in the
            filename, or a List of either cycles or filenames.  If fname
            is None, the data from all cycles is returned.

            If attri is an integer, then fname is not supported.

            The default is None.
        numtype : string, optional
            If attri is a string, numtype determines whether fname is
            the name of a file or, the cycle number.  If numtype is
            'file' it will then  interpret fname as a file.  If numtype
            is 'cycNum' it will then interpret fname as a cycle number.

            If attri is an Integer, then numtype is not supported.

            The default is "cycNum".
        decayed : boolean, optional
            If attri is a string, then decayed is not supported.

            If attri is an integer, then get instantaneously decay
            abundance distribution.

            The default is False.

        Returns
        -------
        array
            If attri is a string, data in the form of a numpy array is
            returned.

            If attri is an integer, Nothing is returned.

        Notes
        -----
        If attri is an integer, then the following variables will be
        added to the instance.

        a_iso_to_plot: mass number of plotted range of species.

        isotope_to_plot: corresponding list of isotopes.

        z_iso_to_plot: corresponding charge numbers.

        el_iso_to_plot: corresponding element names.

        abunds: corresponding abundances.

        isom: list of isomers with their abundances.

        '''
        if type(attri) is type(1):
            print("Calling get method in cycle mode, adding a_iso_to_plot, z.. el.. isotope.. isotope... to instance")
            self._getcycle(attri,decayed)
        elif type(attri) is type("string"):
            data=self._getattr(attri,fname,numtype)
            return data

    def _getcycle(self, cycle, decayed=False):
        ''' Private method for getting a cycle, called from get.'''
        yps=self.get('ABUNDANCE_MF', cycle)
        z=self.get('Z', cycle) #charge
        a=self.get('A', cycle) #mass
        isomers=self.get('ISOM', cycle)

        a_iso_to_plot,z_iso_to_plot,abunds,isotope_to_plot,el_iso_to_plot,isom=\
            self._process_abundance_vector(a,z,isomers,yps)
        self.a_iso_to_plot=a_iso_to_plot
        self.isotope_to_plot=isotope_to_plot
        self.z_iso_to_plot=z_iso_to_plot
        self.el_iso_to_plot=el_iso_to_plot
        self.abunds=array(abunds)
        self.isom=isom

        if decayed:
            try:
                self.decay_idp
            except AttributeError:
                print("WARNING: decayed in _getcycle ignores isomers " \
                    "and will decay alpha-unstable p-rich nuclei as if they were beta+ stable.")
                print("Initialising decay index pointers ....")
                self.decay_indexpointer() # provides self.decay_idp and
            ind_tmp=self.idp_to_stables_in_isostoplot

            isotope_decay=array(isotope_to_plot)[ind_tmp]
            z_iso_decay=array(z_iso_to_plot)[ind_tmp]
            a_iso_decay=array(a_iso_to_plot)[ind_tmp]
            el_iso_decay=array(el_iso_to_plot)[ind_tmp]
            abunds_decay=zeros(len(ind_tmp), dtype='float64')
            for i in range(len(isotope_to_plot)):
                idp=where(isotope_decay==isotope_to_plot[self.decay_idp[i]])[0] # points from
                # i on isotope_to_plot scale to decay target_on_decayed array scale
                abunds_decay[idp] += abunds[i]

            if self.debug:
                print("Decayed array:")
                for i in range(len(ind_tmp)):
                    print(isotope_decay[i], z_iso_decay[i], a_iso_decay[i], el_iso_decay[i], abunds_decay[i])

            self.a_iso_to_plot=a_iso_decay
            self.isotope_to_plot=isotope_decay
            self.z_iso_to_plot=z_iso_decay
            self.el_iso_to_plot=el_iso_decay
            self.abunds=abunds_decay


    def _getattr(self, attri, fname=None, numtype='cycNum'):
        ''' Private method for getting an attribute, called from get.'''
        if str(fname.__class__)=="<type 'list'>":
            isList=True
        else:
            isList=False

        data=[]
        if fname==None:
            fname=self.files
            numtype='file'
            isList=True
        if isList:
            for i in range(len(fname)):
                if attri in self.cattrs:
                    data.append(self.getCycleData(attri,fname[i],numtype))
                elif attri in self.dcols:
                    data.append(self.getColData(attri,fname[i],numtype))
                elif attri in self.get('ISOTP',fname,numtype):
                    data.append(self.getElement(attri,fname[i],numtype))
                else:
                    print('Attribute '+attri+ ' does not exist')
                    print('Returning none')
                    return None

        else:
            if attri in self.cattrs:
                return self.getCycleData(attri,fname,numtype)
            elif attri in self.dcols:
                return self.getColData(attri,fname,numtype)
            elif attri in self.get('ISOTP',fname,numtype):
                return self.getElement(attri,fname,numtype)
            else:
                print('Attribute '+attri+ ' does not exist')
                print('Returning none')
                return None

        return data


    def _readPPN(self, fname, sldir):
        '''
        Private method that reads in and organizes the .ppn file
        Loads the data of the .ppn file into the variable cols.

        '''
        if sldir.endswith(os.sep):
                    #Making sure fname will be formatted correctly
            fname = str(sldir)+str(fname)
        else:
            fname = str(sldir)+os.sep+str(fname)
            self.sldir+=os.sep
        f=open(fname,'r')
        lines=f.readlines()
        for i in range(len(lines)):
            lines[i]=lines[i].strip()

        cols = ['ISOTP', 'ABUNDANCE_MF'] #These are constant, .ppn files have no header to read from
        for i in range(len(lines)):
            if not lines[i].startswith('H'):
                index = i-1
                break

        return cols, index


    def _readFile(self, fname, sldir):
        '''
        private method that reads in and organizes the .DAT file
        Loads the data of the .DAT File into the variables cattrs and cols.
        In both these cases they are dictionaries, but in the case of cols,
        it is a dictionary of numpy array exect for the element ,
        element_name where it is just a list

        '''
        cattrs=[]
        if sldir.endswith(os.sep):
            #Making sure fname will be formatted correctly
            fname = str(sldir)+str(fname)
        else:
            fname = str(sldir)+os.sep+str(fname)
            self.sldir+=os.sep
        f=open(fname,'r')
        lines=f.readlines()
        for i in range(len(lines)):
            lines[i]=lines[i].strip()


        cols=lines[0].strip('H')
        cols=cols.strip()
        cols=cols.split()
        for i in range(len(lines)):
            if lines[i].startswith('#'):
                # if it is a cycle attribute line
                lines[i]=lines[i].strip('#')
                tmp=lines[i].split()
                tmp1=[]
                for j in range(len(tmp)):
                    if tmp[j] != '=' or '':
                        tmp1.append(tmp[j])
                tmp=tmp1

                j=0
                while j <len(tmp):
                    cattrs.append(tmp[j])
                    j+=2


            elif not lines[i].startswith('H'):
                index = i-1
                break


        return cattrs,cols, index

    def findFile(self, fname, numtype):
        """
        Function that finds the associated file for fname when Fname is
        time or NDump.

        Parameters
        ----------
        fname : string
            The name of the file we are looking for.
        numType : string
            Designates how this function acts and how it interprets
            fname.  If numType is 'file', this function will get the
            desired attribute from that file.  If numType is 'cycNum',
            this function will get the desired attribute from that file
            with fname's model number.

        """
        numType=numtype.upper()
        if numType == 'FILE':
                #do nothing
            return fname
        elif numType == 'CYCNUM':
            try:
                fname = int(fname)
            except ValueError:
                print('Improper choice:'+ str(fname))
                print('Reselecting as 0')
                fname = 0
                print('Using '+self.files[fname])
        try:
            return self.files[self.indexp_cyc2filels[fname]]
        except IndexError:
            mods = array(self.get('mod'), dtype=int)
            if fname not in mods:
                print('You seem to try to plot a cycle that is not present: '+str(fname))
                fname = mods[-1]
                print('I will assume you want to plot the last cycle in the run: '+str(fname))
                print('[I am not 100% sure this escape is debugged. You better do this again with')
                print('the correct input.]')
                return self.files[fname]
