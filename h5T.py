#!/usr/bin/python

#
# NuGridpy - Tools for accessing and visualising NuGrid data.
#
# Copyright 2007 - 2014 by the NuGrid Team.
# All rights reserved. See LICENSE.
#

"""
h5T.py V1.1 created by Will Hillary

    email: will.hillary@gmail.com

    Date: December 18 2009

Modified to V1.1 by Daniel Conti

    Email: hoshi@uvic.ca

    Date: December 2010

Modified to V1.2 by Samuel Jones and Luke Siemens

    Email: swjones@uvic.ca, lsiemens@uvic.ca

    Date: March 2014

Modified by Alexander Heger for some python3 compatibility

    Email: alexander.heger@monash.edu

    Date May 2016

This is a simple python terminal interface for hdf5 files generated
using NuGrid Stellar evolution code.  It provides as simple as possible
interface by allowing for intuitive, but sparse commands.  Refer to the
users manual and docstring below for more information.

"""
from __future__ import print_function
from __future__ import absolute_import

from builtins import map
from builtins import str
from builtins import range
from builtins import sorted
from past.builtins import basestring

import os
import numpy as np
import gc
import threading
import time
import glob
import sys
import bisect

try:
    from .ascii_table import *
except ImportError:
    print('No module ascii_table')
    print('Please checkout ascii_table.py svn://forum.astro.keele.ac.uk/utils/pylib and add to python path')
try:
    import h5py as mrT
except ImportError :
    os.system('HDF5_DISABLE_VERSION_CHECK=1')
    import h5py as mrT

class Files(threading.Thread):
    '''
    Changes made by Daniel Conti

    A. Txtfile preprocessor: When the Files class is run in a directory
       for the first time, a preprocessor file is written, in the same
       directory.  When loading in a directory of files, this class
       opens up every file to find the ages and the cycles.  This file
       lets the class work faster, by writing the ages and cycles
       contained within each file to the preprocessor file.  This allows
       the class to not open every h5 file in the directory.  There is
       one downside to this however, This Files class has a list of
       h5File threads, if the a prepocessor file is read, these threads
       do not get started.  So to solve this, whenever a method needs
       some data from one of these h5File threads, it then starts the
       thread and adds a True boolean to a list of booleans (that
       correspond to the list of files) of if that particular thread
       has been called.  Currently, to the authors knowledge, this has
       been taken care of.  If in the future someone adds a method that
       works with any of these threads in self.h5s.  For example, if the
       command:

       >>> table=h5.Table

       where h5 is one of these h5File threads, the future author should
       have this aswell:

       >>> if not self.h5sStarted[self.h5s.index(h5)]
       >>>     h5.start()
       >>>     h5.join()
       >>>     table=h5.Table
       >>>     self.h5sStarted[self.h5s.index(h5)]=True
       >>> else:
       >>>     table=h5.Table

       This will run (start) the required thread and adds a True
       boolean, so we only ever need to start a thread once.

    B. Additional functionality to the get method: The get method can
       now be called with the name of the isotope, instead of the
       isotope and 'iso_massf'.  For example if we wanted the abundance
       of Heleium 3 from cycle 20, we would call:

       >>> get(20,'He-3')

       This would yield identical results to

       >>> get(20,'iso_massf','He-3')

       Also this author found and solved numerous bugs

    C. Added a pattern argument to the init method, A user can use a
       Unix, like regular expression to find files in a directory.

    D. Added a findCycle method, that allows the class to find the
       closest Cycle to what the user inputs, for example if the user
       imputs a cycle that DNE, it will find the next closest one.

    Assumptions
    ===========

    The user must have the module ascii_table and in the python path.
    It can be checked out from the svn svn://forum.astro.keele.ac.uk/ut$
    To use the preprocessor, the user need write access to the folder
    he is working in.

    Cycle numbers are 10 digits long

    Cycle numbers are allways whole numbers

    The preprocessor file will become unstable (ie not return accurite
    results), if any of the H5 files, internally, change their cycles
    or ages.

    If any H5 files are added and removed, the method will realize that
    there are additional or missing files and will write a new
    preprocessor file.

    If any of the files are renamed, The program realize that files were
    renamed and will write a new preprocessor file.

    If a new selection of files are slected, The program realize this
    and will write a new preprocessor file.

    If a prprocessor file is removed, The program realize this and will
    write a new preprocessor file.

    In each file name the cycle number is the list of numbers
    .(surf/out.h5)

    '''
    preprocName='h5Preproc.txt' #filename of the preprocessor file
    preprocExists=False         #Boolean of if the preprocessor file exists
    h5files = []
    h5s = []                    #list of h5File threads, for each of the file names
    cycles = []
    ages =     []
    hattrs = []
    cattrs = []
    Tables = []
    dcols =     []
    filepaths = []
    isotopes = []
    isomeric_states = []
    A = []
    Z = []
    groundState=1      #The identifyer for ground state
    excitedState=2     #The
    isomerDelimiter='m'#The chartacter string thet we use to seperate the isomer
                       #identifer from the rest of the isotope name
    h5sStarted=[]      #A list of booleans of if the thread in h5files has had its
                        #start and join methods called
    #    This is a list of isos to match to.
    isos = ['Neutron','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al',
    'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',
    'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc',
    'Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te', 'I','Xe','Cs','Ba','La','Ce',
    'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',
    'W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At']

    def findCycle(self, cycNum):
        '''
        Method that looks through the self.cycles and returns the
        nearest cycle:

        Parameters
        ----------
        cycNum : int
            int of the cycle desired cycle.

        '''
        cycNum=int(cycNum)
        i=0

        while i < len(self.cycles):
            if cycNum < int(self.cycles[i]):
                break
            i+=1

        if i ==0:
            return self.cycles[i]
        elif i == len(self.cycles):
            return self.cycles[i-1]
        lower=int(self.cycles[i-1])
        higher=int(self.cycles[i])

        if higher- cycNum >= cycNum-lower:
            return self.cycles[i-1]
        else:
            return self.cycles[i]




    #    Upon initialization of an h5fileholder, the h5 files are defined (in their own wrapper-see below)
    #    and gathers some important data from the files.
    def __init__(self ,path='.', fName=None, pattern='*', rewrite=False, verbose=False):
        '''
        Init method

        Parameters
        ----------
        path : string, optional
            The path you would like to work in.  The default is '.'.
        fName : string or list, optional
            A File name, or a list of filenames that you would like to
            work with.  The default is None.
        pattern : string, optional
            A string that contains a particular pattern that the user is
            working with. For example '1001.out' would make this method
            only read files that contain that string.  Note one can
            imput any regular expression (ie shell-style wildcards) into
            this, as one would do in A Unix shell.  The default is '*'.
        rewrite : boolean, optional
            If True, forces a rewrite of the preprocessor file.  The
            default is False.
        verbose : boolean, optional
            Set to True for more output.

        '''
        threading.Thread.__init__(self)
        self.h5files = []
        self.h5s = []
        self.cycles = []
        self.ages =     []
        self.hattrs = []
        self.cattrs = []
        self.Tables = []
        self.dcols =     []
        self.filepaths = []
        self.isotopes = []
        self.elements = []
        self.isomeric_states = []
        self.A = []
        self.Z = []
        self.h5sStarted=[]
        self.verbose=verbose

        if fName==None and pattern=='*':
            self.filename = path
            self.files = os.listdir(self.filename)

            for fil in self.files:
                if os.path.isfile(self.filename + os.sep + fil) and fil.count('.h5'):
                    self.h5files.append(self.filename + os.sep + fil)
        elif pattern=='*':
            self.filename = path
            if self.filename[-1] == os.sep:
                self.filename = self.filename[:-1]
            self.files = os.listdir(self.filename)
            self.h5files = []
            temps = fName
            if type(temps).__name__ != 'list':
                temps = [temps]
            for arg in temps:
                self.h5files.append(self.filename + os.sep + arg)
        elif fName==None:
            self.filename = path
            if not pattern.endswith('*'):
                pattern=pattern+'*'
            if not pattern.startswith('*'):
                pattern='*'+pattern
            self.files=glob.glob(self.filename + os.sep +pattern)
            for i in range(len(self.files)):
                self.files[i]=self.files[i].split(os.sep)[-1]

            for fil in self.files:
                if os.path.isfile(self.filename + os.sep + fil) and fil.count('.h5'):
                    self.h5files.append(self.filename + os.sep + fil)

        #preprocessor stuff
        if self.preprocName.endswith('/'):
            preprocName = str(self.filename)+self.preprocName

        else:
            preprocName = str(self.filename)+'/'+self.preprocName

        self.preprocExists=os.path.exists(preprocName)
        if rewrite:
            self.preprocExists=False
        self.preprocExisted=self.preprocExists # to know whether or not the preproc was
                # there in the first place (SJONES)

        if self.h5files == []:
            print('There are no h5Files in: ', self.filename, 'please try a better folder.')
        else:
            print('Searching files, please wait.......')
            for i in range(len(self.h5files)):
                self.h5sStarted.append(False)
            # SJONES: here, i also now pass whether preprocExisted:
            #self.h5s.append(h5File(self.h5files[0],True, True))
            self.h5s.append(h5File(self.h5files[0],True, True, self.preprocExisted))

            self.h5s[0].start()
            self.h5s[0].join()
            self.h5sStarted[0]=True
            self.start()


            self.join()
            #print self.h5s

            #self.connect(self.h5s[-1], qc.SIGNAL('finished()'), self.continue_h5s)

    def run(self):

        if not self.preprocExisted:
            self.cycles.extend(self.h5s[0].cycle)
            self.ages.extend(self.h5s[0].age)
        self.hattrs = self.h5s[0].hattr
        self.cattrs =   self.h5s[0].cattr
        self.Tables    =   self.h5s[0].Table
        self.dcols     =   self.h5s[0].dcol
        self.cycle_header = self.h5s[0].cycle_header    #This string handles the name of the cycle
        try:
            self.A = self.h5s[0].A[0]
        except IndexError:
            print("Sorry, there is no A vector. This can cause problems for reading abundances. Continue...")
        try:
            self.Z = self.h5s[0].Z[0]
        except IndexError:
            print("Sorry, there is no Z vector. This can cause problems for reading abundances. Continue... ")
        try:
            self.isomeric_states = self.h5s[0].isomeric_state[0]
        except IndexError:
            print("Sorry, there is no isomeric state vector. Continue...")

        new = self.h5s[0].new    #    This boolean handles the changes to the cycle nomenclature format

        if self.filename.endswith(os.sep):
            b = str(self.filename)+str(self.preprocName)
        else:
            b = str(self.filename)+os.sep+str(self.preprocName)

        if self.preprocExists:
            preprocTable=ascii_table(self.preprocName,self.filename)
            if int(preprocTable.hattrs[0])<len(self.h5files):
                self.preprocExists=False
                print('A File was added, rewriteing preprocessor file')

            if self.preprocExists:
                for i in range(len(self.h5files)):
                    if os.path.basename(self.h5files[i])+'-cyc' not in preprocTable.dcols and self.preprocExists:
                        print('A File was renamed, rewriteing preprocessor file')
                        if self.verbose:
                            print(preprocTable.dcols[i], os.path.basename(self.h5files[i])+'-cyc')
                        self.preprocExists=False

        if not self.preprocExists and os.path.exists(b):
            os.system('rm '+b)

# create list of isotopes stored in this h5 file
        try:
            for x in range(len(self.Tables[0])):
                if self.isomeric_states[x] ==1:
                    self.isotopes.append(self.isos[int(self.Tables[1][x])]+'-' +str(int(self.Tables[0][x])))
                else:
                    self.isotopes.append(self.isos[int(self.Tables[1][x])]+'-' +str(int(self.Tables[0][x]))+\
                    self.isomerDelimiter+str(self.isomeric_states[x]-1))
        except IndexError:
            print('This file does not contain any tables.  Isotopic data must be contained elsewhere.')
        t1 = time.time()

# create list of elements stored in this h5 file
        asaved=''
        for atmp in self.Z:
            if atmp != asaved:
                self.elements.append(self.isos[int(atmp)])
                asaved=atmp

        thread_list = []
        if self.verbose:
            print(self.h5s[0].filename)
        for x in range(len(self.h5files)-1):
            # SJONES: here, i also now pass whether preprocExisted:
            #thread_list.append(self.h5s.append(h5File(self.h5files[x+1],False, new)))
            thread_list.append(self.h5s.append(h5File(self.h5files[x+1],False, new, self.preprocExisted)))
            if self.verbose:
                print(self.h5s[x+1].filename)
            if not self.preprocExists:
                self.h5s[x+1].start()
                self.h5s[x+1].join()
                self.h5sStarted[x+1]=True
            #self.h5s[x+1].run()
            #print('Active '+ str(threading.active_count()))

        #    Wait until the threads are done.
        #for thread in self.h5s:
        #    print(thread)
        #    thread.join()

        if not self.preprocExists:
            for x in range(len(self.h5files)-1):
                self.cycles.extend(self.h5s[x+1].cycle)
                self.ages.extend(self.h5s[x+1].age)
            header=[str(len(self.h5files)),'This is a preprocessor file for the directory: '+str(self.filename),\
            'At the time of the creation of this file there were '+str(len(self.h5files))+\
            ' h5 files.']

            try:
                self.cycles = sorted(self.cycles, key=int)
            except TypeError:
                print("There was a problem sorting the cycles.\nYou may have problems later.\nPlease consider reloading(h5T) and trying again.\n")

            try:
                self.ages = sorted(self.ages, key=int)
            except TypeError:
                None

            print('Writing preprocessor files')
            data=[]
            dcols=[]
            length=0
            for i in range(len(self.h5s)):
                print(self.h5s[i].filename.rpartition('/')[2])
                dcols.append(os.path.basename(self.h5s[i].filename)+'-cyc')
                dcols.append(os.path.basename(self.h5s[i].filename)+'-age')
                data.append(self.h5s[i].cycle)
                data.append(self.h5s[i].age)
                if len(self.h5s[i].cycle)>length:
                    length=len(self.h5s[i].cycle)
                if len(self.h5s[i].age)>length:
                    length=len(self.h5s[i].age)

            for i in range(len(data)):
                for j in range(length-len(data[i])):
                    data[i].append(3.14159265)

            write(self.preprocName,header,dcols,data,sldir=self.filename)

        else:
#            print 'Reading preprocessor files'
#            preprocTable=ascii_table(self.preprocName,self.filename)
#            for i in xrange(len(self.h5s)-1):
#                dat=preprocTable.get(os.path.basename(self.h5s[i+1].filename)+'-cyc')
#                dat1=[]
#                for j in xrange(len(dat)):
#                    if dat[j]!=3.14159265:
#                        dat1.append(dat[j])
#
#                dat=dat1
#                for j in xrange(len(dat)):
#                    dat[j]=str(int(dat[j]))
#                    for k in xrange(10-len(dat[j])):
#                        dat[j]='0'+dat[j]
#
#                for j in xrange(len(dat)):
#                    self.cycles.append(dat[j])
#                self.h5s[i+1].cycle=dat
#                dat=preprocTable.get(os.path.basename(self.h5s[i+1].filename)+'-age')
#                dat1=[]
#                for j in xrange(len(dat)):
#                    if dat[j]!=3.14159265:
#                        dat1.append(dat[j])
#                dat=dat1
#                self.h5s[i+1].age=dat
#                for j in xrange(len(dat)):
#                    self.ages.append(dat[j])
            print('Reading preprocessor files')
            preprocTable=ascii_table(self.preprocName,self.filename)
            for i in range(len(self.h5s)):
                dat=preprocTable.get(os.path.basename(self.h5s[i].filename)+'-cyc')
                dat1=[]
                for j in range(len(dat)):
                    if dat[j]!=3.14159265:
                        dat1.append(dat[j])

                dat=dat1
                for j in range(len(dat)):
                    dat[j]=str(int(dat[j]))
                    for k in range(10-len(dat[j])):
                        dat[j]='0'+dat[j]

                for j in range(len(dat)):
                    self.cycles.append(dat[j])
                self.h5s[i].cycle=dat
                dat=preprocTable.get(os.path.basename(self.h5s[i].filename) + '-age')
                dat1=[]
                for j in range(len(dat)):
                    if dat[j]!=3.14159265:
                        dat1.append(dat[j])
                dat=dat1
                self.h5s[i].age=dat
                for j in range(len(dat)):
                    self.ages.append(dat[j])
### end of new section ###
            try:
                self.cycles = sorted(self.cycles, key=int)
            except TypeError:
                print("There was a problem sorting the cycles.\nYou may have problems later.\nPlease consider reloading(h5T) and trying again.\n")

            try:
                self.ages = sorted(self.ages, key=int)
            except TypeError:
                None
        print('File search complete.')
        t2 = time.time()
        if self.verbose:
            print("Total duration is " + str(t2-t1) + " seconds.")
        return

    '''
    def startThreads(self,threads,IsConcurrent):
         """
         Method that starts the h5 file threads
         Input: Threads- a list of threads that need to be started
         """

         if (str(threads.__class__())==("<type 'list'>")):
                threads=[threads]
    '''

    # This function determines which cycle, which file, which storage mechanism (cattr or data) and returns it
    def get(self, cycle_list, dataitem=None, isotope=None, sparse=1):
        '''
        Get Data from HDF5 files.

        There are three ways to call this function

        1. get(dataitem)

            Fetches the datatiem for all cycles. If dataitem is a header
            attribute or list of attributes then the data is retured.
            If detaitem an individulal or list of column attributes,
            data columns or isotopes/elements the data is returned for
            all cycles.

        2. get(cycle_list, dataitem)

            Fetches the dataitem or list of dataitems for the cycle
            or list of cycles. The variable dataitems can contain column
            attributes, data columns, and isotopes/elemnts.

        3. get(cycle_list, dataitem, isotope)

            Fetches the dataitems like the seccond method except that
            one of the dataitems must be either "iso_massf" or "yps",
            and in the data returned "iso_massf" and "yps" are replaced
            with the data from the isotopes.  The isotopes must be in
            the form given by se.isotopes or se.elements.

        Parameters
        ----------
        cycle_list : list, integer or string
            If cycle_list is a list or string and all of the entries
            are header attributes then the attributes are returned.

            If cycle_list is a list or string of dataitems then the
            dataitems are fetched for all cycles.

            If cycle_list is a list, integer or string of cycle numbers
            then data is returned for those cycles.
        dataitem: list or string, optional
            If dataitem is not None then the data for each item is
            returned for the cycle or list of cycles. dataitem may be an
            individual or a mixed list of column attributes, column
            data or isotopes/elements. If dataitem is None then
            cycle_list must be a string.  The default is None.
        isotope: list or string, optional
            If one of the dataitems is "iso_massf" or "yps" then it is
            replaced with the data from the individual isotopes/elements
            listed in isotope.  The default is None.
        sparse : int
            Implements a sparsity factor on the fetched data i.e. only
            the i th cycle in cycle_list data is returned,
            where i = sparse.

        '''

#    Check out the inputs
        t1=time.time()
        isotopes_of_interest = []

        nested_list = False
        # if one of cycle_list, dataitem or isotope is given as a string convert it to a list
        if isinstance(cycle_list, basestring):
            cycle_list = [cycle_list]
        else:
            try:
                if len(cycle_list) == 1:
                    nested_list = True
            except TypeError:
                pass #leave nested_list as false
        if isinstance(dataitem, basestring):
            dataitem = [dataitem]
        if isinstance(isotope, basestring):
            isotope = [isotope]



        if dataitem==None and isotope==None:
            option_ind = 1
            dataitem = cycle_list

            if not any([item in self.hattrs for item in dataitem]):
                cycle_list = self.cycles
            else:
                first_file = mrT.File(self.h5s[0].filename,'r')
                dat = []
                # get all dataitems from header attributes
                for item in dataitem:
                    tmp = first_file.attrs.get(item, None)
                    try:
                        if len(tmp) == 1:
                            tmp = tmp[0]
                    except TypeError: #if a scaler is returned do nothing
                        pass
                    dat.append(tmp)
                # if only one header attribute is required dont return as a list
                if (len(dat) == 1) and (not nested_list):
                    dat = dat[0]
                first_file.close()
                return dat
            if any([item.split('-')[0] in self.isos for item in dataitem]):
                return self.get(cycle_list,dataitem,sparse=sparse)
        elif isotope==None:
            option_ind = 2
            cycle_list = cycle_list
            dataitem = dataitem
            # if one dataitem is given as a string convert it to a list
            if isinstance(dataitem, basestring):
                dataitem = [dataitem]
            new_dataitem = []
            new_isotopes = []
            for item in dataitem:
                if item.split('-')[0] in self.isos:
                    new_isotopes.append(item)
                else:
                    new_dataitem.append(item)
            if len(new_isotopes) != 0:
                tmp = []
                try:
                    tmp = self.get(cycle_list,new_dataitem + ['iso_massf'],new_isotopes,sparse=sparse)
                except: # in some old se files there maybe still yps as the name for the abundance arrays
                    tmp =  self.get(cycle_list,new_dataitem + ['yps'],new_isotopes,sparse=sparse)
                # modify the dat list so dat is structured like dataitems
                dat = []
                #make sure tmp containes the data as a list of cycles
                if isinstance(cycle_list, basestring):
                    tmp = [tmp]
                else:
                    try:
                        if len(cycle_list) == 1:
                            tmp = [tmp]
                    except TypeError:
                        tmp = [tmp]
                for cyc in tmp:
                    temp_dataitem = []
                    for item in dataitem:
                        if item in new_dataitem:
                            temp_dataitem.append(cyc[new_dataitem.index(item)])
                        else:
                            if len(new_dataitem) == 0:
                                temp_dataitem = cyc
                            else:
                                if len(new_isotopes) == 1:
                                    temp_dataitem.append(cyc[-1])
                                else:
                                    temp_dataitem.append(cyc[-1][new_isotopes.index(item)])
                    dat.append(temp_dataitem)
                if (len(dat) == 1) and (not nested_list):
                    dat = dat[0]
                return dat
        else:
# there is an implicite rule here that if you want 2D arrays you have
# to give 3 args, or, in other words you have to give a cycle or cycle
# array; there is no good reason for that, except the programmers
# laziness
            option_ind = 3
            cycle_list = cycle_list
            dataitem = dataitem
            isotopes_of_interest = isotope
# we need to find out the shellnb to know if any yps array may just be
# a one row array, as - for example- in the surf.h5 files
            # SJONES: I think here we only need to look at the first shellnb(!)
            #shellnb=self.get(cycle_list,'shellnb')
            try: #check if cycle_list is not a list
                cycle_list[0]
            except (TypeError,IndexError):
                cycle_list = [cycle_list]
            shellnb=self.get(cycle_list[0],'shellnb')

        if sparse <1:
            sparse=1

        #    Just in case the user inputs integers
        try:
            for x in range(len(cycle_list)):
                cycle_list[x] = str(cycle_list[x])
        except TypeError:
            cycle_list = [str(cycle_list)]

        if option_ind != 1:

            try: #if it is a single cycle make sure its formatted correctly
                if cycle_list.isdigit():
                    cycle_list = [cycle_list]
                    for cycle in cycle_list:
                        if len(cycle) != len(self.cycles[0]):
                                #print "a"
                            diff = len(self.cycles[0])-len(cycle)
                            OO = ''
                            while diff >=1:
                                OO+='0'

                            cycle = OO+cycle

            except AttributeError: ##if it is a list of cycles make sure its formatted correctly
                if cycle_list[0].isdigit():

                    for x in range(len(cycle_list)):
                        if len(str(cycle_list[x])) != len(str(self.cycles[0])):
                            #print "b"
                            diff = len(str(self.cycles[0]))-len(str(cycle_list[x]))

                            OO = ''
                            while diff >=1:
                                OO+='0'
                                diff-=1

                            try:
                                cycle_list[x] = OO+cycle_list[x]
                            except TypeError:
                                cycle_list[0] = OO+cycle_list[0]

        dat = []
        cycle_list.sort()

        cyclelist=np.array(list(map(int, cycle_list)))

        # cycles_requested is a list of indices from cyclelist
        # The index of the larges and smallest indices should be stored
        # in sorted order. As new requests are made if the requests
        # border or over lap then only keep the index of the larges and
        # smallest indices.
        cycles_requested = []

        # Sometimes bad data or last restart.h5 files contain no cycles,
        # causing the code to crash. Do a simple try/except here:
        file_min=[]
        file_max=[]
        try:
            for h5 in self.h5s:
                file_min.append(int(h5.cycle[0]))
                file_max.append(int(h5.cycle[-1]))
        except IndexError:
            print('File '+h5.filename+' contains no data, please remove or rename it')
            print('Once the file has been removed or renamed, the preprocessor file must be re-written. Do this by either removing the file h5Preproc.txt from the data directory or by invoking the se instance with rewrite=True')
            print('At present, h5T cannot check for empty files since the overhead using the mounted VOSpace would be too great.')
            raise IOError('Cycle-less file encountered')
        file_min.sort()
        file_max.sort()

        for h5 in self.h5s:
            #initalize file metadata
            min_file = int(h5.cycle[0])
            max_file = int(h5.cycle[-1])
            min_list = int(cyclelist[0])
            max_list = int(cyclelist[-1])
            index_min = None #if None start at begining
            index_max = None #if None finish at end

            # SJONES Now we need to add the case that the set only contains one file:
            if len(file_min) == 1:
                min_file = min_list - 1
                max_file = max_list + 1
            else:
                file_index = file_min.index(min_file)
                if file_index == 0:
                    if min_list - 1 < min_file:
                        min_file = min_list - 1
                    max_file = (file_min[file_index + 1] + max_file)//2
                elif file_index == len(file_min) - 1:
                    min_file = (file_max[file_index - 1] + min_file)//2 + 1
                    if max_list + 1 > max_file:
                        max_file = max_list + 1
                else:
                    min_file = (file_max[file_index - 1] + min_file)//2 + 1
                    max_file = (file_min[file_index + 1] + max_file)//2

            # calculate the left and right limits of the intersection
            # of the lists h5.cycle and cyclelist
            if (max_list < min_file) or (max_file < min_list):
            # the lists do not intersect
                continue
            elif (min_list <= min_file) and (max_file <= max_list):
                # all of h5.cycle is within cyclelist
                index_min = bisect.bisect_left(cyclelist, min_file)
                index_max = bisect.bisect_right(cyclelist, max_file)
            elif (min_file <= min_list) and (max_list <= max_file):
                # all of cyclelist is within h5.cycle
                index_min = None
                index_max = None
            else:
                if min_list > min_file:
                # cyclelist overlaps the right edge of h5.cycle
                    index_min = None
                    index_max = bisect.bisect_right(cyclelist, max_file)
                else:
                    # cyclelist overlaps the left edge of h5.cylce
                    index_min = bisect.bisect_left(cyclelist, min_file)
                    index_max = None

            # maintin list of all requested cycles by keeping trak of
            # the maximum and minimum indices
            imin = index_min
            if index_min == None:
                imin = 0

            imax = index_max
            if index_max == None:
                imax = len(cyclelist)

            request_min = bisect.bisect_left(cycles_requested, imin)
            request_max = bisect.bisect_right(cycles_requested, imax)

            # if the new request overlabs older request remove them
            del cycles_requested[request_min:request_max]
            if ((request_max-request_min) % 2) ==1:
                # new and old request overlaped on one edge only
                if request_min % 2 == 0:
                    # add new starting index
                    cycles_requested.insert(request_min, imin)
                else:
                    # add new ending index
                    cycles_requested.insert(request_min, imax)
            else:
                # new and old requests overlaped on two edges
                if request_min % 2 == 0:
                    # old request was contained with in new request
                    cycles_requested.insert(request_min, imin)
                    cycles_requested.insert(request_min + 1, imax)
                else:
                    # new request wat contained within old request
                    pass

            if not self.h5sStarted[self.h5s.index(h5)]:
                h5.start()
                h5.join()
                temp = h5.fetch_data_sam(dataitem,cycle_list[index_min:index_max],len(cycle_list),len(dat))
                self.h5sStarted[self.h5s.index(h5)]=True
            else:
                temp = h5.fetch_data_sam(dataitem,cycle_list[index_min:index_max],len(cycle_list),len(dat))

            temp_dat = []
            for temp_num, temp_cycle in enumerate(temp):
                temp_dataforcycle = []
                for dataitem_num, temp_dataitem in enumerate(temp_cycle):
                    # identify what cycle the temp data was collected from
                    temp_dataitem=self.red_dim(temp_dataitem)
#                    if option_ind == 3 and isotopes_of_interest != []:
                    if (dataitem[dataitem_num] == 'iso_massf' or dataitem[dataitem_num] == 'yps') and isotopes_of_interest != []:
                    #    Figure out the index
                        index = []
                        iso_tmp = []

                        if 'iso' in dataitem[dataitem_num]: #if we are looking at an isotope
                            iso_tmp = self.isotopes
                        else:
                            iso_tmp = self.elements

                        for iso in isotopes_of_interest: #finds the location of the isotope
                            x = iso_tmp.index(iso)
                            index.append(x)

                        if index == []:
                            # if none of the isotopes of interest are found
                            # then the index defaults to [0], so that the loop
                            # will still try to acess the data in t.
                            index = [0]
                        islist=True

                        if len(cycle_list)==1:
                            islist=False

#                        shellnb_index = 0
#                        if index_min == None:
#                            shellnb_index = temp_num
#                        else:
#                            shellnb_index = index_min + temp_num
                        temp_multicyc = []
                        for i in index:
#                            if islist:
#                                if shellnb[shellnb_index] == 1:    # again take care of 1-row 2D arrays
                                if shellnb == 1:    # again take care of 1-row 2D arrays
                                    temp_multicyc.append(temp_dataitem[i])
                                else:
                                    temp_multicyc.append(temp_dataitem[:,i])
#                            else:
#                                if shellnb == 1:    # again take care of 1-row 2D arrays
#                                    temp_multicyc.append(temp_dataitem[i])
#                                else:
#                                    temp_multicyc.append(temp_dataitem[:,i])
                        if len(temp_multicyc) == 1: # agian take care of 1-row arrays
                            temp_multicyc = temp_multicyc[0]
                        temp_dataitem = temp_multicyc
                    temp_dataforcycle.append(temp_dataitem)
                if len(temp_dataforcycle) == 1: # agian take care of 1-row arrays
                    temp_dataforcycle = temp_dataforcycle[0]
                # Now add the information to the list we pass back
                temp_dat.append(temp_dataforcycle)
            # calculate the proper insertion point for the data colected from
            # the file h5 in self.h5s
            insert_pnt = 0
            if index_min is not None: #alex: in py2: x < None == False
                for i in range(len(cycles_requested)):
                    if i % 2 == 1:
                        if cycles_requested[i] < index_min:
                            insert_pnt += cycles_requested[i] - cycles_requested[i-1]
                        elif cycles_requested[i - 1] < index_min:
                            insert_pnt += index_min - cycles_requested[i - 1]
            # insert the cycle data from the current file into the apropiat place
            # in the output data.
            dat[insert_pnt:insert_pnt] = temp_dat

        #check if cycles were not requested from the file
# SJONES comment
#        missing_cycles = np.array([])
#        if len(cycles_requested) != 2:
#            if len(cycles_requested) == 0:
#                missing_cycles = np.array([cycle_list])
#            else:
#                cycles_requested = [None] + cycles_requested + [None]
#                for i in xrange(0, len(cycles_requested), 2):
#                    min = cycles_requested[i]
#                    max = cycles_requested[i + 1]
#                    missing_cycles = np.append(missing_cycles, cycle_list[min:max])
#            print "The requested cycles: " + str(missing_cycles) + " are not available in this data set"
#        elif (cycles_requested[0] != 0) or (cycles_requested[1] != len(cyclelist)):
#            min = cycles_requested[0]
#            max = cycles_requested[1]
#            missing_cycles = np.append(missing_cycles, cycle_list[0:min])
#            missing_cycles = np.append(missing_cycles, cycle_list[max:])
#            print "The requested cycles: " + str(missing_cycles) + " are not available in this data set"

        if len(dat) < 2 and option_ind != 3 and (not nested_list):
            try:
                dat = dat[0]
            except IndexError:
                None
            except TypeError:
                None
        try:
            if len(dat) < 2 and isotopes_of_interest != []:
                dat = dat[0]
        except TypeError:
            None
        except IndexError:
            None
        t2=time.time()
        return dat

    #    uses the index information to build list of isos from tables A,Z
    def fetch_isos(self):
        isos = []
        try:
            for x in range(len(self.Tables[1])):
                isos.append([self.isos[int(self.Tables[1][x])], self.Tables[0][x]])
        except IndexError:
            None
        return isos


    def red_dim(self, array):
        """
        This function reduces the dimensions of an array until it is
        no longer of length 1.

        """
        while isinstance(array, list) == True or \
        isinstance(array, np.ndarray) == True:
            try:
                if len(array) == 1:
                    array = array[0]
                else:
                    break
            except:
                break

        return array

#    This wrapper class allows some automated activity when an h5file is initialized.
#    upon inmitialization the h5 file is opened and certain bits of data is read.
#    This class also interacts with the h5fileholder class to access needed data.

class h5File(threading.Thread):
    h5 = None
    filename = None
    cycle = []
    age = []
    dcol = []
    data = []
    skipped_nodes = 0
    ver = ''
    classname = ''
    hattr = []
    cattr=[]
    Table = []
    isomeric_state = []
    A = []
    Z = []
    new = True



    #    Initialize the class
    # SJONES added whether proprocExisted here:
    #def __init__(self, filepath,deep_search, new):
    def __init__(self, filepath,deep_search, new, waspreprocthere):
        threading.Thread.__init__(self)
        #    Instantiate
        self.h5 = None
        self.filename = None
        self.cycle = []
        self.age = []
        self.dcol = []
        self.data = []
        self.skipped_nodes = 0
        self.ver = ''
        self.classname = ''
        self.hattr = []
        self.cattr=[]
        self.Table = []
        self.isomeric_state = []
        self.A = []
        self.Z = []
        self.new = True
        self.preprocExisted = waspreprocthere

        #    Build
        self.new = new
        self.filename = filepath
        self.deep_search = deep_search


        if self.new:
            self.cycle_header = 'cycle'
        else:
            self.cycle_header = 'cycle-'

    def run(self):

        if self.deep_search:
            #self.search_deep()
            self.search_deep_sam()
        else:
            if not self.preprocExisted:
                self.search_shallow()
        try:
            self.h5.close()
        except:
            None

        return None

    #    Fetches a single category of information
    def fetch_data_one(self,dataitem,cycle):
        self.h5 = mrT.File(self.filename,'r')

        try:
            data = self.h5[self.cycle_header+str(cycle)]['SE_DATASET'][dataitem]
        except ValueError:
            try:
                data = self.h5[self.cycle_header+str(cycle)].attrs.get(dataitem, None)
            except TypeError:
                data = self.h5[self.cycle_header+str(cycle)][dataitem]

        try:
            while data.shape[0] < 2:
                data = data[0]
        except (IndexError, AttributeError):
            None


        self.h5.close()
        return data

    def fetch_data_sam(self, dataitemlist, cyclelist, total_list_size,
                       current_list_size):
        #like fetch_data_one but it accepts a list of cycles and a list of data tiems
        quiet = True # flag for supressing info about replacing missing cycles. Should
                     # be integrated better eventually.
        self.h5 = mrT.File(self.filename,'r')
        is_error=False
        missing_cycles = np.array([])
        actual_cycles_flt=np.array([float(cyc) for cyc in self.cycle])
        data=[]
        old_percent = 0
        for cycle in cyclelist:
            current_list_size += 1
            percent = int(current_list_size*100/total_list_size)
            if percent >= old_percent + 5:
                sys.stdout.flush()
                sys.stdout.write("\r reading "+ str(dataitemlist) + "...%d%%" % percent)
                old_percent = percent

            try:
                dataitem_data = []
                for dataitem in dataitemlist:
                    try:
                        dataitem_data.append(self.h5[self.cycle_header+str(cycle)]['SE_DATASET'][dataitem])
                    except ValueError:
                        try:
                            dataitem_data.append(self.h5[self.cycle_header+str(cycle)].attrs.get(dataitem, None))
                        except TypeError:
                            dataitem_data.append(self.h5[self.cycle_header+str(cycle)][dataitem])
                data.append(dataitem_data)
            except KeyError:
                if not is_error:
                    is_error = True
                missing_cycles = np.append(missing_cycles, str(cycle))
                nearestidx=np.abs(actual_cycles_flt-float(cycle)).argmin()
                cycledummy=self.cycle[nearestidx]
                dataitem_data = []
                for dataitem in dataitemlist:
                    try:
                        dataitem_data.append(self.h5[self.cycle_header+str(cycledummy)]['SE_DATASET'][dataitem])
                    except ValueError:
                        try:
                            dataitem_data.append(self.h5[self.cycle_header+str(cycledummy)].attrs.get(dataitem, None))
                        except TypeError:
                            dataitem_data.append(self.h5[self.cycle_header+str(cycledummy)][dataitem])
                data.append(dataitem_data)
#                continue

        if is_error:
            if not quiet:
                print("The requested cycles: " + str(missing_cycles) + " are not available in this data set. They have been replaced with the nearest available data.")

        self.h5.close()
        return data


    #    The typical search algirthm when a h5file class is initialized
    def search_shallow(self):
        self.h5 = mrT.File(self.filename,'r')
        temp = list(self.h5.keys())
        self.cycle=[str(k).replace('cycle','').replace('-','') for k in temp if 'cyc' in k]
        try:
            self.age=[float(self.h5[str(asd)].attrs.get("age",None)) for asd in temp if 'cyc' in asd]
        except ValueError:
            self.age=[self.h5[str(asd)].attrs.get("age",None)[0] for asd in temp if 'cyc' in asd]

# Old code that was getting cycles and ages is commented here below;
# now instead the cycles are got much more effectively (above, and subject to testing!)
# but the ages are got in the same way as before, but ONLY if the preprocessor file
# did not exist at the time of invocation of the Files instance.

#        for te in temp:
#            if te[0] == 'c':
#                if te[5:].isdigit():
#                    self.cycle.append(str(te[5:]))
#
#                    try:
#                        self.age.append(self.h5[te].attrs.get("age",None)[0])
#                    except TypeError:
#                        self.age.append(self.h5[te].attrs.get("age",None))
#                else:
#                    self.cycle.append(str(te.split('-')[1]))
#                    try:
#                        self.age.append(self.h5[te].attrs.get("age", None)[0])
#                    except TypeError:
#                        self.age.append(self.h5[te].attrs.get("age",None))
#                        self.cycle.sort()
#                        self.age.sort()



    def search_deep(self):
        self.h5 = mrT.File(self.filename,'r')
        temp = list(self.h5.keys())

        #    Handles the change in cycle nomenclature
        self.new = True
        for te in temp:
            if te.count('-'):
                self.new = False
                break

        if self.new:
            self.cycle_header = 'cycle'
            for te in temp:
                if te[0] == 'c':
                    if te[5:].isdigit():
                        self.cycle.append(str(te[5:]))
                        try:
                            self.age.append(self.h5[te].attrs.get("age",None)[0])
                        except TypeError:
                            self.age.append(self.h5[te].attrs.get("age",None))
                    else:
                        self.isomeric_state.append(self.h5[te]['data'])
                else:

                    obj = self.h5[te].__iter__()
                    if str(te).count('A'):
                        holder = []
                        for ob in obj:
                            holder.append(ob[0])
                        self.Table.append(holder)
                        self.A.append(holder)
                    elif str(te).count('Z'):
                        holder = []
                        for ob in obj:
                            holder.append(ob[0])
                        self.Table.append(holder)
                        self.Z.append(holder)
                    else:
                        holder = []
                        for ob in obj:
                            holder.append(ob[0])
                        self.Table.append(holder)
                        self.isomeric_state.append(holder)

            try:
                temp =  self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).__getitem__('SE_DATASET').dtype.__str__().split(',')
            except ValueError:
                temp =  self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).dtype.__str__().split(',')

        else:
            self.cycle_header = 'cycle-'
            for te in temp:
                try:
                    self.cycle.append(str(te.split('-')[1]))

                    try:
                        self.age.append(self.h5[te].attrs.get("age", None)[0])
                    except TypeError:
                        self.age.append(self.h5[te].attrs.get("age", None))
                except IndexError:

                    obj =  self.h5[te].__iter__()

                    if str(te).count('A'):
                        holder = []
                        for ob in obj:
                            holder.append(ob[0])
                        self.Table.append(holder)
                        self.A.append(holder)
                    elif str(te).count('Z'):
                        holder = []
                        for ob in obj:
                            holder.append(ob[0])
                        self.Table.append(holder)
                        self.Z.append(holder)
                    else:
                        holder = []
                        for ob in obj:
                            holder.append(ob[0])
                        self.Table.append(holder)
                        self.isomeric_state.append(holder)

            self.cycle.sort()

        # This is kind of stupid, but I have not found a way to access this information directly.

            try:
                temp =  self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).__getitem__('SE_DATASET').dtype.__str__().split(',')
            except ValueError:
                temp =  self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).dtype.__str__().split(',')


        for tem in temp:
            if tem.count('<') ==0:
                try:
                    self.dcol.append(tem.split('\'')[1])
                except IndexError:
                    None

        attrs = self.h5.attrs
        for at in attrs:
            self.hattr.append(at)
        self.cattr = list(self.h5[self.cycle_header+str(self.cycle[0])].attrs.keys())

        table = []
        grp = self.h5[self.cycle_header+str(self.cycle[0])]
        for gr in grp:
            try:
                table.append(float(gr[0]))
            except ValueError:
                None

        self.h5.close()
        return None

    def search_deep_sam(self):
        self.h5 = mrT.File(self.filename,'r')
        temp = list(self.h5.keys())

        # SJONES attempt to speed this up and not search through the file:

        cyclenames=[t for t in temp if 'cyc' in t]

        if '-' in cyclenames[0]:
            cycle_header='cycle='
        else:
            cycle_header='cycle'

        self.cycle=[str(k).replace('cycle','').replace('-','') for k in temp if 'cyc' in k]
        # SJONES: Frustratingly, it will take more work to just take the age from the
        # preprocessor file, where it lives... but for now, we just do this for the first
        # file : (
        if self.preprocExisted:
            self.age=None
        else:
            try:
                self.age=[float(self.h5[str(asd)].attrs.get("age",None)) for asd in temp if 'cyc' in asd]
            except ValueError:
                self.age=[self.h5[str(asd)].attrs.get("age",None)[0] for asd in temp if 'cyc' in asd]

        others=[]
        for t in temp:
            if 'cyc' not in t:
                others.append(t)

        for o in others:
            obj=self.h5[o].__iter__()
            if str(o)=='A':
                holder = []
                for ob in obj:
                    holder.append(ob[0])
                self.Table.append(holder)
                self.A.append(holder)
            elif str(o)=='Z':
                holder = []
                for ob in obj:
                    holder.append(ob[0])
                self.Table.append(holder)
                self.Z.append(holder)
            else:
                holder = []
                for ob in obj:
                     holder.append(ob[0])
                self.Table.append(holder)
                self.isomeric_state.append(holder)

        self.cycle.sort()

        # This is kind of stupid, but I have not found a way to access this information directly.

        try:
            temp =  self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).__getitem__('SE_DATASET').dtype.__str__().split(',')
        except ValueError:
            temp =  self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).dtype.__str__().split(',')

        for tem in temp:
            if tem.count('<') ==0:
                try:
                    self.dcol.append(tem.split('\'')[1])
                except IndexError:
                    None

        attrs = self.h5.attrs
        for at in attrs:
            self.hattr.append(at)
        self.cattr = list(self.h5[self.cycle_header+str(self.cycle[0])].attrs.keys())

        table = []
        grp = self.h5[self.cycle_header+str(self.cycle[0])]
        for gr in grp:
            try:
                table.append(float(gr[0]))
            except ValueError:
                None

        self.h5.close()
        return None
