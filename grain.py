'''
grain.py is a collection of routines to analyze presolar grain data
(and can probably be extended to observational data at a later stage).

This class load the current version of the presolar grain database for
further processing.  A private databse can be given as well as
described later.  Several routines (see below and NuGrid book) can be
used to filter, plot, and retrieve grain data.  The presolar grain
database is supported by the group at Washington University, mainly
Frank Gyngard.  The database can be found at
http://presolar.wustl.edu/PGD/Presolar_Grain_Database.html

Important note: This script assumes that you have a full SVN tree
checked out (or actually, that you have at least the utils folder and
the validation folder on the same level checked out.

Usage of these tools
====================
For questions, bug reports, please contact trappitsch@uchicago.edu
Reto Trappitsch for the NuGrid collaboration

Developing notes
================
data header should be safed lowercase always -> use .lower() argument
also for comparison w/ user input therefore!
if updating the database with new files creates a utf16 error while
reading in the database, open it in excel, save the file, and try
again. if this does not work, good luck...

'''
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from builtins import input
from builtins import zip
from builtins import str
from builtins import range
from past.utils import old_div

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
import xlrd

from .utils import *
from .data_plot import *

### find the path to the tree ###



class gdb(DataPlot, Utils):
    '''
    This class provides easy access to the presolar grain databse, as
    described in the header The databse is read in by default, however
    you can choose a private database and if you do so, if you want to
    use the private database exclusively or together with the whole
    database.

    If you use a private datafile, make sure it is in the same
    structure as the presolar grain database.  The most important
    thing is that you have a column named 'Notes' before the data
    start.  Everything on the right there is data, on the left
    description of the data.  You don't have to have all the data
    columns if you don't have data, but the header of the data columns
    needs to be exactly the same as the in the database.  Look at
    validation/grain_data xls files if you want to see an example of
    the formatting / or use it as a template.

    Parameters
    ----------
    fname : string
        filename to your private database, if not in the main tree
        structure with other databse, give full path.
    gdbdir : string
        In case you have not a full svn tree installed, choose the correct
        path to all the data here
    gdbload : boolean, optional
        True or False: Do you want to load the grain database or not?
        The default is True.
    iniabufile : string, optional
        Which initial abundances should be used to calculate delta
        values.  Here we assume complete SVN tree (need USEEPP).  Give
        absolute path otherwise. The default is 'frames/mppnp/USEEPP/iniab2.0E-02GN93.ppn'.

    '''

    def __init__(self, fname=None, gdbdir=None, gdbload=True,
                 iniabufile='frames/mppnp/USEEPP/iniab2.0E-02GN93.ppn'):
        print('Reading in... this takes a little bit')

        if iniabufile[0] != '/':
            iniabufile = get_svnpath() + iniabufile

        # grab data
        header_desc, header_data, desc, data = preprocessor(fname,gdbdir,gdbload)
        # make dictionary
        descdict = dict(list(zip(header_desc,list(range(len(header_desc))))))
        datadict = dict(list(zip(header_data,list(range(len(header_data))))))

        # style definer
        header_style, style = style_creator(desc,descdict)
        styledict = dict(list(zip(header_style,list(range(len(header_style))))))

        # make private instances w/ all the data
        self._header_desc = header_desc
        self._header_data = header_data
        self._header_style = header_style
        self._desc = desc
        self._data = data
        self._style = style
        self._descdict = descdict
        self._datadict = datadict
        self._styledict = styledict
        # make the working data
        self.header_desc = header_desc
        self.header_data = header_data
        self.header_style = header_style
        self.desc = desc
        self.data = data
        self.style = style
        self.descdict = descdict
        self.datadict = datadict
        self.styledict = styledict
        self.inut = iniabu(iniabufile)

    def __del__(self):
        print('Presolar grain database available at: http://presolar.wustl.edu/PGD/Presolar_Grain_Database.html')

    def reset_filter(self):
        '''
        Resets the filter and goes back to initialized value. This
        routine also resets the style if you have changed it.

        '''
        self.header_desc = self._header_desc
        self.header_data = self._header_data
        self.header_style = self._header_style
        self.desc = self._desc
        self.data = self._data
        self.style = self._style
        self.descdict = self._descdict
        self.datadict = self._datadict
        self.styledict = self._styledict

    def info(self, graintype=True, group=True, reference=False,
             phase=True):
        '''
        This routine gives you informations what kind of grains are
        currently available in your filtered version.  It gives you
        the type of grains available.  More to be implemented upon need.

        Parameters
        ----------
        graintype, group, references, phase : boolean
            What do you wanna print for information.  There can be a
            lot of references, hence references default is False.

        '''

        # create a list with all graintypes
        gtype_info = []
        group_info = []
        ref_info = []
        phase_info = []

        # how many grains in database
        print('There are ' + str(len(self.data)) + ' grains in your database.\n')

        # graintypes
        if graintype:
            for i in range(len(self.desc)):
                gtype_tmp = self.desc[i][self.descdict['Type']]
                wrtchk = True
                for j in range(len(gtype_info)):
                    if gtype_info[j] == gtype_tmp:
                        wrtchk = False
                        break
                if wrtchk:
                    gtype_info.append(gtype_tmp)
            print('Available graintypes are:')
            print('-------------------------')
            print(gtype_info)

        # groups
        if group:
            for i in range(len(self.desc)):
                group_tmp = self.desc[i][self.descdict['Group']]
                wrtchk = True
                for j in range(len(group_info)):
                    if group_info[j] == group_tmp:
                        wrtchk = False
                        break
                if wrtchk:
                    group_info.append(group_tmp)
            print('\nAvailable groups of grains (for silicates and oxides) are:')
            print('----------------------------------------------------------')
            print(group_info)

        # Phases
        if phase:
            for i in range(len(self.desc)):
                phase_tmp = self.desc[i][self.descdict['Phase']]
                wrtchk = True
                for j in range(len(phase_info)):
                    if phase_info[j] == phase_tmp:
                        wrtchk = False
                        break
                if wrtchk:
                    phase_info.append(phase_tmp)
            print('\nAvailable Phases of grains are:')
            print('----------------------------------------------------------')
            print(phase_info)


        # references
        if reference:
            for i in range(len(self.desc)):
                ref_tmp = self.desc[i][self.descdict['Reference']]
                wrtchk = True
                for j in range(len(ref_info)):
                    if ref_info[j] == ref_tmp:
                        wrtchk = False
                        break
                if wrtchk:
                    ref_info.append(ref_tmp)
            print('\nReferences for grains:')
            print('----------------------')
            print(ref_info)





    def filter_desc(self, graintype=None, group=None, reference=None,
                    size=None, phase=None):
        '''
        This routine is to filter for description elements.  You can
        check what is available in the description by running,

        >>> i.header_desc()

        where i is the instance you loaded.

        You can run the filter multiple times!  You can filter for the
        following types:

        Parameters
        ----------
        graintype : string or list
            Give graintypes as either 'M' for only mainstream or more
            than one ['M','Z'].
        group : integer or list
            Group of graintypes, important for oxides and silicates,
            since they are split into groups and not into types.
            Example 1, or give a list [1,3].
        reference : string or list
            Give the reference you want to filter for, try an i.info()
            to pick the right name!  You can select a single
            referennce as string or multiple references in as a list.
        size : string
            Filter for grain sizes, give '<5.0' or '>5.0' as a string
            for larger or smaller than a given grainsize in um.  Only
            data with known grainsizes are chosen.  Often grain sizes
            are given in a times b, where a and b are the minumum and
            maximum measurements from an image.  If you give a >5.0
            then grains with the smaller dimension >5um are taken into
            account.  If you want <5.0 then grains with the upper
            dimension <5um are taken into account.

        '''

        # filter for graintype
        if graintype != None:
            indexing = []
            # index file on which lines to pick
            if type(graintype) == str:
                graintype = [graintype]
            # filter
            for typ in graintype:
                for i in range(len(self.desc)):
                    if self.desc[i][self.descdict['Type']] == typ:
                        indexing.append(i)
            # filter:
            self._filter_desc(indexing)

        # filter for graintype
        if phase != None:
            indexing = []
            # index file on which lines to pick
            if type(phase) == str:
                phase = [phase]
            # filter
            for typ in phase:
                for i in range(len(self.desc)):
                    if self.desc[i][self.descdict['Phase']] == typ:
                        indexing.append(i)
            # filter:
            self._filter_desc(indexing)


        # filter for group (oxides and silicates)
        if group != None:
            indexing = []
            # index file on which lines to pick
            if type(group) != list:
                group = [group]
            # filter
            for grp in group:
                for i in range(len(self.desc)):
                    if self.desc[i][self.descdict['Group']] == str(int(grp)):
                        indexing.append(i)
            # filter:
            self._filter_desc(indexing)


        # filter for reference
        if reference != None:
            indexing = []
            # index file on which lines to pick
            if type(reference) != list:
                reference = [reference]
            # filter
            for ri in range(len(reference)):
                for i in range(len(self.desc)):
                    if self.desc[i][self.descdict['Reference']] == reference[ri]:
                        indexing.append(i)
            # filter:
            self._filter_desc(indexing)


        # filter for grainzise
        if size != None:
            indexing = []
            # index file on which lines to pick
            # filter
            operator = size[0:1]
            size = float(size[1:len(size)])
            for i in range(len(self.desc)):
                if  self.desc[i][self.descdict['Size (microns)']] != '':
                    try:
                        # print self.desc[i][self.descdict['Size (microns)']]
                        comperator1 = self.desc[i][self.descdict['Size (microns)']].split('x')[0]
                        comperator2 = self.desc[i][self.descdict['Size (microns)']].split('x')[1]
                        comperator = [float(comperator1),float(comperator2)]
                        if operator == '<':
                            comperator = np.min(comperator)
                        else:
                            comperator = np.max(comperator)
                    except IndexError or AttributeError:
                        try:
                            comperator = float(self.desc[i][self.descdict['Size (microns)']])
                        except ValueError:
                            continue

                    if operator == '>':
                        if comperator > size:
                            indexing.append(i)
                    elif operator == '<':
                        if comperator < size:
                            indexing.append(i)
                    else:
                        continue

            # filter:
            self._filter_desc(indexing)


    def _filter_desc(self, indexing):
        '''
        Private function to filter data, goes with filter_desc

        '''

        # now filter data
        if len(indexing) > 0:
            desc_tmp = np.zeros((len(indexing),len(self.header_desc)),dtype='|S1024')
            data_tmp = np.zeros((len(indexing),len(self.header_data)))
            style_tmp= np.zeros((len(indexing),len(self.header_style)),dtype='|S1024')
            for i in range(len(indexing)):
                for j in range(len(self.header_desc)):
                    desc_tmp[i][j] = self.desc[indexing[i]][j]
                for k in range(len(self.header_data)):
                    data_tmp[i][k] = self.data[indexing[i]][k]
                for l in range(len(self.header_style)):
                    style_tmp[i][l]= self.style[indexing[i]][l]
            self.desc = desc_tmp
            self.data = data_tmp
            self.style= style_tmp
        else:
            print('No filter selected or no data found!')



    def filter_single_grain(self):
        '''
        This subroutine is to filter out single grains. It is kind of
        useless if you have tons of data still in the list. To work on
        there, you have other filters (filter_desc and filter_data)
        available!  This filter gives an index to every grain, plots
        the most important information, and then asks you to pick a
        filter.  No input necessary, input is given during the routine

        '''

        my_index = 0
        my_grains = [['Index','Label','Type','Group','Meteorite','Mineralogy','C12/C13','d(Si29/Si30)','d(Si30/Si29)']]

        # add the data to this grain list
        for it in range(len(self.data)):
            my_grains.append([my_index,self.desc[it][self.descdict['Grain Label']], self.desc[it][self.descdict['Type']], self.desc[it][self.descdict['Group']], self.desc[it][self.descdict['Meteorite']], self.desc[it][self.descdict['Mineralogy']], self.data[it][self.datadict['12c/13c']], self.data[it][self.datadict['d(29si/28si)']], self.data[it][self.datadict['d(30si/28si)']]])
            my_index += 1

        for prt_line in my_grains:
            print(prt_line)

        # now write the selector for the index of the grains to select which one should be
        # available and which ones should be dumped
        usr_input = ''
        usr_input = input('Select the grains by index that you want to use. Please separate the indeces by a comma, e.g., 1 or 0,2,3,4\n')

        # process user index
        if usr_input == '':
            print('No data selected to filter.')
            return None
        elif len(usr_input) == 1:
            usr_index = [usr_input]
        else:
            usr_index = usr_input.split(',')
        for it in range(len(usr_index)):
            usr_index[it] = int(usr_index[it])

        # filter
        desc_tmp = np.zeros((len(usr_index),len(self.header_desc)),dtype='|S1024')
        data_tmp = np.zeros((len(usr_index),len(self.header_data)))
        style_tmp= np.zeros((len(usr_index),len(self.header_style)),dtype='|S1024')
        for i in range(len(usr_index)):
            for j in range(len(self.header_desc)):
                desc_tmp[i][j] = self.desc[usr_index[i]][j]
            for k in range(len(self.header_data)):
                data_tmp[i][k] = self.data[usr_index[i]][k]
            for l in range(len(self.header_style)):
                style_tmp[i][l]= self.style[usr_index[i]][l]
        self.desc = desc_tmp
        self.data = data_tmp
        self.style= style_tmp




    def filter_data(self, isos, limit, delta=True):
        '''
        This subroutine filters isotopic values according to the limit
        you give.  You can filter in ratio or in delta space.

        Parameters
        ----------
        isos : list
            isotopes you want to filter for, e.g., give as
            ['Si-28', 'Si-30'] for the 28/30 ratio.
        limit : string
            what do you want to filter for, e.g., ratio or delta > 100,
            then give '>100'.
        delta : boolean, optional
            do you wanna filter in delta space, then set to True,
            otherwise to False.  The default is True.

        '''

        # check availability
        dat_index, delta_b, ratio_b = self.check_availability(isos)
        if dat_index == -1:
            print('Isotopes selected are not available. Check i.datadict (where i is your instance) for availability of isotopes.')
            return None

        # select if larger or smaller and define limit
        if limit[0:1] == '>':
            comperator = 'gt'
        elif limit[0:1] == '<':
            comperator = 'st'
        else:
            print('Comperator not specified. Limit must be given as \'>5.\' for example.')
            return None

        try:
            limit = float(limit[1:len(limit)])
        except ValueError:
            print('Limit must be given as \'>5.\' for example.')
            return None

        # now calculate the actual limit to compare with, depending on if it delta or not or whatsoever
        if delta == delta_b:   # input and available same
            if ratio_b:   # one over
                if delta:
                    tmp = self.delta_to_ratio(isos,limit,oneover=True)
                    comp_lim = self.ratio_to_delta(isos,tmp)   # check
                else:
                    comp_lim = old_div(1.,limit)   # check

            else:   # all fine
                comp_lim = limit

        else:   # input and availability not the same
            if ratio_b:   # one over
                if delta:   # delta given, ratio one over wanted
                    comp_lim = self.delta_to_ratio(isos,limit,oneover=True)
                else:   # ratio given, delta one over wanted
                    comp_lim = self.ratio_to_delta(isos,limit,oneover=True)
            else:   # not one over
                if delta:   # delta given, ratio wanted
                    comp_lim = self.delta_to_ratio(isos,limit)
                else:
                    comp_lim = self.ratio_to_delta(isos,limit)

        # indexing vector
        indexing = []
        for i in range(len(self.data)):
            dat_val = self.data[i][dat_index]
            if comperator == 'st':
                if dat_val < comp_lim:
                    indexing.append(i)
            else:
                if dat_val > comp_lim:
                    indexing.append(i)

        # now filter data
        if len(indexing) > 0:
            desc_tmp = np.zeros((len(indexing),len(self.header_desc)),dtype='|S1024')
            data_tmp = np.zeros((len(indexing),len(self.header_data)))
            for i in range(len(indexing)):
                for j in range(len(self.header_desc)):
                    desc_tmp[i][j] = self.desc[indexing[i]][j]
                for k in range(len(self.header_data)):
                    data_tmp[i][k] = self.data[indexing[i]][k]
            self.desc = desc_tmp
            self.data = data_tmp
        else:
            print('No filter selected!')




    def filter_uncertainty(self, isos, limit, delta=True):
        '''
        This subroutine filters isotopic values according to the limit
        you give.  You can filter in ratio or in delta space.  This
        routine is based on the uncertainties, e.g., if you want to
        select only high quality data.

        Parameters
        ----------
        isos : list
            Isotopes you want to filter for, e.g., give as
            ['Si-28', 'Si-30'] for the 28/30 ratio.
        limit : string
            What do you want to filter for, e.g., ratio or delta > 100,
            then give '>100'.
        delta : boolean, optional
            Do you wanna filter in delta space, then set to True,
            otherwise to False.  The default is True.
        '''

        # check availability
        dat_index, delta_b, ratio_b = self.check_availability(isos)
        if dat_index == -1:
            print('Isotopes selected are not available. Check i.datadict (where i is your instance) for availability of isotopes.')
            return None

        # select if larger or smaller and define limit
        if limit[0:1] == '>':
            comperator = 'gt'
        elif limit[0:1] == '<':
            comperator = 'st'
        else:
            print('Comperator not specified. Limit must be given as \'>5.\' for example.')
            return None

        try:
            limit = float(limit[1:len(limit)])
        except ValueError:
            print('Limit must be given as \'>5.\' for example.')
            return None

        # now calculate the actual limit to compare with, depending on if it delta or not or whatsoever
        if delta == delta_b:   # input and available same
            if ratio_b:   # one over
                if delta:
                    tmp = self.delta_to_ratio(isos,limit,oneover=True)
                    comp_lim = self.ratio_to_delta(isos,tmp)   # check
                else:
                    comp_lim = old_div(1.,limit)   # check

            else:   # all fine
                comp_lim = limit

        else:   # input and availability not the same
            if ratio_b:   # one over
                if delta:   # delta given, ratio one over wanted
                    comp_lim = self.delta_to_ratio(isos,limit,oneover=True)
                else:   # ratio given, delta one over wanted
                    comp_lim = self.ratio_to_delta(isos,limit,oneover=True)
            else:   # not one over
                if delta:   # delta given, ratio wanted
                    comp_lim = self.delta_to_ratio(isos,limit)
                else:
                    comp_lim = self.ratio_to_delta(isos,limit)

        # indexing vector
        indexing = []
        for i in range(len(self.data)):
            dat_val = self.data[i][dat_index+1]
            if comperator == 'st':
                if dat_val < comp_lim:
                    indexing.append(i)
            else:
                if dat_val > comp_lim:
                    indexing.append(i)

        # now filter data
        if len(indexing) > 0:
            desc_tmp = np.zeros((len(indexing),len(self.header_desc)),dtype='|S1024')
            data_tmp = np.zeros((len(indexing),len(self.header_data)))
            for i in range(len(indexing)):
                for j in range(len(self.header_desc)):
                    desc_tmp[i][j] = self.desc[indexing[i]][j]
                for k in range(len(self.header_data)):
                    data_tmp[i][k] = self.data[indexing[i]][k]
            self.desc = desc_tmp
            self.data = data_tmp
        else:
            print('No filter selected!')


    def style_chg_label(self,type,symb=None,edc=None,fac=None,smbsz=None,edw=None,lab=None):
        '''
        This routine changes the plotting style that is set by default.
        The style is changed according the the label that you choose.
        Changing according to reference, use style_chg_ref() function!
        You can change it back to default by resetting the filter using
        g.reset_filter() routine, assuming that g is your instance. The
        format that is used here is:

        ['Symbol', 'Edge color', 'Face color', 'Symbol size', 'Edge width'
        ,'Label']

        You can see the current styles by running

        Attention: You have to give values to all variables that are
        compatible with the python mathplotlib. If not, it's your fault
        if nothing works.

        g.style

        Parameters
        ----------
        type : string
            Select the label of the grains you want to change.
        symb : string, optional
            Select new symbol. None for no change.
        edc : string, optional
            Select new edge color. None for no change.
        fac : string, optional
            Select new face color. None for no change.
        smbsz : string, optional
            Select new symbol size. None for no change.
        edw : string, optional
            Select new edge width. None for no change.
        lab : string, optional
            Select new label. None for no change. Watch out, if you
            want to do more specifications later, the type will
            have changed to the new label.
        '''

        # do stuff for selected type
        for i in range(len(self.style)):
            # check if type is correct, otherwise continue directly
            if self.style[i][self.styledict['Label']] == type:
                # change symbol:
                if symb != None:
                    self.style[i][self.styledict['Symbol']] = symb
                # change edge color
                if edc != None:
                    self.style[i][self.styledict['Edge color']] = edc
                # change face color
                if fac != None:
                    self.style[i][self.styledict['Face color']] = fac
                # change symbol size
                if smbsz != None:
                    self.style[i][self.styledict['Symbol size']] = smbsz
                # change edge width
                if edw != None:
                    self.style[i][self.styledict['Edge width']] = edw
                # change label
                if lab != None:
                    self.style[i][self.styledict['Label']] = lab

    def style_chg_ref(self,ref,symb=None,edc=None,fac=None,smbsz=None,edw=None,lab=None):
        '''
        This routine changes the plotting style that is set by default.
        The style is changed according the the reference of the paper
        as given in the grain database. For change according to type of
        grain, use the routine syle_chg_label().

        ['Symbol', 'Edge color', 'Face color', 'Symbol size', 'Edge width'
        ,'Label']

        You can see the current styles by running

        Attention: You have to give values to all variables that are
        compatible with the python mathplotlib. If not, it's your fault
        if nothing works.

        g.style

        Parameters
        ----------
        ref : string
            Select the reference of the grains you want to change.
        symb : string, optional
            Select new symbol. None for no change.
        edc : string, optional
            Select new edge color. None for no change.
        fac : string, optional
            Select new face color. None for no change.
        smbsz : string, optional
            Select new symbol size. None for no change.
        edw : string, optional
            Select new edge width. None for no change.
        lab : string, optional
            Select new label. None for no change.
        '''

        # do stuff for selected reference
        for i in range(len(self.style)):
            # check if reference is correct, otherwise continue directly
            if self.desc[i][self.descdict['Reference']] == ref:
                # change symbol:
                if symb != None:
                    self.style[i][self.styledict['Symbol']] = symb
                # change edge color
                if edc != None:
                    self.style[i][self.styledict['Edge color']] = edc
                # change face color
                if fac != None:
                    self.style[i][self.styledict['Face color']] = fac
                # change symbol size
                if smbsz != None:
                    self.style[i][self.styledict['Symbol size']] = smbsz
                # change edge width
                if edw != None:
                    self.style[i][self.styledict['Edge width']] = edw
                # change label
                if lab != None:
                    self.style[i][self.styledict['Label']] = lab



    ##### PLOTTING PREPARATOR #####

    def plot_ratio_return(self, isox, isoy, deltax=True, deltay=True):
        '''
        This routine returns data isotopic data to plot from the
        filtered list of data.

        Parameters
        ----------
        isox : list
            Isotopes for x axis in standard format ['Si-28','Si-30'].
        isoy : list
            Same as isox but for y axis.
        deltax : boolean, optional
            If true then x-axis values are in delta format.  The default
            is True.
        deltay : boolean, optional
            Same as for x-axis but for y-axis.  The default is True.

        Returns
        -------
        grpl_xdata
            grain plot x-axis data.
        grpl_xerr
            x-axis error bars.
        grpl_ydata
            grain plot y-axis data.
        grpl_yerr
            y-axis error bars.
        grpl_style
            style data for the different symbols.

        '''

        # check availability
        index_x, delta_b_x, ratio_b_x = self.check_availability(isox)
        index_y, delta_b_y, ratio_b_y = self.check_availability(isoy)
        if index_x == -1 or index_y == -1:
            print('Following input data are not available in the database. Revise your input.')
            if index_x == -1:
                print('x axis data not available')
            if index_y == -1:
                print('y axis data not available')
            return None

        # create x and y data as 1d vectors, also error bars
        xdata_vec = np.zeros((len(self.data)))
        ydata_vec = np.zeros((len(self.data)))
        xdata_err = np.zeros((len(self.data)))
        ydata_err = np.zeros((len(self.data)))
        for it in range(len(self.data)):
            xdata_vec[it] = self.data[it][index_x]
            ydata_vec[it] = self.data[it][index_y]
            xdata_err[it] = self.data[it][index_x+1]
            ydata_err[it] = self.data[it][index_y+1]

        # index data that are nan
        index_nan = []
        for it in range(len(xdata_vec)):
            if np.isnan(xdata_vec[it]) or np.isnan(ydata_vec[it]):
                index_nan.append(it)

        # make range of all incides
        index_filtered = list(range(len(xdata_vec)))
        for it in range(len(index_nan)):
            index_filtered.remove(index_nan[it])

        xdata_tmp = np.zeros((len(index_filtered)))
        ydata_tmp = np.zeros((len(index_filtered)))
        xerr_tmp  = np.zeros((len(index_filtered)))
        yerr_tmp  = np.zeros((len(index_filtered)))
        style_plt = np.zeros((len(index_filtered),len(self.header_style)),dtype='|S1024')

        for i in range(len(index_filtered)):
            xdata_tmp[i] = xdata_vec[index_filtered[i]]
            ydata_tmp[i] = ydata_vec[index_filtered[i]]
            xerr_tmp[i]  = xdata_err[index_filtered[i]]
            yerr_tmp[i]  = ydata_err[index_filtered[i]]
            for j in range(len(style_plt[i])):
                style_plt[i][j] = self.style[index_filtered[i]][j]
        xdata_vec = xdata_tmp
        ydata_vec = ydata_tmp
        xdata_err = xerr_tmp
        ydata_err = yerr_tmp

        # loop through error and set nans to 0
        for i in range(len(xdata_err)):
            if np.isnan(xdata_err[i]):
                xdata_err[i] = 0.
            if np.isnan(ydata_err[i]):
                ydata_err[i] = 0.

        # make start stop index for groups
        start_stop = []
        start = 0
        for it in range(len(xdata_vec)-1):
            if (style_plt[it] == style_plt[it+1]).all():
                continue
            else:
                stop = it + 1
                start_stop.append([start,stop])
                start = stop
        # last entry
        if start_stop == []:
            start_stop.append([0,len(xdata_vec)])
        else:
            start_stop.append([start_stop[len(start_stop)-1][1],len(xdata_vec)+1])

        # now append things to return variables
        grain_plt_xdata = []
        grain_plt_ydata = []
        grain_plt_xerr  = []
        grain_plt_yerr  = []
        grain_plt_style = []
        for i in range(len(start_stop)):
            grain_plt_xdata.append(xdata_vec[start_stop[i][0]:start_stop[i][1]])
            grain_plt_ydata.append(ydata_vec[start_stop[i][0]:start_stop[i][1]])
            grain_plt_xerr.append(xdata_err[start_stop[i][0]:start_stop[i][1]])
            grain_plt_yerr.append(ydata_err[start_stop[i][0]:start_stop[i][1]])
            grain_plt_style.append(style_plt[start_stop[i][0]])


        return [grain_plt_xdata,grain_plt_xerr,grain_plt_ydata,grain_plt_yerr,grain_plt_style]


    def plot_pattern_return(self, isos, delta=True):
        '''
        This routine returns data isotopic data to plot from the
        filtered list of data.

        Parameters
        ----------
        isos : list
            Isotopes for x axis in standard format
            [['Si-30','Si-28'],['Si-29','Si-30'],...]
        isoy : list
            Same as isox but for y axis.
        deltay: boolean, optional
            Same as for x-axis but for y-axis.  The default is True.

        Returns
        -------
        grpl_data
            grain plot x-axis data.
        grpl_err
            x-axis error bars.
        grpl_style
            style data for the different symbols.

        '''
        # check availability
        index = []
        delta_b = []
        ratio_b = []
        for i in range(len(isos)):
            tmpi,tmpd,tmpr = self.check_availability(isos[i])
            index.append(tmpi)
            delta_b.append(tmpd)
            ratio_b.append(tmpd)

        for i in range(len(index)):
            if index[i] == -1:
                print('Input not available for: ' + isos[i] + '. Revise!')
                return None

        # create x and y data as 1d vectors, also error bars
        data_vec = np.zeros((len(self.data),len(isos)))
        data_err = np.zeros((len(self.data),len(isos)))
        for it in range(len(self.data)):
            for jt in range(len(isos)):
                data_vec[it][jt] = self.data[it][index[jt]]
                data_err[it][jt] = self.data[it][index[jt]+1]

        # index data that are nan
        index_nan = []
        for it in range(len(data_vec)):
            for jt in range(len(data_vec[it])):
                if np.isnan(data_vec[it][jt]):
                    index_nan.append(it)

        # make range of all incides
        index_filtered = list(range(len(data_vec)))
        for it in range(len(index_nan)):
            index_filtered.remove(index_nan[it])

        data_tmp = np.zeros((len(index_filtered),len(isos)))
        err_tmp  = np.zeros((len(index_filtered),len(isos)))
        style_plt = np.zeros((len(index_filtered),len(self.header_style)),dtype='|S1024')

        for i in range(len(index_filtered)):
            data_tmp[i] = data_vec[index_filtered[i]]
            err_tmp[i]  = data_err[index_filtered[i]]
            for j in range(len(style_plt[i])):
                style_plt[i][j] = self.style[i][j]
        xdata_vec = xdata_tmp
        ydata_vec = ydata_tmp
        xdata_err = xerr_tmp
        ydata_err = yerr_tmp

        # loop through error and set nans to 0
        for i in range(len(xdata_err)):
            if np.isnan(xdata_err[i]):
                xdata_err[i] = 0.
            if np.isnan(ydata_err[i]):
                ydata_err[i] = 0.
        # FIXME here
        # make start stop index for groups
        start_stop = []
        start = 0
        for it in range(len(xdata_vec)-1):
            if (style_plt[it] == style_plt[it+1]).all():
                continue
            else:
                stop = it
                start_stop.append([start,stop])
                start = stop+1
        # last entry
        if start_stop == []:
            start_stop.append([0,len(xdata_vec)])
        else:
            start_stop.append([start_stop[len(start_stop)-1][1]+1,len(xdata_vec)])

        # now append things to return variables
        grain_plt_xdata = []
        grain_plt_ydata = []
        grain_plt_xerr  = []
        grain_plt_yerr  = []
        grain_plt_style = []
        for i in range(len(start_stop)):

            grain_plt_xdata.append(xdata_vec[start_stop[i][0]:start_stop[i][1]])
            grain_plt_ydata.append(ydata_vec[start_stop[i][0]:start_stop[i][1]])
            grain_plt_xerr.append(xdata_err[start_stop[i][0]:start_stop[i][1]])
            grain_plt_yerr.append(ydata_err[start_stop[i][0]:start_stop[i][1]])
            grain_plt_style.append(style_plt[start_stop[i][0]])


        return [grain_plt_xdata,grain_plt_xerr,grain_plt_ydata,grain_plt_yerr,grain_plt_style]






    ##### SMALL HELPER ROUTINES #####
    def check_availability(self, isos):
        '''
        This routine checks if the requested set of isotopes is
        available in the dataset.

        Parameters
        ----------
        isos : list
            set of isotopes in format ['Si-28','Si-30'].

        Returns
        -------
        list
           [index, delta_b, ratio_b].
           index: where is it.
           delta_b: is it a delta value or not?
           ratio_ib:  True if ratio is inverted, false if not

        '''
        # make names
        iso1name =  iso_name_converter(isos[0])
        iso2name =  iso_name_converter(isos[1])

        ratio = iso1name + '/' + iso2name
        ratio_inv = iso2name + '/' + iso1name
        delta = 'd(' + iso1name + '/' + iso2name + ')'
        delta_inv = 'd(' + iso2name + '/' + iso1name + ')'

        index = -1
        # search for data entry
        try:
            index = self.datadict[ratio]
            delta_b = False
            ratio_b = False
        except KeyError:
            try:
                index = self.datadict[ratio_inv]
                delta_b = False
                ratio_b = True
            except KeyError:
                try:
                    index = self.datadict[delta]
                    delta_b = True
                    ratio_b = False
                except KeyError:
                    try:
                        index = self.datadict[delta_inv]
                        delta_b = True
                        ratio_b = True
                    except KeyError:
                        index = -1
                        delta_b = None
                        ratio_b = None
        return index, delta_b, ratio_b



    def ratio_to_delta(self, isos_ss, ratio, oneover=False):
        '''
        Transforms an isotope ratio into a delta value

        Parameters
        ----------
        isos_ss: list or float
            list w/ isotopes, e.g., ['N-14','N-15'] OR the solar
            system ratio.
        ratio : float
            ratio of the isotopes to transform.
        oneover : boolean
            take the inverse of the ratio before transforming (never
            inverse of delta value!).  The default is False.

        Returns
        -------
        float
            delta value

        '''
        # define if isos_ss is the ratio or the isotopes
        if type(isos_ss) == float:
            ss_ratio = isos_ss
        elif type(isos_ss) == list:
            ss_ratio = self.inut.isoratio_init(isos_ss)
        else:
            print('Check input of isos_ss into ratio_to_delta routine')
            return None

        # check if one over is necessary or not
        if oneover:
            ratio = old_div(1,ratio)

        # calculate delta value
        delta = (old_div(ratio, ss_ratio) - 1.) * 1000.

        return delta

    def delta_to_ratio(self, isos_ss, delta, oneover=False):
        '''
        Transforms a delta value into an isotopic ratio

        Parameters
        ----------
        isos_ss : list or float
            list w/ isotopes, e.g., ['N-14','N-15'] OR the solar
            system ratio.
        delta : float
            delta value of the isotopes to transform.
        oneover : boolean
            take the inverse of the ratio before returning it (never
            of the delta value!).  The default is False.

        Returns
        -------
        float
            delta value

        '''
        # define if isos_ss is the ratio or the isotopes
        if type(isos_ss) == float:
            ss_ratio = isos_ss
        elif type(isos_ss) == list:
            ss_ratio = self.inut.isoratio_init(isos_ss)
        else:
            print('Check input of isos_ss into ratio_to_delta routine')
            return None

        # transform to ratio
        ratio = (old_div(delta, 1000.) + 1) * ss_ratio

        # one over necessary or not?
        if oneover:
            ratio = old_div(1,ratio)

        return ratio



##### SMALL SUBROUTINES THAT DO NOT NEED TO BE INSIDE THE CLASS, KIND OF GRAIN.PY SPECIFIC UTILS #####
def iso_name_converter(iso):
    '''
    Converts the name of the given isotope (input), e.g., 'N-14' to
    14N as used later to compare w/ grain database.

    '''
    sp = iso.split('-')
    output = sp[1] + sp[0]

    return output.lower()


def get_svnpath():
    '''
    This subroutine gives back the path of the whole svn tree
    installation, which is necessary for the script to run.

    '''
    svnpathtmp = __file__
    splitsvnpath = svnpathtmp.split('/')
    if len(splitsvnpath) == 1:
        svnpath = os.path.abspath('.') + '/../../'
    else:
        svnpath = ''
        for i in range(len(splitsvnpath)-3):
            svnpath += splitsvnpath[i] + '/'
    return svnpath



############################
##### BIG PREPROCESSOR #####
############################
# subroutine that reads in data and splits into nice numpy arrays
def preprocessor(fname,gdbdir,gdbload, wb_sic=None):
    if gdbdir == None:
        # path to svn
        svnpathtmp = get_svnpath()
        gdbdir = svnpathtmp + 'validation/grain_data/'   # grain data directory

    # Initialize private file if available
    if fname != None:
        wb_pri = xlrd.open_workbook(gdbdir + fname)
        sh_pri = wb_pri.sheet_by_index(0)
        print('Private file ' + fname + ' initialized.')


    # Initialize grain database
    if gdbload:
        # SiC
        wb_sic = xlrd.open_workbook(gdbdir + 'SiC-All.xls')
        sh_sic = wb_sic.sheet_by_index(0)
        # Graphites
        wb_gra = xlrd.open_workbook(gdbdir + 'graphite-All.xls')
        sh_gra = wb_gra.sheet_by_index(0)
        # Oxides and Silicates
        wb_oxi = xlrd.open_workbook(gdbdir + 'oxide-silicate-all.xls')
        sh_oxi = wb_oxi.sheet_by_index(0)
        # Misc grains
        wb_mis = xlrd.open_workbook(gdbdir + 'miscellaneous-SiN.xls')
        sh_mis = wb_mis.sheet_by_index(0)

    # now bring all files together into one database (if private is not the only file specified)
    header_data = list()   # header for data
    header_desc = list()   # header for description

    if gdbload:
        # SiC - first file
        head = sh_sic.row_values(0)
        # now split up to header_data and header_desc
        headswtch=True   # switch from description to data
        for head_i in head:
            if headswtch:   # write description header
                writeswtch = True
                for i in range(len(header_desc)):
                    if header_desc[i] == head_i:
                        writeswtch = False
                if writeswtch:
                    if head_i.replace(' ','') != '':
                        header_desc.append(head_i)
                if len(head_i) >= 5:
                    if head_i[0:5] == 'Notes':
                        headswtch = False
            else:   # write data header
                writeswtch = True
                for i in range(len(header_data)):
                    if header_data[i] == head_i:
                        writeswtch = False
                if writeswtch:
                    if head_i.replace(' ','') != '':
                        header_data.append(head_i)
        # Graphites
        head = sh_gra.row_values(0)
        headswtch=True   # switch from description to data
        for head_i in head:
            if headswtch:   # write description header
                writeswtch = True
                for i in range(len(header_desc)):
                    if header_desc[i] == head_i:
                        writeswtch = False
                if writeswtch:
                    if head_i.replace(' ','') != '':
                        header_desc.append(head_i)
                if len(head_i) >= 5:
                    if head_i[0:5] == 'Notes':
                        headswtch = False
            else:   # write data header
                writeswtch = True
                for i in range(len(header_data)):
                    if header_data[i] == head_i:
                        writeswtch = False
                if writeswtch:
                    if head_i.replace(' ','') != '':
                        header_data.append(head_i)
        # Silicates and Oxides
        head = sh_oxi.row_values(0)
        headswtch=True   # switch from description to data
        for head_i in head:
            if headswtch:   # write description header
                writeswtch = True
                for i in range(len(header_desc)):
                    if header_desc[i] == head_i:
                        writeswtch = False
                if writeswtch:
                    if head_i.replace(' ','') != '':
                        header_desc.append(head_i)
                if len(head_i) >= 5:
                    if head_i[0:5] == 'Notes':
                        headswtch = False
            else:   # write data header
                writeswtch = True
                for i in range(len(header_data)):
                    if header_data[i] == head_i:
                        writeswtch = False
                if writeswtch:
                    if head_i.replace(' ','') != '':
                        header_data.append(head_i)
        # Misc
        head = sh_mis.row_values(0)
        headswtch=True   # switch from description to data
        for head_i in head:
            if headswtch:   # write description header
                writeswtch = True
                for i in range(len(header_desc)):
                    if header_desc[i] == head_i:
                        writeswtch = False
                if writeswtch:
                    if head_i.replace(' ','') != '':
                        header_desc.append(head_i)
                if len(head_i) >= 5:
                    if head_i[0:5] == 'Notes':
                        headswtch = False
            else:   # write data header
                writeswtch = True
                for i in range(len(header_data)):
                    if header_data[i] == head_i:
                        writeswtch = False
                if writeswtch:
                    if head_i.replace(' ','') != '':
                        header_data.append(head_i)

    # Private file
    if fname != None:
        head = sh_pri.row_values(0)
        headswtch=True   # switch from description to data
        for head_i in head:
            if headswtch:   # write description header
                writeswtch = True
                for i in range(len(header_desc)):
                    if header_desc[i] == head_i:
                        writeswtch = False
                if writeswtch:
                    if head_i.replace(' ','') != '':
                        header_desc.append(head_i)
                if len(head_i) >= 5:
                    if head_i[0:5] == 'Notes':
                        headswtch = False
            else:   # write data header
                writeswtch = True
                for i in range(len(header_data)):
                    if header_data[i] == head_i:
                        writeswtch = False
                if writeswtch:
                    if head_i.replace(' ','') != '':
                        header_data.append(head_i)

    # Raise error if nothing is specified
    if gdbload == False and fname == None:
        print('Nothing to load is specified!')
        return [],[],[],[]


    if gdbload:
        # Prepare the data -> description and data, fill it into appropriate forms
        # total amount of data entries
        sic_len = len(sh_sic.col_values(0))
        gra_len = len(sh_gra.col_values(0))
        oxi_len = len(sh_oxi.col_values(0))
        mis_len = len(sh_mis.col_values(0))
        totdata = sic_len + gra_len + oxi_len + mis_len - 4
        if fname != None:
            pri_len = len(sh_pri.col_values(0))
            totdata += pri_len - 1
        totdata = int(totdata)
    else:
        pri_len = len(sh_pri.col_values(0))
        totdata = pri_len - 1

    # initialize the description list
    descr = np.zeros((totdata,len(header_desc)+1),dtype='|S1024')   # string with 1024 characters


    # initialize data array
    data = np.zeros((totdata,len(header_data)))
    for i in range(len(data)):
        for j in range(len(data[i])):
            data[i][j] = nan

    jadder = 0
    if gdbload:
        # SiC
        sic_hdict = dict(list(zip(sh_sic.row_values(0),np.arange(len(sh_sic.row_values(0))))))
        # description data
        for i in range(len(header_desc)):
            try:
                dat_tmp = sh_sic.col_values(sic_hdict[header_desc[i]])
            except KeyError:
                dat_tmp = False
            if dat_tmp != False:
                for j in range(1,len(dat_tmp)):
                    if type(dat_tmp[j]) == float:
                        dat_tmp[j] = str(dat_tmp[j])
                    if dat_tmp[j].replace(' ','') != '':
                        try:
                            descr[j+jadder-1][i] = dat_tmp[j]
                        except UnicodeEncodeError:
                            descr[j+jadder-1][i] == ''      # actual data
                    descr[j+jadder-1][len(header_desc)] = 'SiC'

        for i in range(len(header_data)):
            try:
                dat_tmp = sh_sic.col_values(sic_hdict[header_data[i]])
            except KeyError:
                dat_tmp = False
            if dat_tmp!= False:
                for j in range(1,len(dat_tmp)):
                    if type(dat_tmp[j]) != float:
                        dat_tmp_append = nan
                    else:
                        dat_tmp_append = float(dat_tmp[j])
                    data[j-1][i] = dat_tmp_append

        jadder += sic_len - 1

        # Graphites
        gra_hdict = dict(list(zip(sh_gra.row_values(0),np.arange(len(sh_gra.row_values(0))))))
        # description data
        for i in range(len(header_desc)):
            try:
                dat_tmp = sh_gra.col_values(gra_hdict[header_desc[i]])
            except KeyError:
                dat_tmp = False
            if dat_tmp != False:
                for j in range(1,len(dat_tmp)):
                    if type(dat_tmp[j]) == float:
                        dat_tmp[j] = str(dat_tmp[j])
                    if dat_tmp[j].replace(' ','') != '':
                        try:
                            descr[j+jadder-1][i] = dat_tmp[j]
                        except UnicodeEncodeError:
                            descr[j+jadder-1][i] == ''
                    descr[j+jadder-1][len(header_desc)] = 'Graphites'
        # actual data
        for i in range(len(header_data)):
            try:
                dat_tmp = sh_gra.col_values(gra_hdict[header_data[i]])
            except KeyError:
                dat_tmp = False
            if dat_tmp!= False:
                for j in range(1,len(dat_tmp)):
                    if type(dat_tmp[j]) != float:
                        dat_tmp_append = nan
                    else:
                        dat_tmp_append = float(dat_tmp[j])
                    data[j+jadder-1][i] = dat_tmp_append

        jadder += gra_len - 1

        # Oxides
        oxi_hdict = dict(list(zip(sh_oxi.row_values(0),np.arange(len(sh_oxi.row_values(0))))))
        # description data
        for i in range(len(header_desc)):
            try:
                dat_tmp = sh_oxi.col_values(oxi_hdict[header_desc[i]])
            except KeyError:
                dat_tmp = False
            if dat_tmp != False:
                for j in range(1,len(dat_tmp)):
                    if type(dat_tmp[j]) == float:
                        dat_tmp[j] = str(dat_tmp[j])
                    if dat_tmp[j].replace(' ','') != '':
                        try:
                            descr[j+jadder-1][i] = dat_tmp[j]
                        except UnicodeEncodeError:
                            descr[j+jadder-1][i] == ''
                    descr[j+jadder-1][len(header_desc)] = 'Oxides, Silicates'
        # actual data
        for i in range(len(header_data)):
            try:
                dat_tmp = sh_oxi.col_values(oxi_hdict[header_data[i]])
            except KeyError:
                dat_tmp = False
            if dat_tmp!= False:
                for j in range(1,len(dat_tmp)):
                    if type(dat_tmp[j]) != float:
                        dat_tmp_append = nan
                    else:
                        dat_tmp_append = float(dat_tmp[j])
                    data[j+jadder-1][i] = dat_tmp_append

        jadder += oxi_len - 1

        # Misc
        mis_hdict = dict(list(zip(sh_mis.row_values(0),np.arange(len(sh_mis.row_values(0))))))
        # description data
        for i in range(len(header_desc)):
            try:
                dat_tmp = sh_mis.col_values(mis_hdict[header_desc[i]])
            except KeyError:
                dat_tmp = False
            if dat_tmp != False:
                for j in range(1,len(dat_tmp)):
                    if type(dat_tmp[j]) == float:
                        dat_tmp[j] = str(dat_tmp[j])
                    if dat_tmp[j].replace(' ','') != '':
                        try:
                            descr[j+jadder-1][i] = dat_tmp[j]
                        except UnicodeEncodeError:
                            descr[j+jadder-1][i] == ''
                    descr[j+jadder-1][len(header_desc)] = 'Misc'
        # actual data
        for i in range(len(header_data)):
            try:
                dat_tmp = sh_mis.col_values(mis_hdict[header_data[i]])
            except KeyError:
                dat_tmp = False
            if dat_tmp!= False:
                for j in range(1,len(dat_tmp)):
                    if type(dat_tmp[j]) != float:
                        dat_tmp_append = nan
                    else:
                        dat_tmp_append = float(dat_tmp[j])
                    data[j+jadder-1][i] = dat_tmp_append

        jadder += mis_len - 1




    # Private file
    if fname != None:
        pri_hdict = dict(list(zip(sh_pri.row_values(0),np.arange(len(sh_pri.row_values(0))))))
        # description data
        for i in range(len(header_desc)):
            try:
                dat_tmp = sh_pri.col_values(pri_hdict[header_desc[i]])
            except KeyError:
                dat_tmp = False
            if dat_tmp != False:
                for j in range(1,len(dat_tmp)):
                    if type(dat_tmp[j]) == float:
                        dat_tmp[j] = str(dat_tmp[j])
                    if dat_tmp[j].replace(' ','') != '':
                        try:
                            descr[j+jadder-1][i] = dat_tmp[j]
                        except UnicodeEncodeError:
                            descr[j+jadder-1][i] == ''      # actual data
                    descr[j+jadder-1][len(header_desc)] = 'Private'
        for i in range(len(header_data)):
            try:
                dat_tmp = sh_pri.col_values(pri_hdict[header_data[i]])
            except KeyError:
                dat_tmp = False
            if dat_tmp!= False:
                for j in range(1,len(dat_tmp)):
                    if type(dat_tmp[j]) != float:
                        dat_tmp_append = nan
                    else:
                        dat_tmp_append = float(dat_tmp[j])
                    data[j+jadder-1][i] = dat_tmp_append




    # make the data header lower case
    for i in range(len(header_data)):
        header_data[i] = header_data[i].lower()

    header_desc.append('Database')

    return header_desc, header_data, descr, data


### style creator ###
def style_creator(desc,descdict):
    # make style definitions for plotting and everything
    header_style = ['Symbol', 'Edge color', 'Face color', 'Symbol size', 'Edge width','Label']
    style = np.zeros((len(desc),len(header_style)),dtype='|S1024')
    # fill the styles in
    for i in range(len(style)):
        style[i][4] = '1'   # edge width
        # SiC grains
        if desc[i][descdict['Database']] == 'SiC':
            style[i][0] = 'o'   # symbol
            style[i][3] = '8'   # symbol size
            if desc[i][descdict['Type']] == 'M':
                style[i][1] = '0.4'
                style[i][5] = 'SiC M'
            elif desc[i][descdict['Type']] == 'X':
                style[i][1] = 'b'
                style[i][5] = 'SiC X'
            elif desc[i][descdict['Type']] == 'Y':
                style[i][1] = 'g'
                style[i][5] = 'SiC Y'
            elif desc[i][descdict['Type']] == 'Z':
                style[i][1] = 'r'
                style[i][5] = 'SiC Z'
            elif desc[i][descdict['Type']] == 'AB':
                style[i][1] = 'c'
                style[i][5] = 'SiC AB'
            elif desc[i][descdict['Type']] == 'C' or desc[i][descdict['Type']] == 'U/C':
                style[i][1] = 'y'
                style[i][5] = 'SiC C'
            elif desc[i][descdict['Type']] == 'N':
                style[i][1] = 'm'
                style[i][5] = 'SiC nova'
            else:
                style[i][1] = '0.7'
                style[i][5] = 'SiC unclassified'
            style[i][2] = style[i][1]

        elif desc[i][descdict['Database']] == 'Graphites':
            style[i][0] = 's'   # symbol
            style[i][3] = '8'   # symbol size
            if desc[i][descdict['Type']] == 'HD':
                style[i][1] = '0.4'
                style[i][5] = 'Graphite HD'
            elif desc[i][descdict['Type']] == 'LD':
                style[i][1] = 'b'
                style[i][5] = 'Graphite LD'
            else:
                style[i][1] = '0.7'
                style[i][5] = 'Graphite'
            style[i][2] = style[i][1]

        elif desc[i][descdict['Database']] == 'Oxides, Silicates':
            style[i][0] = '^'   # symbol
            style[i][3] = '8'   # symbol size
            if desc[i][descdict['Group']] == '1':
                style[i][1] = '0.4'
                style[i][5] = 'Oxide / Silicate Group 1'
            elif desc[i][descdict['Group']] == '2':
                style[i][1] = 'b'
                style[i][5] = 'Oxide / Silicate Group 2'
            elif desc[i][descdict['Group']] == '3':
                style[i][1] = 'g'
                style[i][5] = 'Oxide / Silicate Group 3'
            elif desc[i][descdict['Group']] == '4':
                style[i][1] = 'r'
                style[i][5] = 'Oxide / Silicate Group 4'
            else:
                style[i][1] = '0.7'
                style[i][5] = 'Oxide / Silicate'
            style[i][2] = style[i][1]

        elif desc[i][descdict['Database']] == 'Misc':
            style[i][0] = 'v'
            style[i][3] = '8'
            style[i][1] = '0.4'
            style[i][2] = style[i][1]
            style[i][5] = 'Misc'

        elif desc[i][descdict['Database']] == 'Private':
            style[i][0] = '>'
            style[i][1] = '0.4'
            style[i][3] = '8'
            style[i][2] = style[i][1]
            style[i][5] = 'Private'

        else:
            style[i][0] = '+'
            style[i][1] = '0.4'
            style[i][3] = '8'
            style[i][2] = style[i][1]
            style[i][5] = 'Unknown'

    return header_style, style
