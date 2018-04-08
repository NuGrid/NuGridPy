
#
# NuGridpy - Tools for accessing and visualising NuGrid data.
#
# Copyright 2007 - 2014 by the NuGrid Team.
# All rights reserved. See LICENSE.
#

"""
data_plot.py

SuperClass module for the YProfile and Table Classes.  It contains
numerous plot function for the YProfile and Table Classes.

If one in the future wants their class to inherit this superclasses
methods this is what is required:

A. Place ``from data_table import *`` at the top of the module.
B. If the class is defined like 'class MyClass:', change that to
   'class MyClass(DataTable):'
C. To properly use DataTable's methods properly one will need these
   methods: a get(atri) that returns a numpy array of Data, or a list
   of numpy arrays of data.  The arguments of this function would need
   to be atri which is the name of the data one is looking for.

"""
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from builtins import zip
from builtins import str
from builtins import input
from builtins import range
from past.builtins import basestring
from builtins import object
from past.utils import old_div
from numpy import *
from math import *
import matplotlib.pylab as pyl
import matplotlib.pyplot as pl
#from matplotlib.mpl import colors,cm # depreciated in mpl ver 1.3
                                      # use line below instead
from matplotlib import colors,cm
import matplotlib
from matplotlib.patches import Rectangle, Arrow
from matplotlib.collections import PatchCollection
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea
from matplotlib.lines import Line2D
from matplotlib.ticker import *
from collections import OrderedDict
import numpy as np
import os
import os.path
import threading
import time
import sys

from . import astronomy as ast


def _padding_model_number(number, max_num):
    '''
    This method returns a zero-front padded string

    It makes out of str(45) -> '0045' if 999 < max_num < 10000. This is
    meant to work for reasonable integers (maybe less than 10^6).

    Parameters
    ----------
    number : integer
        number that the string should represent.
    max_num : integer
        max number of cycle list, implies how many 0s have be padded

    '''

    cnum = str(number)
    clen = len(cnum)

    cmax = int(log10(max_num)) + 1

    return (cmax - clen)*'0' + cnum


class DataPlot(object):

    _classTest_data = {
        'ppm.yprofile': 'YProfile',
        'ascii_table.ascii_table': 'AsciiTable',
        'nugridse.se': 'se',
        'mesa.mesa_profile': 'mesa_profile',
        'mesa.star_log': 'mesa.star_log',
        'mesa.history_data': 'mesa.star_log',
        'ppn.xtime': 'xtime',
        'ppn.abu_vector': 'PPN',
        'starobs.plot': 'starobs',
        'grain.gdb': 'grain',
        }

    def _classTest(self):
        '''
        Determines what the type of class instance the subclass is, so
        we can dynamically determine the behaviour of methods.

        The data this method uses (_classTest_data) NEEDS to be
        modified if any names of files or classes are changed.

        TODO - The entire use of this class needs to be refactored to use
        derived classes instead.
        '''
        c = '.'.join(str(self.__class__)[:-2].rsplit('.', 2)[-2:])
        return self._classTest_data.get(c, '')

    def _which(self, program):
        '''
        Mimics which in the unix shell.

        '''
        def is_exe(fpath):
            return os.path.exists(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file

        return None

    def _logarithm(self, tmpX, tmpY, logX, logY, base):
        logXER=False
        logYER=False
        for i in range(len(tmpX)):
            if tmpX[i]<=0. and logX:
                #print 'We can not log a number less than or equal to zero'
                #print 'Attempting to remove incompatible values from X'
                logXER=True
            if tmpY[i]<=0. and logY:
                #print 'We can not log a number less than or equal to zero'
                #print 'Attempting to remove incompatible values from Y'
                logYER=True
        tmX=[]
        tmY=[]

        if logXER:
            for i in range(len(tmpX)):
                if tmpX[i]>0.:
                    tmX.append( tmpX[i])
                    tmY.append(tmpY[i])
            tmpX=tmX
            tmpY=tmY
        elif logYER:
            for i in range(len(tmpY)):
                if tmpY[i]>0.:
                    tmX.append( tmpX[i])
                    tmY.append(tmpY[i])
            tmpX=tmX
            tmpY=tmY

        tmX=tmpX
        tmY=tmpY

        if logX:
            tmX=tmpX
            try:
                for i in range(len(tmpX)):
                    tmX[i]=log(tmpX[i],base)
            except ValueError:
                #print 'We can not log a number less than or equal to zero'
                #print 'Attempting to remove incompatible values from X'
                logXER=True
        if logY:
            tmY=tmpY
            try:
                for i in range(len(tmpY)):
                    tmY[i]=log(tmpY[i],base)
            except ValueError:
                #print 'We can not log a number less than or equal to zero'
                #print 'Attempting to remove incompatible values from Y'
                logYER=True

        if logX:
            tmpX=tmX
        if logY:
            tmpY=tmY

        return tmpX,tmpY

    def _sparse(self, x, y, sparse):
        """
        Method that removes every non sparse th element.

        For example:
        if this argument was 5, This method would plot the 0th, 5th,
        10th ... elements.

        Parameters
        ----------
        x : list
            list of x values, of length j.
        y : list
            list of y values, of length j.
        sparse : integer
            Argument that skips every so many data points.

        """
        tmpX=[]
        tmpY=[]

        for i in range(len(x)):
            if sparse == 1:
                return x,y
            if (i%sparse)==0:
                tmpX.append(x[i])
                tmpY.append(y[i])
        return tmpX, tmpY

    def plotMulti(self, atrix, atriy, cyclist, title, path='/',
                  legend=None, labelx=None, labely=None, logx=False,
                  logy=False, base=10, sparse=1, pdf=False,
                  limits=None):
        '''
        Method for plotting multiple plots and saving it to multiple
        pngs or PDFs.

        Parameters
        ----------
        atrix : string
            The name of the attribute you want on the x axis.
        atriy : string
            The name of the attribute you want on the Y axis.
        cyclist : list
            List of cycles that you would like plotted.
        title : string
            The title of the graph and the name of the file.
        path : string, optional
            The file path. The default is '/'
        Legend : list or intager, optional
            A list of legends for each of your cycles, or one legend for
            all of the cycles. The default is None.
        labelx : string, optional
            The label on the X axis. The default is None.
        labely : string, optional
            The label on the Y axis. The default is None.
        logx : boolean, optional
            A boolean of whether the user wants the x axis
            logarithmically. The default is False.
        logy : boolean, optional
            A boolean of whether the user wants the Y axis
            logarithmically. The default is False.
        base : integer, optional
            The base of the logarithm. The default is 10.
        sparse : integer, optional
            Argument that skips every so many data points.  For example
            if this argument was 5, This method would plot the 0th,
            5th, 10th ... elements. The default is 1.
        pdf : boolean, optional
            A boolean of if the image should be saved to a pdf file.
            xMin, xMax, yMin, YMax:  plot coordinates.  The default is
            False.
        limits : list, optional
            The length four list of the x and y limits.  The order of
            the list is xmin, xmax, ymin, ymax. The default is None.

        '''
        if str(legend.__class__)!="<type 'list'>":# Determines the legend is a list
            legendList=False
        else:
            legendList=True

        if legendList and len(cyclist) !=len(legend): #if it is a list, make sure there is an entry for each cycle
            print('Please input a proper legend, with correct length, aborting plot')
            return None
        for i in range(len(cyclist)):
            if legendList:
                self.plot(atrix,atriy,cyclist[i],'ndump',legend[i],labelx,labely,base=base,sparse=sparse, \
                                  logx=logx,logy=logy,show=False,limits=limits)
            else:
                self.plot(atrix,atriy,cyclist[i],'ndump',legend,labelx,labely,base=base,sparse=sparse, \
                                  logx=logx,logy=logy,show=False,limits=limits)

            pl.title(title)
            if not pdf:
                currentDir = os.getcwd()
                os.chdir(path)
                pl.savefig(title+str(cyclist[i])+'.png', dpi=400)
                os.chdir(currentDir)
            else:
                currentDir = os.getcwd()
                os.chdir(path)
                pl.savefig(title+str(cyclist[i])+'.pdf', dpi=400)
                os.chdir(currentDir)
            pl.clf()
        return None

    def plot(self, atrix, atriy, fname=None, numtype='ndump',
             legend=None, labelx=None, labely=None, indexx=None,
             indexy=None, title=None, shape='.', logx=False,
             logy=False, path='/', base=10, sparse=1, show=True, pdf=False,limits=None,
             markevery=None, linewidth=1):
        """
        Simple function that plots atriy as a function of atrix

        This method will automatically find and plot the requested data.

        Parameters
        ----------
        atrix : string
            The name of the attribute you want on the x axis.
        atriy : string
            The name of the attribute you want on the Y axis.
        fname : optional
            Be the filename, Ndump or time, or cycle, If fname is a
            list, this method will then save a png for each cycle in the
            list.  Warning, this must be a list of cycles and not a
            list of filenames. The default is None.
        numtype : string, optional
            designates how this function acts and how it interprets
            fname. if numtype is 'file', this function will get the
            desird attribute from that file.  if numtype is 'NDump'
            function will look at the cycle with that nDump.  if numtype
            is 't' or 'time' function will find the _cycle with the
            closest time stamp. The default is 'ndump'.
        legend : list or intager, optional
            A list of legends for each of your cycles, or one legend for
            all of the cycles. The default is None.
        labelx : string, optional
            The label on the X axis. The default is None.
        labely : string, optional
            The label on the Y axis. The default is None.
        indexx : optional
            Depreciated: If the get method returns a list of lists,
            indexx would be the list at the index indexx in the list.
            The default is None.
        indexy : optional
            Depreciated: If the get method returns a list of lists,
            indexy would be the list at the index indexx in the list.
            The default is None.
        title : string, optional
            The Title of the Graph. The default is None.
        shape : string, optional
            What shape and colour the user would like their plot in.
            Please see
            http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
            for all possible choices. The default is '.'.
        logx : boolean, optional
            A boolean of weather the user wants the x axi
            logarithmically. The default is False.
        logy : boolean, optional
            A boolean of weather the user wants the Y axis
            logarithmically. The default is False.
        path : string, optional
            Usef for PlotMulti, give the path where to save the Figures
        base : integer, optional
            The base of the logarithm. The Default is 10.
        sparse : integer, optional
            Argument that skips every so many data points. For example
            if this argument was 5, This method would plot the 0th, 5th,
            10th ... elements. The default is 1.
        show : boolean, optional
            A boolean of if the plot should be displayed useful with the
            multiPlot method. The default is True.
        pdf : boolean, optional
            PDF for PlotMulti? Default: False
        limits : list, optional
            The length four list of the x and y limits. The order of the
            list is xmin, xmax, ymin, ymax. The defautl is .
        markevery : integer or tupler, optional
            Set the markevery property to subsample the plot when
            using markers.  markevery can be None, very point will be
            plotted. It can be an integer N, Every N-th marker will be
            plotted starting with marker 0. It can be a tuple,
            markevery=(start, N) will start at point start and plot
            every N-th marker. The default is None.
        linewidth : integer, optional
            Set linewidth. The default is 1.

        Notes
        -----
        WARNING: Unstable if get returns a list with only one element (x=[0]).

        parameters: indexx and indexy have been deprecated.
        """
        t1=time.time()
        #Setting the axis labels

        if labelx== None :
            labelx=atrix
        if labely== None :
            labely=atriy

        if title!=None:
            title=title
        else:
            title=labely+' vs '+labelx

        if str(fname.__class__)=="<type 'list'>":
            self.plotMulti(atrix,atriy,fname,title,path,legend,labelx,labely,logx, logy, 10,1,pdf,limits)
            return
        tmpX=[]
        tmpY=[]
        singleX=False
        singleY=False
        #Getting data
        plotType=self._classTest()
        if plotType=='YProfile':
            if fname==None:
                fname=self.cycles[-1]

            listY=self.get(atriy,fname, numtype,resolution='a')
            listX=self.get(atrix,fname, numtype,resolution='a')
        elif plotType=='se':
            if fname==None:
                listY=self.get( atriy,sparse=sparse)
                listX=self.get(atrix,sparse=sparse)
            else:
                listY=self.get(fname, atriy,sparse=sparse)
                listX=self.get(fname, atrix,sparse=sparse)

            t2= time.time()
            print(t2 -t1)
        elif plotType=='PPN' :
            if fname==None and atrix not in self.cattrs and atriy not in self.cattrs:
                fname=len(self.files)-1
            if numtype=='ndump':
                numtype='cycNum'
            listY=self.get(atriy,fname,numtype)
            listX=self.get(atrix,fname,numtype)
        elif plotType=='xtime' or plotType=='mesa_profile' or plotType=='AsciiTable' or plotType=='mesa.star_log' or plotType=='starobs':
            listY=self.get(atriy)
            listX=self.get(atrix)
        else:
            listY=self.get(atriy)
            listX=self.get(atrix)
        tmpX=[]
        tmpY=[]
        if isinstance(listX[0], basestring) or isinstance(listY[0], basestring):
            for i in range(len(listX)):
                if '*****' == listX[i] or '*****' == listY[i]:
                    print('There seems to be a string of * in the lists')
                    print('Cutting out elements in both the lists that have an index equal to or greater than the index of the location of the string of *')
                    break
                tmpX.append(float(listX[i]))
                tmpY.append(float(listY[i]))

            listX=tmpX
            listY=tmpY




        #Determining if listX is a list or a list of lists
        try:
            j=listX[0][0]
        except:
            singleX = True

        if len(listX) == 1:  # If it is a list of lists with one element.
            tmpX=listX[0]
        elif singleX == True:# If it is a plain list of values.
            tmpX=listX
        elif indexx==None and len(listX)>1: # If it is a list of lists of values.
                                            # take the largest
            tmpX=listX[0]
            for i in range(len(listX)):
                if len(tmpX)<len(listX[i]):
                    tmpX=listX[i]
        elif indexx<len(listX): # If an index is specified, use that index
            tmpX=listX[indexx]
        else:
            print('Sorry that indexx does not exist, returning None')
            return None

        #Determining if listY is a list or a list of lists
        try:
            j=listY[0][0]
        except:
            singleY = True

        if len(listY) == 1: # If it is a list of lists with one element.
            #print 'hello'
            tmpY=listY[0]
        elif singleY == True: # If it is a plain list of values.
            #print 'world'
            tmpY=listY
        elif indexy==None and len(listY)>1:# If it is a list of lists of values.
                                            # take the largest
            #print 'fourth'
            tmpY=listY[0]
            for i in range(len(listY)):
                if len(tmpY)<len(listY[i]):
                    tmpY=listY[i]
        elif indexy<len(listY): # If an index is specified, use that index
            #print 'sixth'
            tmpY=listY[indexy]
        else:
            print('Sorry that indexy does not exist, returning None')
            return None
        '''
        elif indexy==None and len(listY)==1:
                #print 'fifth'
                tmpY=listY
        '''




        #Here, if we end up with different sized lists to plot, it
        #searches for a list that is of an equal length
        if len(tmpY)!=len(tmpX):
            found=False
            print("It seems like the lists are not of equal length")
            print("Now attempting to find a compatible list for ListX")
            for i in range(len(listY)):
                if not singleY and len(tmpX)==len(listY[i]):
                    tmpY=listY[i]
                    found=True

            if not found:
                print("Now attempting to find a compatible list for ListY")
                for i in range(len(listX)):

                    if not singleX and len(tmpY)==len(listX[i]):
                        tmpX=listX[i]
                        found=True

            if found:
                print("Suitable list found")
            else:

                print("There is no suitalble list, returning None")
                return None
        if len(tmpY)!=len(tmpX) and single == True:
            print('It seems that the selected lists are of different\nsize, now returning none')
            return None
        # Sparse stuff
        if plotType!='se':
            tmpX,tmpY=self._sparse(tmpX,tmpY, sparse)

        # Logarithm stuff
        if logy or logx:
            tmpX,tmpY=self._logarithm(tmpX,tmpY,logx,logy,base)

        # Here it ensures that if we are plotting ncycle no values of '*' will be plotted
        tmX=[]
        tmY=[]
        for i in range(len(tmpX)):
            tmX.append(str(tmpX[i]))
            tmY.append(str(tmpY[i]))

        tmpX=[]
        tmpY=[]
        for i in range(len(tmX)):
            if '*' in tmX[i] or '*' in tmY[i]:
                print('There seems to be a string of * in the lists')
                print('Cutting out elements in both the lists that have an index equal to or greater than the index of the location of the string of *')
                break
            tmpX.append(float(tmX[i]))
            tmpY.append(float(tmY[i]))
        listX=tmpX
        listY=tmpY

        #Setting the axis labels

        if logx:
            labelx='log '+labelx
        if logy:
            labely='log '+labely

        if legend!=None:
            legend=legend
        else:
            legend=labely+' vs '+labelx



        pl.plot(listX,listY,shape,label=legend,markevery=markevery,linewidth=linewidth)
        pl.legend()
        pl.title(title)
        pl.xlabel(labelx)
        pl.ylabel(labely)
        if show:
            pl.show()

        if limits != None and len(limits)==4:

            pl.xlim(limits[0],limits[1])
            pl.ylim(limits[2],limits[3])

    def plot_isoratios(self,xiso,yiso,fign=1,spec=None,deltax=True,deltay=True,logx=False,logy=False,
                        title=None,legend=None,legloc='lower right',errbar=True,dcycle=500,addiso=None,
                        co_toggle='c',cust_toggle=None,shift=0,weighting=None,zoneselect=None,iniabufile='iniab2.0E-02GN93.ppn',
                        plt_sparse=1,plt_symb='o',plt_col='k',plt_lt='-',plt_lw=1.,alpha_dum=1.,plt_massrange=False,plt_show=True,
                        figsave=False):
        '''
        This is the new routine to plot isotopic ratios for ALL input. rt, June 2014

        Parameters:
        -----------
        xiso : np.array
            x data to plot. This can be an array or a list of arrays, depending on who calls the routine
        yiso : np.array
            y data to plot. This can be an array or a list of arrays, depending on who calls the routine
        fign : integer, optional
            Figure number
        spec : string, optional
            What specifications do you want to do when coming from nugridse models. Choose 'surf' for
            surface models or 'exp' for explosions (out files)
        deltax : boolean, optional
            X axis in delta values?
        deltay : boolean, optional
            Y axis in delta values?
        logx : boolean, optional
            Logarithmic x axis?
        logy : boolean, optional
            Logarithmic y axis?
        title : string, optional
            Title for your plot
        legend : string, optional
            Legend for your model / grains. For grains the legend is automatically taken from the
            grain class
        legloc : string / integer, optional
            Location of the legend, use matplotlib standard. Use None to not plot legend if plotted
            by default, e.g., from grain class routine.
        errbar : boolean, optional
            Error bars on grain data?
        dcycle : integer, optional
            Difference between cycles to take for thermal pulse searching, if searching is
            deactivated, dcycle describes how often cycles are sampled.  The default is 500.
        addiso : list, optional
            For explosive models. Add an isotope.  Format ['C-12', 0.5 ,'N-12'] to add N12
            to C12 and multiply it with a factor of 0.5.  Multiple isotopes can be added, the
            factor is optional and does not have to be given.  Isotopes can be added to other
            isotopes as well, i.e., [['C-12', 'N-12'], ['C-13', 'N-13']].  The default None.
            Notice that while addiso = [['N-14','O-14'],['N-14',fractionation,'C-14']]
            works, other options like addiso = [['N-14','O-14',fractionation,'C-14']] or
            addiso = [['N-14',1,'O-14',fractionation,'C-14']] are not working and give Typerror.
            CAREFUL, that for the option addiso = [['N-14',fractionation,'C-14','O-14']] there is
            no error message, but the fractionation is applied to both O14 and C14!
        co_toggle : string, optional
            For explosive models, choose what shells you want to look for! Select 'c' for
            selecting zones with C/O >= 1.  Select 'o' for C/O <= 1.  If 'a' takes the
            whole star.  The defalut is 'c'. See cust_toggle (below) for an alternative!
        cust_toggle : list, optional
            This option is like co_toggle (and overwrites it when chosen) but lets you choose
            your own comparison. For example you want to find zones that have a 10 fold
            overabundance of Ti-46 and Ti-47 over O-16 and Zr-96, you can choose here
            [['Ti-46','Ti-47'],['O-16','Zr-96'],100.] Assuming the first list is is x, the
            second list y, and the comparator number is f, the statement only plots shells
            in which the condition x/y>f is fulfilled. x and y are number sums of the chosen
            isotopes, f has to be given as a float. This is only for explosive shells. Please
            note, if this toggle is NOT None, then co_toggle is overwritten!
        shift : integer, optional
            For explosive models, how much do you want to shift the models back from the
            last cycle? By default (0) the last cycle is taken.
        weighting : string, optional
            For explosive models. If None then, plot every profile separately.  If 'zone'
            then, average each zone.  If 'all' then average all selected zones.  The
            default is None.
        zoneselect: string, optional
            For explosive models. Select if you want to plot 'all' zones or outer most zone.
            Arguments are 'all' and 'top', respectively. Default is None, then the user is
            asked to provide this information during the routine as input.
        iniabufile : string, optional
            Initial abundance file. Use absolute path for your file or filename to choose a
            given file in USEEPP. Attention: You need a standard tree checked out from SVN
        plt_sparse : integer, optional
            Every how many datapoints is the plot done? Not used for some routines!
        plt_symb : string, optional
            Symbol for the plot. In case of grains, this is handled automatically.
        plt_col : string / float, optional
            Color for plotted curve. In case of grains, this is handled automatically.
        plt_lt : string, optional
            line type for plot.
        plt_lw : float, optional
            Line width for plot.
        alpha_dum : trasparency to apply to grains data, in case of many data are plotted.
            This may be allpied also for theoretical curves.
        plt_massrange : boolean, optional
            For explosive models. Plot mass of shell with first and last datapoint of
            each zone.  If list given, label those zones.  The default is False.
        plt_show : boolean, optional
            Do you want to show the plot or not?
        figsave : string, optional
            Give path and filename here, if you want to save the figure.
        '''

        from . import utils as u

        ### WORK ON PATH ###
        # define svn path form path where script runs, depending on standard input or not
        if len(iniabufile.split('/')) == 1 :   # means not an absolute path!
            scriptpathtmp = __file__
            if len(scriptpathtmp.split('/')) == 1:   # in folder where nugridse is
                scriptpathtmp = os.path.abspath('.') + '/nugridse.py'   # to get the current dir
            svnpathtmp = '/'
            for i in range(len(scriptpathtmp.split('/'))-3):   # -3 to go to folders up!
                if scriptpathtmp.split('/')[i] != '':
                    svnpathtmp += scriptpathtmp.split('/')[i] + '/'
            iniabufile = svnpathtmp + 'frames/mppnp/USEEPP/' + iniabufile   # make absolute path for iniabufile


        ### get solar system ratios for the isotopes that are specified in the input file ###
        inut = u.iniabu(iniabufile)
        try:
            xrat_solsys = inut.isoratio_init(xiso)
        except KeyError:   # if isotope not available, e.g., if plotting Ti-44 / Ti-48 ratio
            xrat_solsys = 0.
        try:
            yrat_solsys = inut.isoratio_init(yiso)
        except KeyError:
            yrat_solsys = 0.

        # number ratio for solar system ratio
        xrat_solsys *= (old_div(float(xiso[1].split('-')[1]), float(xiso[0].split('-')[1])))
        yrat_solsys *= (old_div(float(yiso[1].split('-')[1]), float(yiso[0].split('-')[1])))


        # initialize xdataerr and ydataerr as None
        xdataerr = None
        ydataerr = None



        ### DO PLOTS FROM NUGRIDSE CLASS ###
        if self._classTest() == 'se':
            if spec==None:
                spec = str(eval(input('Please specify \'surf\' for surface models (AGB stars) or \'exp\' for explosive'
                             'models and zone finding, etc., and press enter: ')))

            ### SURFACE MODELS - PLOT AGB STAR STUFF ###
            if spec == 'surf':
                print('Plotting AGB star stuff')

                # read in thermal pulse position and co ratio
                tp_pos, co_return = self._tp_finder(dcycle)
                tp_pos_tmp = []
                co_return_tmp  = []
                tp_pos_tmp.append(1)
                co_return_tmp.append(co_return[0])
                for i in range(len(tp_pos)):
                    tp_pos_tmp.append(tp_pos[i])
                    co_return_tmp.append(co_return[i])
                tp_pos = tp_pos_tmp
                co_return = co_return_tmp

                # read in data
                iso_alldata = self.get(tp_pos,[xiso[0],xiso[1],yiso[0],yiso[1]])
                xrat = np.zeros(len(iso_alldata))
                yrat = np.zeros(len(iso_alldata))
                for i in range(len(iso_alldata)):
                    xrat[i] = old_div(iso_alldata[i][0], iso_alldata[i][1])
                    yrat[i] = old_div(iso_alldata[i][2], iso_alldata[i][3])

                # make number ratios
                for i in range(len(xrat)):
                    xrat[i] *= old_div(float(xiso[1].split('-')[1]), float(xiso[0].split('-')[1]))
                    yrat[i] *= old_div(float(yiso[1].split('-')[1]), float(yiso[0].split('-')[1]))

                # if delta values are requested, need to calculate those now
                if deltax:
                    xrat = (old_div(xrat,xrat_solsys) -1.) * 1000.
                if deltay:
                    yrat = (old_div(yrat,yrat_solsys) -1.) * 1000.

                # now we might have o and c rich zones. prepare stuff for plotting
                xdata_o = []
                ydata_o = []
                xdata_c = []
                ydata_c = []
                for i in range(len(co_return)):
                    if co_return[i] <= 1.:
                        xdata_o.append(xrat[i])
                        ydata_o.append(yrat[i])
                    else:
                        xdata_c.append(xrat[i])
                        ydata_c.append(yrat[i])
                if xdata_o != [] and xdata_c != []:
                    xdata_o.append(xdata_c[0])
                    ydata_o.append(ydata_c[0])
                # now make the styles
                style_o = [plt_symb + '--', plt_col, '1.', '4', '2',None]
                style_c = [plt_symb + plt_lt,  plt_col, plt_col, '7.', '1', legend]

                # now make data for plotting
                xdata = []
                ydata = []
                style = []
                if xdata_o != []:
                    xdata.append(xdata_o)
                    ydata.append(ydata_o)
                    style.append(style_o)
                if xdata_c != []:
                    xdata.append(xdata_c)
                    ydata.append(ydata_c)
                    style.append(style_c)




            ### EXPLOSIVE MODELS ###
            elif spec == 'exp':
                print('explosive models')
                # compatibility
                co_toggle = co_toggle.lower()
                isotope_list = [xiso[0],xiso[1],yiso[0],yiso[1]]

                # cycle
                cyc_no = self.se.cycles[len(self.se.cycles)-1-shift]
                mco_data = self.get(cyc_no,['mass','C-12','C-13','O-16','O-17','O-18',xiso[0],xiso[1],yiso[0],yiso[1]])
                mass = mco_data[0]
                # if no custom toggle for enrichment
                if cust_toggle == None:
                    c_elem = mco_data[1]+mco_data[2]
                    o_elem = mco_data[3]+mco_data[4]+mco_data[5]
                    co_ratio = c_elem / o_elem * (old_div(16., 12.))
                    co_comp_val = 1.
                else:
                    co_data1 = self.get(cyc_no,cust_toggle[0])
                    co_data2 = self.get(cyc_no,cust_toggle[1])
                    for i in range(len(co_data1)):
                        if i == 0:
                            c_elem = co_data1[i]
                        else:
                            c_elem += co_data1[i]
                    for i in range(len(co_data2)):
                        if i == 0:
                            o_elem = co_data2[i]
                        else:
                            o_elem += co_data2[i]
                    # now we need to make the mass number of everything in here to make number ratios
                    massn1 = 0.
                    for i in range(len(co_data1)):
                        massn1 += sum(co_data1[i]) * float(cust_toggle[0][i].split('-')[1])
                    massn1 /= sum(c_elem)
                    massn2 = 0.
                    for i in range(len(co_data2)):
                        massn2 += sum(co_data2[i]) * float(cust_toggle[1][i].split('-')[1])
                    massn2 /= sum(o_elem)
                    co_ratio = c_elem / o_elem * (old_div(massn2, massn1))   # this has nothing to do with a C/O ratio anymore! but keep name
                    # comparator value
                    co_comp_val = float(cust_toggle[2])

                # get the data now
                isotope_profile = []
                for i in range(6,10):   # in mco_data
                    isotope_profile.append(mco_data[i])

                # add radioactive isotopes (if given)
                if addiso != None:
                    if type(addiso[0] == list):   # then list of lists
                        for i in range(len(addiso)):
                            for j in range(4):
                                if isotope_list[j] == addiso[i][0]:
                                    multiplicator_addiso = 1.
                                    try:
                                        multiplicator_addiso = float(addiso[i][1])
                                        starter = 2
                                    except ValueError:
                                        starter = 1
                                    for k in range(starter,len(addiso[i])):
                                        isotope_profile[j] += array(self.get(cyc_no,addiso[i][k])) * multiplicator_addiso
                    else:
                        for j in range(4):
                            if isotope_list[j] == addiso[0]:
                                multiplicator_addiso = 1.
                                try:
                                    multiplicator_addiso = float(addiso[1])
                                    starter = 2
                                except ValueError:
                                    starter = 1
                                for k in range(starter,len(addiso)):
                                    isotope_profile[j] += array(self.get(cyc_no,addiso[k])) * multiplicator_addiso

                # search for carbon / oxygen rich layers
                crich = []   # alternating start stop values. if odd number, then surface is c-rich, but add stop number
                dumb = True
                if cust_toggle != None:
                    for i in range(len(co_ratio)):
                        if dumb:
                            if co_ratio[i] >= co_comp_val:
                                crich.append(i)
                                dumb = False
                                continue
                        else:
                            if co_ratio[i] < co_comp_val:
                                crich.append(i)
                                dumb = True
                elif co_toggle != 'a':
                    for i in range(len(co_ratio)):
                        if co_toggle == 'c':   # carbon rich
                            if dumb:
                                if co_ratio[i] >= 1:
                                    crich.append(i)
                                    dumb = False
                                    continue
                            else:
                                if co_ratio[i] < 1:
                                    crich.append(i)
                                    dumb = True
                        elif co_toggle == 'o':   # oxygen rich
                            if dumb:
                                if co_ratio[i] <= 1:
                                    crich.append(i)
                                    dumb = False
                                    continue
                            else:
                                if co_ratio[i] > 1:
                                    crich.append(i)
                                    dumb = True
                        else:
                            print('Select your enrichment!')
                            return None
                else:   # take whole star
                    print('Using all profiles to mix')
                    crich.append(0)
                    crich.append(len(co_toggle))

                if len(crich)%2 == 1:
                    crich.append(len(co_ratio)-1)

                if len(crich) == 0:
                    print('Star did not get rich in C or O, depending on what you specified')
                    return None

                # make isotope_profile into array and transpose
                isotope_profile = array(isotope_profile).transpose()


                # Ask user which zones to use
                if co_toggle != 'a':
                    if cust_toggle != None:
                        print('\n\nI have found the following zones:\n')
                    elif co_toggle == 'c':
                        print('\n\nI have found the following carbon rich zones:\n')
                    elif co_toggle == 'o':
                        print('\n\nI have found the following oxygen rich zones:\n')

                    mass_tmp = zeros((len(crich)))
                    for i in range(len(crich)):
                        mass_tmp[i] = mass[crich[i]]

                    j = 1
                    for i in range(old_div(len(crich),2)):
                        print('Mass range (' + str(j) + '):\t' + str(mass_tmp[2*i]) + ' - ' + str(mass_tmp[2*i+1]))
                        j += 1
                    print('\n')
                    if zoneselect == 'all':
                        usr_zones = 0
                    elif zoneselect == 'top':
                        usr_zones = [j-1]
                    else:
                        usr_zones = eval(input('Please select which mass range you want to use. Select 0 for all zones. Otherwise give one zone or a list of zones separated by comma (e.g.: 1, 2, 4): '))

                    crich_dumb = crich
                    if usr_zones == 0:
                        print('I continue w/ all zones then')
                    elif type(usr_zones) == int:   # only one zone selected
                        tmp = int(usr_zones)-1
                        crich = crich_dumb[2*tmp:2*tmp+2]
                    else:
                        crich = []
                        for i in range(len(usr_zones)):
                            tmp = int(usr_zones[i])-1
                            crich.append(crich_dumb[2*tmp])
                            crich.append(crich_dumb[2*tmp + 1])

                # weight profiles according to weighting factor using the selected crich
                # define isos_to_use variable for later
                if weighting == None:
                    isos_to_use = []
                    for i in range(old_div(len(crich),2)):
                        isos_dumb = []
                        n = crich[2*i]
                        while n <= crich[2*i+1]:
                            isos_dumb.append(isotope_profile[n])
                            n += 1
                        isos_to_use.append(array(isos_dumb))


                elif weighting.lower() == 'zone' or weighting.lower() == 'zones':
                    # make array w/ mass weigted isotope ratio (4) for all mass zones
                    isotope_profile_cweight = zeros((old_div(len(crich),2),4))
                    mass_tot = []
                    for i in range(len(isotope_profile_cweight)):   # 2*i is start, 2*i+1 is stop value
                        if crich[2*i] == 0:
                            print('C- / O-rich in first shell (core).')
                        else:
                            dumb = crich[2*i + 1]
                            j = crich[2*i]
                            mass_tmp = 0
                            while j <= dumb:
                                mass_shell = mass[j] - mass[j-1]
                                mass_tmp += mass_shell
                                for k in range(4):
                                    isotope_profile_cweight[i][k] += isotope_profile[j][k]*mass_shell
                                j += 1
                            mass_tot.append(mass_tmp)
                    for i in range(len(isotope_profile_cweight)):
                        for j in range(4):
                            isotope_profile_cweight[i][j] /= mass_tot[i]
                    isos_to_use = [array(isotope_profile_cweight)]

                elif weighting.lower() == 'all':   # average all zones by mass
                    isos_tmp = zeros((1, len(isotope_profile[0])))
                    for i in range(len(isotope_profile)-1):   # neglect surface effects
                        for j in range(len(isos_tmp[0])):
                            mass_shell = mass[i+1] - mass[i]
                            isos_tmp[0][j] += isotope_profile[i][j]*mass_shell
                    # weight all
                    isos_tmp /= sum(mass)
                    isos_to_use = [isos_tmp]


                # change to isotope numbers from mass!
                for i in range(len(isos_to_use)):
                    for j in range(len(isos_to_use[i])):
                        for k in range(len(isos_to_use[i][j])):
                            # here we just divide 'iso_massf' output with the mass number
                            # this means that in the end, the isotope ratios in number space are correc
                            # but have to use ratios from here on for meaningful stuff
                            isos_to_use[i][j][k] /= float(isotope_list[k].split('-')[1])

                # do the ratios and stuff
                ratiox = []
                ratioy = []
                for i in range(len(isos_to_use)):
                    ratiox_dumb = []
                    ratioy_dumb = []
                    for j in range(len(isos_to_use[i])):
                        ratiox_dumb.append(old_div(isos_to_use[i][j][0], isos_to_use[i][j][1]))
                        ratioy_dumb.append(old_div(isos_to_use[i][j][2], isos_to_use[i][j][3]))
                    ratiox.append(array(ratiox_dumb))
                    ratioy.append(array(ratioy_dumb))

                # make arrays for ratiox and ratioy
                ratiox = array(ratiox)
                ratioy = array(ratioy)

                # make number ratio out of everything
                for i in range(len(ratiox)):
                    for j in range(len(ratiox[i])):
                        ratiox[i][j] *= (old_div(float(xiso[1].split('-')[1]), float(xiso[0].split('-')[1])))
                        ratioy[i][j] *= (old_div(float(yiso[1].split('-')[1]), float(yiso[0].split('-')[1])))

                if deltax:
                    ratiox_tmp = []
                    for i in range(len(ratiox)):
                        ratiox_tmp_tmp = []
                        for j in range(len(ratiox[i])):
                            ratiox_tmp_tmp.append((old_div(ratiox[i][j], xrat_solsys) - 1.) * 1000.)
                        ratiox_tmp.append(ratiox_tmp_tmp)
                    ratiox = array(ratiox_tmp)
                if deltay:
                    ratioy_tmp = []
                    for i in range(len(ratioy)):
                        ratioy_tmp_tmp = []
                        for j in range(len(ratioy[i])):
                            ratioy_tmp_tmp.append((old_div(ratioy[i][j], yrat_solsys) - 1.) * 1000.)
                        ratioy_tmp.append(ratioy_tmp_tmp)
                    ratioy = array(ratioy_tmp)

                # create massrange array if necessary
                plt_massrange_lst = []
                if plt_massrange==True:   # use == True because otherwise the list enters here too... why?
                    for i in range(len(ratiox)):
                        plt_massrange_lst.append([ratiox[i][0], ratioy[i][0], mass[crich[2*i]]])
                        plt_massrange_lst.append([ratiox[i][len(ratiox[i])-1], ratioy[i][len(ratioy[i])-1], mass[crich[2*i+1]]])   # start: x-ratio, y-ratio, mass label

                elif plt_massrange != False:
                    for plt_mr_val in plt_massrange:
                        mrng_i = 0
                        while plt_mr_val > mass[mrng_i] and mrng_i < len(mass):
                            mrng_i += 1
                        if mrng_i > 0:
                            mrng_i -= 1
                        mratx = old_div(isotope_profile[mrng_i][0],isotope_profile[mrng_i][1])
                        mraty = old_div(isotope_profile[mrng_i][2],isotope_profile[mrng_i][3])
                        if deltax:
                            mratx = (old_div(mratx, ratiox_solsys) - 1.) * 1000.
                        if deltay:
                            mraty = (old_div(mraty, ratioy_solsys) - 1.) * 1000.
                        plt_massrange_lst.append([mratx,mraty,mass[mrng_i]])

                # make style and prepare for plotting here
                xdata = ratiox
                ydata = ratioy
                style_tmp0= [plt_symb + plt_lt,  plt_col, plt_col, '13.', '1', legend]
                style_tmp = [plt_symb + plt_lt,  plt_col, plt_col, '13.', '1', None]
                style = []
                for i in range(len(xdata)):
                    if i == 0:
                        style.append(style_tmp0)
                    else:
                        style.append(style_tmp)

            else:
                print('You did not specify a useful spec argument -> abort.')
                return None

        ### PLOTS FROM GRAIN CLASS ###
        if self._classTest() == 'grain':
            print('Presolar grains are cool!')
            xdata,xdataerr,ydata,ydataerr,style = self.plot_ratio_return(xiso,yiso,deltax,deltay)
            legend=True
            plt_sparse=1   # to avoid monkey input
            plt_lw = 0.

        ### PLOT ###
        # data is prepared now, make the plots. data must be in format
        # [[data1],[data2],[data3],...]
        # three arrays like this, for xdata, ydata, and style.
        # style format: symbol, edge color, face color, symbol size, edge width, label
        # this is then compatible with grain.py style definitions

        # Size of font etc.
        params = {'axes.labelsize':  20,
                  'text.fontsize':   14,
                  'legend.fontsize': 14,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14}
        pl.rcParams.update(params)

        pl.figure(fign)
        for i in range(len(xdata)):
            if errbar:
                if xdataerr != None or ydataerr != None:
                    pl.errorbar(xdata[i],ydata[i],xerr=xdataerr[i],yerr=ydataerr[i],marker=style[i][0],
                                color=style[i][1],linestyle='',lw=2,markevery=plt_sparse,alpha=alpha_dum)
            pl.plot(xdata[i],ydata[i],style[i][0],c=style[i][1],mfc=style[i][2],ms=float(style[i][3]),
                    mew=float(style[i][4]),label=style[i][5],markevery=plt_sparse,linewidth=plt_lw,alpha=alpha_dum)

        # plot text labels if necessary
        if plt_massrange != False:
            for mrng_ind in range(len(plt_massrange_lst)):
                pl.text(plt_massrange_lst[mrng_ind][0], plt_massrange_lst[mrng_ind][1],
                        str(round(plt_massrange_lst[mrng_ind][2], 2)), ha='right', va='bottom', color=plt_col,fontsize=15.)


        # log?
        if logx and logy == False:
            pl.semilogx()
        elif logx == False and logy:
            pl.semilogy()
        elif logx and logy:
            pl.loglog()

        # legend
        if legend != None and legloc != None:
            pl.legend(loc=legloc)

        # title and labels
        if title != None:
            pl.title(title)
        if deltax:
            pl.xlabel('$\delta$($^{' + xiso[0].split('-')[1] + '}$' +xiso[0].split('-')[0] + '/$^{' + xiso[1].split('-')[1] + '}$' +xiso[1].split('-')[0] + ')' )
        else:
            pl.xlabel('$^{' + xiso[0].split('-')[1] + '}$' +xiso[0].split('-')[0] + '/$^{' + xiso[1].split('-')[1] + '}$' +xiso[1].split('-')[0])
        if deltay:
            pl.ylabel('$\delta$($^{' + yiso[0].split('-')[1] + '}$' +yiso[0].split('-')[0] + '/$^{' + yiso[1].split('-')[1] + '}$' +yiso[1].split('-')[0] + ')' )
        else:
            pl.ylabel('$^{' + yiso[0].split('-')[1] + '}$' +yiso[0].split('-')[0] + '/$^{' + yiso[1].split('-')[1] + '}$' +yiso[1].split('-')[0])

        # plot horizontal and vertical lines
        print(xrat_solsys, yrat_solsys)
        if deltay:
            pl.axhline(0,color='k')
        else:
            pl.axhline(yrat_solsys,color='k')
        if deltax:
            pl.axvline(0,color='k')
        else:
            pl.axvline(xrat_solsys,color='k')

        # borders of plot
        pl.gcf().subplots_adjust(bottom=0.15)
        pl.gcf().subplots_adjust(left=0.15)

        # save and show
        if figsave != False:
            pl.savefig(figsave)

        if plt_show:
            pl.show()

    def plot_isopattern(self,isos,normiso,spec=None,tpulse='all',dcycle=500,fign=1,deltay=False,logy=False,
                        iniabufile='iniab2.0E-02GN93.ppn',legend=None,plt_symb='o',plt_col='k',plt_lt='-',
                        plt_lw=1.,plt_show=True):
        '''
        This routine plots isotopic pattern plots for different input along with grain data.

        Parameters
        ----------
        isos : list / string
            Enter the list of isotopes that you want to consider here. The list should
            be given in the standard format, e.g., ['Fe-54','Fe-56','Fe-57','Fe-58'] or
            give the element as a string if you want to consider all stable isotopes,
            e.g., 'Fe'
        normiso : string
            Give the isotope all the ratios should be normalized to here, e.g., 'Fe-56'
        spec : string, optional
            What specifications do you want to do when coming from nugridse models. Choose 'surf' for
            surface models or 'exp' for explosions (out files)
        tpulse : string, optional
            In case you have an AGB star, here you decide which thermal pulse to plot. You can choose
            'all' (default) to plot all TPs, 'c' or 'o' for all Carbon (C/O > 1) or all Oxygen
            (C/O < 1) rich, respectively, or 'last' for the last thermal pulse only
        dcycle : integer, optional
            Difference between cycles to take for thermal pulse searching, if searching is
            deactivated, dcycle describes how often cycles are sampled.  The default is 500.
        fign : integer, optional
            Number of the figure
        deltay : boolean, optional
            Do you want to do delta values on y axis or regular ratios?
        logy : boolean, optional
            Y axis logarithmic?
        iniabufile : string, optional
            Initial abundance file. Use absolute path for your file or filename to choose a
            given file in USEEPP. Attention: You need a standard tree checked out from SVN
        legend : string, optional
            Legend for your model / grains. For grains the legend is automatically taken from the
            grain class
        plt_symb : string, optional
            Symbol for the plot. In case of grains, this is handled automatically.
        plt_col : string / float, optional
            Color for plotted curve. In case of grains, this is handled automatically.
        plt_lt : string, optional
            line type for plot.
        plt_lw : float, optional
            Line width for plot.
        plt_show : boolean, optional
            Show plot?
        '''

        from . import utils as u

        ### WORK ON PATH ###
        # define svn path form path where script runs, depending on standard input or not
        if len(iniabufile.split('/')) == 1 :   # means not an absolute path!
            scriptpathtmp = __file__
            if len(scriptpathtmp.split('/')) == 1:   # in folder where nugridse is
                scriptpathtmp = os.path.abspath('.') + '/nugridse.py'   # to get the current dir
            svnpathtmp = '/'
            for i in range(len(scriptpathtmp.split('/'))-3):   # -3 to go to folders up!
                if scriptpathtmp.split('/')[i] != '':
                    svnpathtmp += scriptpathtmp.split('/')[i] + '/'
            iniabufile = svnpathtmp + 'frames/mppnp/USEEPP/' + iniabufile   # make absolute path for iniabufile


        ### make a list of isotopes in case an element is specified and not a list of isotopes ###
        if type(isos) == str:
            tmp = u.iniabu.stable_el
            tmp2 = -1
            for i in range(len(tmp)):
                if isos == tmp[i][0]:
                    tmp2 = tmp[i]
                    break

            if tmp2==-1:
                print('No valid element selected. Abort.')
                return None

            isos_tmp = []
            for i in range(1,len(tmp2)):
                isos_tmp.append([isos + '-' + str(int(tmp2[i])),normiso])
            isos = isos_tmp
            # isos is now a list in the sense of [['Fe-54','Fe-56'],['Fe-56','Fe-56'],...]
        elif type(isos) == list:
            isos_tmp = []
            for i in range(len(isos)):
                isos_tmp.append([isos[i], normiso])
            isos = isos_tmp
        else:
            print('Specify a valid isotope, see docstring.')
            return None

        # make a list with all isotopes just as a list (not the fraction as isos)
        isoslist = []
        for i in range(len(isos)):
            isoslist.append(isos[i][0])

        # find out where normisotope sits
        for i in range(len(isos)):
            if isos[i][0] == isos[i][1]:
                posnorm = i
                break

        ### get solar system ratios for the isotopes that are specified in the input file ###
        inut = u.iniabu(iniabufile)
        ss_rat = []
        for i in range(len(isos)):
            ss_rat.append(inut.isoratio_init(isos[i]))
        # make number ratios
        for i in range(len(isos)):
            ss_rat[i] *= (old_div(float(isos[i][1].split('-')[1]), float(isos[i][0].split('-')[1])))

        ### DO PLOTS FROM NUGRIDSE CLASS ###
        if self._classTest() == 'se':
            if spec==None:
                spec = str(eval(input('Please specify \'surf\' for surface models (AGB stars) or \'exp\' for explosive'
                             'models and zone finding, etc., and press enter: ')))

            ### SURFACE MODELS - PLOT AGB STAR STUFF ###
            if spec == 'surf':
                print('Plotting AGB star stuff')

                # read in thermal pulse position and co ratio
                tp_pos, co_return = self._tp_finder(dcycle)
                tp_pos_tmp = []
                co_return_tmp  = []
                tp_pos_tmp.append(1)
                co_return_tmp.append(co_return[0])
                for i in range(len(tp_pos)):
                    tp_pos_tmp.append(tp_pos[i])
                    co_return_tmp.append(co_return[i])
                tp_pos = tp_pos_tmp
                co_return = co_return_tmp

                # read in data
                iso_alldata = self.get(tp_pos,isoslist)
                norm_isotope = np.zeros(len(iso_alldata))
                for i in range(len(iso_alldata)):
                    norm_isotope[i] = iso_alldata[i][posnorm]
                # now make ratios iso_ratios out of iso_alldata
                iso_ratios = np.zeros((len(iso_alldata),len(iso_alldata[0])))
                for i in range(len(iso_ratios)):
                    for j in range(len(iso_ratios[i])):
                        iso_ratios[i][j] = old_div(iso_alldata[i][j], norm_isotope[i])

                # make number ratios
                for i in range(len(iso_ratios)):
                    for j in range(len(iso_ratios[i])):
                        iso_ratios[i][j] *= old_div(float(isos[j][1].split('-')[1]), float(isos[j][0].split('-')[1]))

                # if delta values are requested, need to calculate those now
                if deltay:
                    for i in range(len(iso_ratios)):
                        for j in range(len(iso_ratios[i])):
                            iso_ratios[i][j] = (old_div(iso_ratios[i][j], ss_rat[j]) - 1.) * 1000.

                ### now prepare data to plot ###

                 # make style list
                style = []

                if tpulse == 'last':
                    ydata = [iso_ratios[len(iso_ratios)-1]]
                    style = [[plt_symb + plt_lt,  plt_col, plt_col, '7.', '1', legend]]
                elif tpulse == 'o':
                    ydata = []
                    style = []
                    for i in range(len(co_return)):
                        if co_return[i] < 1:
                            ydata.append(iso_ratios[i])
                            style.append([plt_symb + plt_lt,  plt_col, plt_col, '7.', '1', legend])
                    if ydata == []:
                        print('No O rich thermal pulses found.')
                        return None
                else:   # this means carbon rich, either only or all, but make marker size first
                    msizelst = []   # marker size list for all subsequent ones
                    crich_list = []
                    for i in range(len(co_return)):
                        if co_return[i] >= 1:
                            crich_list.append(co_return[i])
                    crich_max = np.max(crich_list)
                    crich_min = np.min(crich_list)
                    slope_tmp = old_div(9., (crich_max-crich_min))
                    b_tmp = 3. - slope_tmp * crich_min
                    for i in range(len(crich_list)):
                        msizelst.append(crich_list[i] * slope_tmp + b_tmp)
                    # make ydata and style
                    ydata = []
                    style = []
                    j=0
                    for i in range(len(co_return)):
                        if co_return[i] < 1.:
                            j += 1
                            if tpulse == 'all':
                                ydata.append(iso_ratios[i])
                                style.append([plt_symb + '--', plt_col, plt_col, 1., 1., None])
                        else:
                            ydata.append(iso_ratios[i])
                            style.append([plt_symb + plt_lt, plt_col, plt_col, msizelst[i-j], 1., legend])

            ### EXPLOSIVE MODELS ###
            elif spec=='exp':
                print('Explosive models not yet implemented... sorry')
                return None
            else:
                print('You did not specify a useful spec argument -> abort.')
                return None

        ###### PLOT ######

        ### make ratios for data to plot, first find where the normiso sits ###
        # make x axis vector
        xdata = []
        for i in range(len(isos)):
            xdata.append(int(isos[i][0].split('-')[1]))

        # Size of font etc.
        params = {'axes.labelsize':  20,
                  'text.fontsize':   14,
                  'legend.fontsize': 14,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14}
        pl.rcParams.update(params)

        # plot
        pl.figure(fign)
        for i in range(len(ydata)):
            pl.plot(xdata,ydata[i],style[i][0],c=style[i][1],mfc=style[i][2],ms=float(style[i][3]),
                    label=style[i][5],linewidth=plt_lw,mew=float(style[i][4]))

        # limits and x axis scale and handling
        pl.xlim([xdata[0]-0.5,xdata[len(xdata)-1]+0.5])

        # labels and axis scaling
        pl.xlabel('Mass number')
        if deltay:
            pl.ylabel('$\delta (^i$' + normiso.split('-')[0] + ' / solar)')
        else:
            pl.ylabel('$^{i}$'+ normiso.split('-')[0] + ' / solar')

        # borders of plot
        pl.gcf().subplots_adjust(bottom=0.15)
        pl.gcf().subplots_adjust(left=0.15)

        # log
        if logy:
            pl.semilogy()

        if plt_show:
            pl.show()


    def _clear(self, title=True, xlabel=True, ylabel=True):
        '''
        Method for removing the title and/or xlabel and/or Ylabel.

        Parameters
        ----------
        Title : boolean, optional
            Boolean of if title will be cleared.  The default is True.
        xlabel : boolean, optional
            Boolean of if xlabel will be cleared.  The default is True.
        ylabel : boolean, optional
            Boolean of if ylabel will be cleared.  The default is True.

        '''
        if title:
            pyl.title('')
        if xlabel:
            pyl.xlabel('')
        if ylabel:
            pyl.ylabel('')

    # From mesa.py
    def _xlimrev(self):
        ''' reverse xrange'''
        xmax,xmin=pyl.xlim()
        pyl.xlim(xmin,xmax)

    def abu_chartMulti(self, cyclist, mass_range=None, ilabel=True,
                       imlabel=True, imlabel_fontsize=12, imagic=False,
                       boxstable=True, lbound=20, plotaxis=[0,0,0,0],
                       color_map='jet', pdf=False, title=None, path=None):
        '''
        Method that plots abundence chart and saves those figures to a
        .png file (by default). Plots a figure for each cycle in the
        argument cycle.

        Parameters
        ----------
        cyclist : list
            The list of cycles we are plotting.
        mass_range : list, optional
            A 1x2 array containing the lower and upper mass range.  If
            this is an instance of abu_vector this will only plot
            isotopes that have an atomic mass within this range.  This
            will throw an error if this range does not make sence ie
            [45,2] if None, it will plot over the entire range.  The
            default is None.
        ilabel : boolean, optional
            Elemental labels off/on.  The default is True.
        imlabel : boolean, optional
            Label for isotopic masses off/on.  The efault is True.
        imlabel_fontsize : intager, optional
            Fontsize for isotopic mass labels.  The default is 12.
        imagic : boolean, optional
            Turn lines for magic numbers off/on.  The default is False.
        boxstable : boolean, optional
            Plot the black boxes around the stable elements.  The
            default is True.
        lbound : tuple, optional
            Boundaries for colour spectrum ploted.  The defaults is 20.
        plotaxis : list, optional
            Set axis limit: If [0,0,0,0] the complete range in (N,Z)
            will be plotted.  The default is [0,0,0,0].
        color_map : string, optional
            Color map according to choices in matplotlib
            (e.g. www.scipy.org/Cookbook/Matplotlib/Show_colormaps).
            The default is 'jet'.
        pdf : boolean, optional
            What format will this be saved in pdf/png.  The default is
            True.
        title : string, optional
            The title of the plots and the saved images.  The default is
            None.
        '''

        if self._which('dvipng')==None:
            print("This method may need the third party program dvipng to operate")
            print('It is located at http://sourceforge.net/projects/dvipng/')

        max_num = max(cyclist)
        for i in range(len(cyclist)):
            self.abu_chart( cyclist[i], mass_range ,ilabel,imlabel,imlabel_fontsize,imagic,\
                            boxstable,lbound,plotaxis,False,color_map)
            if title !=None:
                pl.title(title)
            else:
                name='AbuChart'
            if path is not None:
                name = os.path.join(path, name)
            number_str=_padding_model_number(cyclist[i],max_num)
            if not pdf:
                pl.savefig(name+number_str+'.png', dpi=100)
            else:
                pl.savefig(name+number_str+'.pdf', dpi=200)
            pl.close()

        return None

    #from mppnp.se
    def abu_chart(self, cycle, mass_range=None ,ilabel=True,
                  imlabel=True, imlabel_fontsize=12, imagic=False,
                  boxstable=True, lbound=(-12, 0),
                  plotaxis=[0, 0, 0, 0], show=True, color_map='jet',
                  ifig=None,data_provided=False,thedata=None,
                  savefig=False,drawfig=None,drawax=None,mov=False,
                  path=None):
        '''
        Plots an abundance chart

        Parameters
        ----------
        cycle : string, integer or list
            The cycle we are looking in. If it is a list of cycles,
            this method will then do a plot for each of these cycles
            and save them all to a file.
        mass_range : list, optional
            A 1x2 array containing the lower and upper mass range.  If
            this is an instance of abu_vector this will only plot
            isotopes that have an atomic mass within this range.  This
            will throw an error if this range does not make sence ie
            [45,2] if None, it will plot over the entire range.  The
            default is None.
        ilabel : boolean, optional
            Elemental labels off/on.  The default is True.
        imlabel : boolean, optional
            Label for isotopic masses off/on.  The default is True.
        imlabel_fontsize : integer, optional
            Fontsize for isotopic mass labels.  The default is 12.
        imagic : boolean, optional
            Turn lines for magic numbers off/on.  The default is False.
        boxstable : boolean, optional
            Plot the black boxes around the stable elements.  The
            defaults is True.
        lbound : tuple, optional
            Boundaries for colour spectrum ploted.  The default is
            (-12,0).
        plotaxis : list, optional
            Set axis limit.  If [0, 0, 0, 0] the complete range in (N,Z)
            will be plotted.  It equates to [xMin, xMax, Ymin, Ymax].
            The default is [0, 0, 0, 0].
        show : boolean, optional
            Boolean of if the plot should be displayed.  Useful with
            saving multiple plots using abu_chartMulti.  The default is
            True.
        color_map : string, optional
            Color map according to choices in matplotlib
            (e.g. www.scipy.org/Cookbook/Matplotlib/Show_colormaps).
            The default is 'jet'.
        ifig : integer, optional
            Figure number, if ifig is None it wiil be set to the cycle
            number.  The defaults is None.
        savefig : boolean, optional
            Whether or not to save the figure.
            The default is False
        drawfig, drawax, mov : optional, not necessary for user to set these variables
            The figure and axes containers to be drawn on, and whether or not a movie is
            being made (only True when se.movie is called, which sets mov to True
            automatically
        path: path where to save figure

        '''

        if ifig == None and not mov:
            ifig=cycle

        if type(cycle)==type([]):
            self.abu_chartMulti(cycle, mass_range,ilabel,imlabel,imlabel_fontsize,imagic,boxstable,\
                                lbound,plotaxis,color_map, path=path)
            return
        plotType=self._classTest()

        if mass_range!=None and mass_range[0]>mass_range[1]:
            raise IOError("Please input a proper mass range")

        if plotType=='se':
            if not data_provided:
                cycle=self.se.findCycle(cycle)
#                nin=zeros(len(self.se.A))
#                zin=zeros(len(self.se.Z))
                yin=self.get(cycle, 'iso_massf')
                isom=self.se.isomeric_states
                masses = self.se.get(cycle,'mass')
            else:
                cycle=cycle # why so serious?
                yin=thedata[0]
                isom=self.se.isomeric_states
                masses = thedata[1]

#            for i in xrange(len(nin)):
#                zin[i]=self.se.Z[i]
#                nin[i]=self.se.A[i]-zin[i]
            # SJONES implicit loop instead:
            zin=array([el for el in self.se.Z])
            nin=array([el for el in self.se.A])-zin

            #Test if the mass cell order is inverted
            #and hence mass[-1] the center.
            if masses[0]>masses[-1]:
                #invert
                print('Inverted order of mass cells will be taken into account.')
                yin=yin[::-1]
                masses=masses[::-1]

            if mass_range != None:
                # trim out only the zones needed:
                tmpyps=[]
                masses.sort() # SJ: not sure why this sort if necessary
#                for i in xrange(len(masses)):
#                    if (masses[i] >mass_range[0] and masses[i]<mass_range[1]) or\
#                            (masses[i]==mass_range[0] or masses[i]==mass_range[1]):
#                        tmpyps.append(yin[i])
#                yin=tmpyps
                # find lower and upper indices and slice instead:
                idxl=np.abs(masses-mass_range[0]).argmin()
                if masses[idxl] < mass_range[0]: idxl+=1
                idxu=np.abs(masses-mass_range[1]).argmin()
                if masses[idxu] > mass_range[1]: idxu-=1
                yin=yin[idxl:idxu+1]




            #tmp=zeros(len(yin[0]))
            #for i in xrange(len(yin)):
            #    for j in xrange(len(yin[i])):
            #        tmp[j]+=yin[i][j]

            tmp2=sum(yin,axis=0) # SJONES sum along axis instead of nested loop
            tmp=old_div(tmp2,len(yin))

            yin=tmp

        elif plotType=='PPN':

            ain=self.get('A',cycle)
            zin=self.get('Z',cycle)
            nin=ain-zin
            yin=self.get('ABUNDANCE_MF',cycle)
            isom=self.get('ISOM',cycle)

            if mass_range != None:
                tmpA=[]
                tmpZ=[]
                tmpIsom=[]
                tmpyps=[]
                for i in range(len(nin)):
                    if (ain[i] >mass_range[0] and ain[i]<mass_range[1])\
                    or (ain[i]==mass_range[0] or ain[i]==mass_range[1]):
                        tmpA.append(nin[i])
                        tmpZ.append(zin[i])
                        tmpIsom.append(isom[i])
                        tmpyps.append(yin[i])
                zin=tmpZ
                nin=tmpA
                yin=tmpyps
                isom=tmpIsom

        else:
            raise IOError("This method, abu_chart, is not supported by this class")

        # in case we call from ipython -pylab, turn interactive on at end again
        turnoff=False
        if not show:
            try:
                ioff()
                turnoff=True
            except NameError:
                turnoff=False

        nnmax = int(max(nin))+1
        nzmax = int(max(zin))+1
        nzycheck = zeros([nnmax,nzmax,3])

        for i in range(len(nin)):
            if isom[i]==1:
                ni = int(nin[i])
                zi = int(zin[i])

                nzycheck[ni,zi,0] = 1
                nzycheck[ni,zi,1] = yin[i]

        #######################################################################
        # elemental names: elname(i) is the name of element with Z=i

        elname=self.elements_names

        #### create plot
        ## define axis and plot style (colormap, size, fontsize etc.)
        if plotaxis==[0,0,0,0]:
            xdim=10
            ydim=6
        else:
            dx = plotaxis[1]-plotaxis[0]
            dy = plotaxis[3]-plotaxis[2]
            ydim = 6
            xdim = ydim*dx/dy


        params = {'axes.labelsize':  12,
                  'text.fontsize':   12,
                  'legend.fontsize': 12,
                  'xtick.labelsize': 12,
                  'ytick.labelsize': 12,
                  'text.usetex': True}
        #pl.rcParams.update(params) #May cause Error, someting to do with tex
        if mov:
            fig=drawfig
            fig.set_size_inches(xdim,ydim)
            artists=[]
        else:
            fig=pl.figure(ifig,figsize=(xdim,ydim),dpi=100)
        axx = 0.10
        axy = 0.10
        axw = 0.85
        axh = 0.8
        if mov:
            ax=drawax
        else:
            ax=pl.axes([axx,axy,axw,axh])
        # Tick marks
        xminorlocator = MultipleLocator(1)
        xmajorlocator = MultipleLocator(5)
        ax.xaxis.set_major_locator(xmajorlocator)
        ax.xaxis.set_minor_locator(xminorlocator)
        yminorlocator = MultipleLocator(1)
        ymajorlocator = MultipleLocator(5)
        ax.yaxis.set_major_locator(ymajorlocator)
        ax.yaxis.set_minor_locator(yminorlocator)

        # color map choice for abundances

        cmapa = cm.get_cmap(name=color_map)
        # color map choice for arrows
        cmapr = cm.autumn
        # if a value is below the lower limit its set to white
        cmapa.set_under(color='w')
        cmapr.set_under(color='w')
        # set value range for abundance colors (log10(Y))
        norma = colors.Normalize(vmin=lbound[0],vmax=lbound[1])
        # set x- and y-axis scale aspect ratio to 1
        ax.set_aspect('equal')
        #print time,temp and density on top
        temp = ' '#'%8.3e' %ff['temp']
        time = ' '#'%8.3e' %ff['time']
        dens = ' '#'%8.3e' %ff['dens']

        #May cause Error, someting to do with tex
        '''
        #box1 = TextArea("t : " + time + " s~~/~~T$_{9}$ : " + temp + "~~/~~$\\rho_{b}$ : " \
        #             + dens + ' g/cm$^{3}$', textprops=dict(color="k"))
        anchored_box = AnchoredOffsetbox(loc=3,
                        child=box1, pad=0.,
                        frameon=False,
                        bbox_to_anchor=(0., 1.02),
                        bbox_transform=ax.transAxes,
                        borderpad=0.,
                        )
        ax.add_artist(anchored_box)
        '''
        ## Colour bar plotted

        patches = []
        color = []
        for i in range(nzmax):
            for j in range(nnmax):
                if nzycheck[j,i,0]==1:
                    xy = j-0.5,i-0.5

                    rect = Rectangle(xy,1,1,)

                    # abundance
                    yab = nzycheck[j,i,1]
                    if yab == 0:

                        yab=1e-99


                    col =log10(yab)

                    patches.append(rect)
                    color.append(col)

        p = PatchCollection(patches, cmap=cmapa, norm=norma)
        p.set_array(array(color))
        p.set_zorder(1)
        if mov:
            artist1=ax.add_collection(p)
            artists.append(artist1)
        else:
            ax.add_collection(p)
        if not mov:
            cb = pl.colorbar(p)

            # colorbar label
            cb.set_label('log$_{10}$(X)',fontsize='x-large')

        # plot file name
        graphname = 'abundance-chart'+str(cycle)

        # Add black frames for stable isotopes
        if boxstable:
            for i in range(len(self.stable_el)):
                if i == 0:
                    continue


                tmp = self.stable_el[i]
                try:
                    zz= self.elements_names.index(tmp[0]) #charge
                except:
                    continue

                for j in range(len(tmp)):
                    if j == 0:
                        continue

                    nn = int(tmp[j]) #atomic mass
                    nn=nn-zz

                    xy = nn-0.5,zz-0.5
                    rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=3.)
                    rect.set_zorder(2)
                    ax.add_patch(rect)


        # decide which array to take for label positions
        iarr = 0

        # plot element labels
        if ilabel:
            for z in range(nzmax):
                try:
                    nmin = min(argwhere(nzycheck[:,z,iarr]))[0]-1
                    ax.text(nmin,z,elname[z],horizontalalignment='center',verticalalignment='center',\
                            fontsize='x-small',clip_on=True)
                except ValueError:
                    continue

        # plot mass numbers
        if imlabel:
            for z in range(nzmax):
                for n in range(nnmax):
                    a = z+n
                    if nzycheck[n,z,iarr]==1:
                        ax.text(n,z,a,horizontalalignment='center',verticalalignment='center',\
                                fontsize=imlabel_fontsize,clip_on=True)

        # plot lines at magic numbers
        if imagic:
            ixymagic=[2, 8, 20, 28, 50, 82, 126]
            nmagic = len(ixymagic)
            for magic in ixymagic:
                if magic<=nzmax:
                    try:
                        xnmin = min(argwhere(nzycheck[:,magic,iarr]))[0]
                        xnmax = max(argwhere(nzycheck[:,magic,iarr]))[0]
                        line = ax.plot([xnmin,xnmax],[magic,magic],lw=3.,color='r',ls='-')
                    except ValueError:
                        dummy=0
                if magic<=nnmax:
                    try:
                        yzmin = min(argwhere(nzycheck[magic,:,iarr]))[0]
                        yzmax = max(argwhere(nzycheck[magic,:,iarr]))[0]
                        line = ax.plot([magic,magic],[yzmin,yzmax],lw=3.,color='r',ls='-')
                    except ValueError:
                        dummy=0

        # set axis limits
        if plotaxis==[0,0,0,0]:

            xmax=max(nin)
            ymax=max(zin)
            ax.axis([-0.5,xmax+0.5,-0.5,ymax+0.5])
        else:
            ax.axis(plotaxis)

        # set x- and y-axis label
        ax.set_xlabel('neutron number (A-Z)',fontsize=14)
        ax.set_ylabel('proton number Z',fontsize=14)
        if not mov:
            pl.title('Isotopic Chart for cycle '+str(int(cycle)))
        if savefig:
            if path is not None:
                graphname = os.path.join(path, graphname)
            fig.savefig(graphname)
            print(graphname,'is done')
        if show:
            pl.show()
        if turnoff:
            ion()

        if mov:
            return p,artists
        else:
            return

    def abu_flux_chart(self, cycle, ilabel=True, imlabel=True,
                       imagic=False, boxstable=True, lbound=(-12,0),
                       plotaxis=[0,0,0,0], which_flux=None, prange=None,
                       profile='charged', show=True):
        '''
        Plots an abundance and flux chart

        Parameters
        ----------
        cycle : string, integer or list
            The cycle we are looking in. If it is a list of cycles,
            this method will then do a plot for each of these cycles
            and save them all to a file.
        ilabel : boolean, optional
            Elemental labels off/on.  The default is True.
        imlabel : boolean, optional
            Label for isotopic masses off/on.  The default is True.
        imagic : boolean, optional
            Turn lines for magic numbers off/on.  The default is False.
        boxstable : boolean, optional
            Plot the black boxes around the stable elements.  The
            defaults is True.
        lbound : tuple, optional
            Boundaries for colour spectrum ploted.  The default is
            (-12,0).
        plotaxis : list, optional
            Set axis limit.  If [0, 0, 0, 0] the complete range in (N,Z)
            will be plotted.  It equates to [xMin, xMax, Ymin, Ymax].
            The default is [0, 0, 0, 0].
        which_flux : integer, optional
            Set to 0 for nucleosynthesis flux plot.  Set to 1 for
            energy flux plot.  Setting wich_flux to 0 is equivelent to
            setting it to 0.  The default is None.
        prange : integer, optional
            Range of fluxes to be considered, if prange is None then
            the plot range is set to 8.  The default is None.
        profile : string, optional
            'charged' is ideal setting to show charged particle
            reactions flow.  'neutron' is ideal setting for neutron
            captures flows.  The default is 'charged'.
        show : boolean, optional
            Boolean of if the plot should be displayed.  Useful with
            saving multiple plots using abu_chartMulti.  The default is
            True.

        '''
        #######################################################################
        #### plot options
        # Set axis limit: If default [0,0,0,0] the complete range in (N,Z) will
        # be plotted, i.e. all isotopes, else specify the limits in
        # plotaxis = [xmin,xmax,ymin,ymax]

        #######################################################################

        # read data file
        #inpfile = cycle
        #ff = fdic.ff(inpfile)
        # with the flux implementation I am not using mass range for now.
        # It may be introduced eventually.
        mass_range = None
        if str(cycle.__class__)=="<type 'list'>":
            self.abu_chartMulti(cycle, mass_range,ilabel,imlabel,imlabel_fontsize,imagic,boxstable,\
                                lbound,plotaxis)
            return
        plotType=self._classTest()

        #if mass_range!=None and mass_range[0]>mass_range[1]:
            #print 'Please input a proper mass range'
            #print 'Returning None'
            #return None

        if plotType=='se':
            cycle=self.se.findCycle(cycle)
            nin=zeros(len(self.se.A))
            zin=zeros(len(self.se.Z))
            for i in range(len(nin)):
                nin[i]=self.se.A[i]
                zin[i]=self.se.Z[i]
            for i in range(len(nin)):
                nin[i]=nin[i]-zin[i]
            yin=self.get(cycle, 'iso_massf')
            isom=self.se.isomeric_states

            masses = self.se.get(cycle,'mass')
            if mass_range != None:
                masses = self.se.get(cycle,'mass')
                masses.sort()

            if mass_range != None:
                tmpyps=[]
                masses = self.se.get(cycle,'mass')
                masses = self.se.get(cycle,'mass')
                masses.sort()
                for i in range(len(masses)):
                    if (masses[i] >mass_range[0] and masses[i]<mass_range[1]) or\
                            (masses[i]==mass_range[0] or masses[i]==mass_range[1]):
                        tmpyps.append(yin[i])
                yin=tmpyps


            tmp=zeros(len(yin[0]))
            for i in range(len(yin)):
                for j in range(len(yin[i])):
                    tmp[j]+=yin[i][j]

            tmp=old_div(tmp,len(yin))

            yin=tmp

        elif plotType=='PPN':

            ain=self.get('A',cycle)
            zin=self.get('Z',cycle)
            nin=ain-zin
            yin=self.get('ABUNDANCE_MF',cycle)
            isom=self.get('ISOM',cycle)

            if mass_range != None:
                tmpA=[]
                tmpZ=[]
                tmpIsom=[]
                tmpyps=[]
                for i in range(len(nin)):
                    if (ain[i] >mass_range[0] and ain[i]<mass_range[1])\
                    or (ain[i]==mass_range[0] or ain[i]==mass_range[1]):
                        tmpA.append(nin[i])
                        tmpZ.append(zin[i])
                        tmpIsom.append(isom[i])
                        tmpyps.append(yin[i])
                zin=tmpZ
                nin=tmpA
                yin=tmpyps
                isom=tmpIsom

        else:
            print('This method, abu_chart, is not supported by this class')
            print('Returning None')
            return None
        # in case we call from ipython -pylab, turn interactive on at end again
        turnoff=False
        if not show:
            try:
                ioff()
                turnoff=True
            except NameError:
                turnoff=False

        nnmax = int(max(nin))+1
        nzmax = int(max(zin))+1
        nnmax_plot = nnmax
        nzmax_plot = nzmax
        nzycheck = zeros([nnmax,nzmax,3])
        nzycheck_plot = zeros([nnmax,nzmax,3])
        for i in range(len(nin)):
            if isom[i]==1:
                ni = int(nin[i])
                zi = int(zin[i])

                nzycheck[ni,zi,0] = 1
                nzycheck[ni,zi,1] = yin[i]
                nzycheck_plot[ni,zi,0] = 1



        #######################################################################
        # elemental names: elname(i) is the name of element with Z=i

        elname=self.elements_names

        #### create plot
        ## define axis and plot style (colormap, size, fontsize etc.)
        if plotaxis==[0,0,0,0]:
            xdim=10
            ydim=6
        else:
            dx = plotaxis[1]-plotaxis[0]
            dy = plotaxis[3]-plotaxis[2]
            ydim = 6
            xdim = ydim*dx/dy


        params = {'axes.labelsize':  15,
                  'text.fontsize':   12,
                  'legend.fontsize': 15,
                  'xtick.labelsize': 15,
                  'ytick.labelsize': 15,
                  'text.usetex': True}
        #pl.rcParams.update(params) #May cause Error, someting to do with tex
        #fig=pl.figure(figsize=(xdim,ydim),dpi=100)
        fig=pl.figure()
        if profile == 'charged':
            ax1 = fig.add_subplot(1, 2, 1)
        elif profile == 'neutron':
            ax1 = fig.add_subplot(2, 1, 1)
        #axx = 0.10
        #axy = 0.10
        #axw = 0.85
        #axh = 0.8
        #ax1=pl.axes([axx,axy,axw,axh])
        # Tick marks
        xminorlocator = MultipleLocator(1)
        xmajorlocator = MultipleLocator(5)
        ax1.xaxis.set_major_locator(xmajorlocator)
        ax1.xaxis.set_minor_locator(xminorlocator)
        yminorlocator = MultipleLocator(1)
        ymajorlocator = MultipleLocator(5)
        ax1.yaxis.set_major_locator(ymajorlocator)
        ax1.yaxis.set_minor_locator(yminorlocator)

        # color map choice for abundances
        #cmapa = cm.jet
        cmapa = cm.summer
        # color map choice for arrows
        cmapr = cm.summer
        # if a value is below the lower limit its set to white
        cmapa.set_under(color='w')
        cmapr.set_under(color='w')
        # set value range for abundance colors (log10(Y))
        norma = colors.Normalize(vmin=lbound[0],vmax=lbound[1])
        # set x- and y-axis scale aspect ratio to 1
        #ax1.set_aspect('equal')
        #print time,temp and density on top
        temp = ' '#'%8.3e' %ff['temp']
        time = ' '#'%8.3e' %ff['time']
        dens = ' '#'%8.3e' %ff['dens']

        #May cause Error, someting to do with tex
        '''
        #box1 = TextArea("t : " + time + " s~~/~~T$_{9}$ : " + temp + "~~/~~$\\rho_{b}$ : " \
        #             + dens + ' g/cm$^{3}$', textprops=dict(color="k"))
        anchored_box = AnchoredOffsetbox(loc=3,
                        child=box1, pad=0.,
                        frameon=False,
                        bbox_to_anchor=(0., 1.02),
                        bbox_transform=ax.transAxes,
                        borderpad=0.,
                        )
        ax.add_artist(anchored_box)
        '''
        ## Colour bar plotted

        patches = []
        color = []

        for i in range(nzmax):
            for j in range(nnmax):
                if nzycheck[j,i,0]==1:
                    xy = j-0.5,i-0.5

                    rect = Rectangle(xy,1,1,)

                    # abundance
                    yab = nzycheck[j,i,1]
                    if yab == 0:

                        yab=1e-99


                    col =log10(yab)

                    patches.append(rect)
                    color.append(col)


        p = PatchCollection(patches, cmap=cmapa, norm=norma)
        p.set_array(array(color))
        p.set_zorder(1)
        ax1.add_collection(p)
        cb = pl.colorbar(p)

        # colorbar label
        if profile == 'neutron':
            cb.set_label('log$_{10}$(X)',fontsize='x-large')

        # plot file name
        graphname = 'abundance-flux-chart'+str(cycle)

        # Add black frames for stable isotopes
        if boxstable:
            for i in range(len(self.stable_el)):
                if i == 0:
                    continue


                tmp = self.stable_el[i]
                try:
                    zz= self.elements_names.index(tmp[0]) #charge
                except:
                    continue

                for j in range(len(tmp)):
                    if j == 0:
                        continue

                    nn = int(tmp[j]) #atomic mass
                    nn=nn-zz

                    xy = nn-0.5,zz-0.5
                    rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=4.)
                    rect.set_zorder(2)
                    ax1.add_patch(rect)




        # decide which array to take for label positions
        iarr = 0

        # plot element labels
        if ilabel:
            for z in range(nzmax):
                try:
                    nmin = min(argwhere(nzycheck[:,z,iarr]))[0]-1
                    nmax = max(argwhere(nzycheck[:,z,iarr]))[0]+1
                    ax1.text(nmin,z,elname[z],horizontalalignment='center',verticalalignment='center',\
                            fontsize='medium',clip_on=True)
                    ax1.text(nmax,z,elname[z],horizontalalignment='center',verticalalignment='center',\
                            fontsize='medium',clip_on=True)
                except ValueError:
                    continue

        # plot mass numbers
        if imlabel:
            for z in range(nzmax):
                for n in range(nnmax):
                    a = z+n
                    if nzycheck[n,z,iarr]==1:
                        ax1.text(n,z,a,horizontalalignment='center',verticalalignment='center',\
                                fontsize='xx-small',clip_on=True)


        # plot lines at magic numbers
        if imagic:
            ixymagic=[2, 8, 20, 28, 50, 82, 126]
            nmagic = len(ixymagic)
            for magic in ixymagic:
                if magic<=nzmax:
                    try:
                        xnmin = min(argwhere(nzycheck[:,magic,iarr]))[0]
                        xnmax = max(argwhere(nzycheck[:,magic,iarr]))[0]
                        line = ax1.plot([xnmin,xnmax],[magic,magic],lw=3.,color='r',ls='-')
                    except ValueError:
                        dummy=0
                if magic<=nnmax:
                    try:
                        yzmin = min(argwhere(nzycheck[magic,:,iarr]))[0]
                        yzmax = max(argwhere(nzycheck[magic,:,iarr]))[0]
                        line = ax1.plot([magic,magic],[yzmin,yzmax],lw=3.,color='r',ls='-')
                    except ValueError:
                        dummy=0

        # set axis limits
        if plotaxis==[0,0,0,0]:

            xmax=max(nin)
            ymax=max(zin)
            ax1.axis([-0.5,xmax+0.5,-0.5,ymax+0.5])
        else:
            ax1.axis(plotaxis)

        # set x- and y-axis label
        ax1.set_ylabel('Proton number',fontsize='xx-large')
        if profile == 'charged':
            ax1.set_xlabel('Neutron number',fontsize='xx-large')
        #pl.title('Isotopic Chart for cycle '+str(int(cycle)))

        #
        # here below I read data from the flux_*****.DAT file.
        #
        file_name = 'flux_'+str(cycle).zfill(5)+'.DAT'
        print(file_name)
        f = open(file_name)
        lines = f.readline()
        lines = f.readlines()
        f.close()

        print_max_flux_in_plot =  False
        # color map choice for fluxes
        #cmapa = cm.jet
        cmapa = cm.autumn
        # color map choice for arrows
        cmapr = cm.autumn
        # starting point of arrow
        coord_x_1 = []
        coord_y_1 = []
        # ending point of arrow (option 1)
        coord_x_2 = []
        coord_y_2 = []
        # ending point of arrow (option 2)
        coord_x_3 = []
        coord_y_3 = []
        # fluxes
        flux_read = []
        flux_log10 = []

        if which_flux == None or which_flux == 0:
            print('chart for nucleosynthesis fluxes [dYi/dt]')
            line_to_read = 9
        elif which_flux == 1:
            print('chart for energy fluxes')
            line_to_read = 10
        elif which_flux > 1:
            print("you have only option 0 or 1, not larger than 1")

        single_line = []
        for i in range(len(lines)):
            single_line.append(lines[i].split())
            coord_y_1.append(int(single_line[i][1]))
            coord_x_1.append(int(single_line[i][2])-coord_y_1[i])
            coord_y_2.append(int(single_line[i][5]))
            coord_x_2.append(int(single_line[i][6])-coord_y_2[i])
            coord_y_3.append(int(single_line[i][7]))
            coord_x_3.append(int(single_line[i][8])-coord_y_3[i])
            try:
                flux_read.append(float(single_line[i][line_to_read]))
            except ValueError: # this is done to avoid format issues like 3.13725-181...
                flux_read.append(1.0E-99)
            flux_log10.append(log10(flux_read[i]+1.0e-99))

        print(file_name,' read!')

        # I need to select smaller sample, with only fluxes inside plotaxis.
        if plotaxis!=[0,0,0,0]:
            coord_y_1_small=[]
            coord_x_1_small=[]
            coord_y_2_small=[]
            coord_x_2_small=[]
            coord_y_3_small=[]
            coord_x_3_small=[]
            flux_log10_small = []
            for i in range(len(flux_log10)):
                I_am_in = 0
                if coord_y_1[i] > plotaxis[2] and coord_y_1[i] < plotaxis[3] and coord_x_1[i] > plotaxis[0] and coord_x_1[i] < plotaxis[1]:
                    I_am_in = 1
                    coord_y_1_small.append(int(coord_y_1[i]))
                    coord_x_1_small.append(int(coord_x_1[i]))
                    coord_y_2_small.append(int(coord_y_2[i]))
                    coord_x_2_small.append(int(coord_x_2[i]))
                    coord_y_3_small.append(int(coord_y_3[i]))
                    coord_x_3_small.append(int(coord_x_3[i]))
                    flux_log10_small.append(flux_log10[i])
                if coord_y_3[i] > plotaxis[2] and coord_y_3[i] < plotaxis[3] and coord_x_3[i] > plotaxis[0] and coord_x_3[i] < plotaxis[1] and I_am_in == 0:
                    I_am_in = 1
                    coord_y_1_small.append(int(coord_y_1[i]))
                    coord_x_1_small.append(int(coord_x_1[i]))
                    coord_y_2_small.append(int(coord_y_2[i]))
                    coord_x_2_small.append(int(coord_x_2[i]))
                    coord_y_3_small.append(int(coord_y_3[i]))
                    coord_x_3_small.append(int(coord_x_3[i]))
                    flux_log10_small.append(flux_log10[i])



        # elemental labels off/on [0/1]
        ilabel = 1

        # label for isotopic masses off/on [0/1]
        imlabel = 1

        # turn lines for magic numbers off/on [0/1]
        imagic = 0

        # flow is plotted over "prange" dex. If flow < maxflow-prange it is not plotted
        if prange == None:
            print('plot range given by default')
            prange = 8.

        #############################################
        #print flux_log10_small
        # we should scale prange on plot_axis range, not on max_flux!
        max_flux = max(flux_log10)
        ind_max_flux = flux_log10.index(max_flux)
        if plotaxis!=[0,0,0,0]:
            max_flux_small = max(flux_log10_small)

        if plotaxis==[0,0,0,0]:
            nzmax = int(max(max(coord_y_1),max(coord_y_2),max(coord_y_3)))+1
            nnmax = int(max(max(coord_x_1),max(coord_x_2),max(coord_x_3)))+1
            coord_x_1_small = coord_x_1
            coord_x_2_small = coord_x_2
            coord_x_3_small = coord_x_3
            coord_y_1_small = coord_y_1
            coord_y_2_small = coord_y_2
            coord_y_3_small = coord_y_3
            flux_log10_small= flux_log10
            max_flux_small  = max_flux
        else:
            nzmax = int(max(max(coord_y_1_small),max(coord_y_2_small),max(coord_y_3_small)))+1
            nnmax = int(max(max(coord_x_1_small),max(coord_x_2_small),max(coord_x_3_small)))+1

        for i in range(nzmax):
            for j in range(nnmax):
                if nzycheck[j,i,0]==1:
                    xy = j-0.5,i-0.5
                    rect = Rectangle(xy,1,1,)
                    patches.append(rect)


        nzycheck = zeros([nnmax_plot,nzmax,3])
        coord_x_out = zeros(len(coord_x_2_small), dtype='int')
        coord_y_out = zeros(len(coord_y_2_small),dtype='int')
        for i in range(len(flux_log10_small)):
            nzycheck[coord_x_1_small[i],coord_y_1_small[i],0] = 1
            nzycheck[coord_x_1_small[i],coord_y_1_small[i],1] = flux_log10_small[i]
            if coord_x_2_small[i] >= coord_x_3_small[i]:
                coord_x_out[i] = coord_x_2_small[i]
                coord_y_out[i] = coord_y_2_small[i]
                nzycheck[coord_x_out[i],coord_y_out[i],0] = 1
                nzycheck[coord_x_out[i],coord_y_out[i],1] = flux_log10_small[i]
            elif coord_x_2_small[i] < coord_x_3_small[i]:
                coord_x_out[i] = coord_x_3_small[i]
                coord_y_out[i] = coord_y_3_small[i]
                nzycheck[coord_x_out[i],coord_y_out[i],0] = 1
                nzycheck[coord_x_out[i],coord_y_out[i],1] = flux_log10_small[i]
            if flux_log10_small[i]>max_flux_small-prange:
                nzycheck[coord_x_1_small[i],coord_y_1_small[i],2] = 1
                nzycheck[coord_x_out[i],coord_y_out[i],2] = 1

        #### create plot
        if profile == 'charged':
            ax2 = fig.add_subplot(1, 2, 2)
        elif profile == 'neutron':
            ax2 = fig.add_subplot(2, 1, 2)
        # Tick marks
        xminorlocator = MultipleLocator(1)
        xmajorlocator = MultipleLocator(5)
        ax2.xaxis.set_major_locator(xmajorlocator)
        ax2.xaxis.set_minor_locator(xminorlocator)
        yminorlocator = MultipleLocator(1)
        ymajorlocator = MultipleLocator(5)
        ax2.yaxis.set_major_locator(ymajorlocator)
        ax2.yaxis.set_minor_locator(yminorlocator)
        ## define axis and plot style (colormap, size, fontsize etc.)
        if plotaxis==[0,0,0,0]:
            xdim=10
            ydim=6
        else:
            dx = plotaxis[1]-plotaxis[0]
            dy = plotaxis[3]-plotaxis[2]
            ydim = 6
            xdim = ydim*dx/dy

        format = 'pdf'
        # set x- and y-axis scale aspect ratio to 1
        #ax2.set_aspect('equal')

        # Add black frames for stable isotopes
        # Add black frames for stable isotopes
        if boxstable:
            for i in range(len(self.stable_el)):
                if i == 0:
                    continue


                tmp = self.stable_el[i]
                try:
                    zz= self.elements_names.index(tmp[0]) #charge
                except:
                    continue

                for j in range(len(tmp)):
                    if j == 0:
                        continue

                    nn = int(tmp[j]) #atomic mass
                    nn=nn-zz

                    xy = nn-0.5,zz-0.5
                    rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=4.)
                    rect.set_zorder(2)
                    ax2.add_patch(rect)


        apatches = []
        acolor = []
        m = old_div(0.8,prange)
        vmax=ceil(max(flux_log10_small))
        vmin=max(flux_log10_small)-prange
        b=-vmin*m+0.1
        normr = colors.Normalize(vmin=vmin,vmax=vmax)
        ymax=0.
        xmax=0.

        for i in range(len(flux_log10_small)):
            x = coord_x_1_small[i]
            y = coord_y_1_small[i]
            dx = coord_x_out[i]-coord_x_1_small[i]
            dy = coord_y_out[i]-coord_y_1_small[i]
            if flux_log10_small[i]>=vmin:
                arrowwidth = flux_log10_small[i]*m+b
                arrow = Arrow(x,y,dx,dy, width=arrowwidth)
                if xmax<x:
                    xmax=x
                if ymax<y:
                    ymax=y
                acol = flux_log10_small[i]
                apatches.append(arrow)
                acolor.append(acol)
            xy = x-0.5,y-0.5
            rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=1.)
            patches.append(rect)
            xy = x+dx-0.5,y+dy-0.5
            rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=1.)
            patches.append(rect)


        p = PatchCollection(patches,norm=0,facecolor='w')
        p.set_zorder(1)
        ax2.add_collection(p)




        a = PatchCollection(apatches, cmap=cmapr, norm=normr)
        a.set_array(array(acolor))
        a.set_zorder(3)
        ax2.add_collection(a)
        cb = pl.colorbar(a)

        # colorbar label
        cb.set_label('log$_{10}$($x$)',fontsize='x-large')
        if profile == 'neutron':
            cb.set_label('log$_{10}$(f)',fontsize='x-large')

        # decide which array to take for label positions
        iarr = 2

        # plot element labels
        for z in range(nzmax):
            try:
                nmin = min(argwhere(nzycheck_plot[:,z,iarr-2]))[0]-1
                nmax = max(argwhere(nzycheck_plot[:,z,iarr-2]))[0]+1
                ax2.text(nmin,z,elname[z],horizontalalignment='center',verticalalignment='center',fontsize='medium',clip_on=True)
                ax2.text(nmax,z,elname[z],horizontalalignment='center',verticalalignment='center',fontsize='medium',clip_on=True)
            except ValueError:
                continue

        # plot mass numbers
        if imlabel:
            for z in range(nzmax):
                for n in range(nnmax_plot):
                    a = z+n
                    if nzycheck_plot[n,z,iarr-2]==1:
                        ax2.text(n,z,a,horizontalalignment='center',verticalalignment='center',fontsize='xx-small',clip_on=True)

        # plot lines at magic numbers
        if imagic==1:
            ixymagic=[2, 8, 20, 28, 50, 82, 126]
            nmagic = len(ixymagic)
            for magic in ixymagic:
                if magic<=nzmax:
                    try:
                        xnmin = min(argwhere(nzycheck[:,magic,iarr-2]))[0]
                        xnmax = max(argwhere(nzycheck[:,magic,iarr-2]))[0]
                        line = ax2.plot([xnmin,xnmax],[magic,magic],lw=3.,color='r',ls='-')
                    except ValueError:
                        dummy=0
                if magic<=nnmax:
                    try:
                        yzmin = min(argwhere(nzycheck[magic,:,iarr-2]))[0]
                        yzmax = max(argwhere(nzycheck[magic,:,iarr-2]))[0]
                        line = ax2.plot([magic,magic],[yzmin,yzmax],lw=3.,color='r',ls='-')
                    except ValueError:
                        dummy=0

        # set axis limits
        if plotaxis==[0,0,0,0]:
            ax2.axis([-0.5,xmax+0.5,-0.5,ymax+0.5])
        else:
            ax2.axis(plotaxis)

        # set x- and y-axis label
        ax2.set_xlabel('Neutron number',fontsize='xx-large')
        if profile == 'neutron':
            ax2.set_ylabel('Proton number',fontsize='xx-large')
        if which_flux == None or which_flux == 0:
            max_flux_label="max flux = "+str('{0:.4f}'.format(max_flux))
        elif which_flux == 1:
            max_flux_label="max energy flux = "+str('{0:.4f}'.format(max_flux))
        if print_max_flux_in_plot:
            ax2.text(plotaxis[1]-1.8,plotaxis[2]+0.1,max_flux_label,fontsize=10.)

        #fig.savefig(graphname)
        print(graphname,'is done')
        if show:
            pl.show()
        if turnoff:
            ion()
        return

    def iso_abundMulti(self, cyclist, stable=False, amass_range=None,
                       mass_range=None, ylim=[0,0], ref=-1,
                       decayed=False, include_title=False, title=None,
                       pdf=False, color_plot=True, grid=False,
                       point_set=1):
        '''
        Method that plots figures and saves those figures to a .png
        file.  Plots a figure for each cycle in the argument cycle.
        Can be called via iso_abund method by passing a list to cycle.

        Parameters
        ----------
        cycllist : list
            The cycles of interest.  This method will do a plot for
            each cycle and save them to a file.
        stable : boolean, optional
            A boolean of whether to filter out the unstables.  The
            defaults is False.
        amass_range : list, optional
            A 1x2 array containing the lower and upper atomic mass
            range.  If None plot entire available atomic mass range.
            The default is None.
        mass_range : list, optional
            A 1x2 array containing the lower and upper mass range.  If
            this is an instance of abu_vector this will only plot
            isotopes that have an atominc mass within this range.  This
            will throw an error if this range does not make sense ie
            [45,2].  If None, it will plot over the entire range.  The
            defaults is None.
        ylim : list, optional
            A 1x2 array containing the lower and upper Y limits.  If
            it is [0,0], then ylim will be determined automatically.
            The default is [0,0].
        ref : integer or list, optional
            reference cycle.  If it is not -1, this method will plot
            the abundences of cycle devided by the cycle of the same
            instance given in the ref variable.  If ref is a list it
            will be interpreted to have two elements:
            ref=['dir/of/ref/run',cycle] which uses a refernece cycle
            from another run.  If any abundence in the reference cycle
            is zero, it will replace it with 1e-99.  The default is -1.
        decayed : boolean, optional
            If True plot decayed distributions, else plot life
            distribution.  The default is False.
        include_title : boolean, optional
            Include a title with the plot.  The default is False.
        title : string, optional
            A title to include with the plot.  The default is None.
        pdf : boolean, optional
            Save image as a [pdf/png].  The default is False.
        color_plot : boolean, optional
            Color dots and lines [True/False].  The default is True.
        grid : boolean, optional
            print grid.  The default is False.
        point_set : integer, optional
            Set to 0, 1 or 2 to select one of three point sets, useful
            for multiple abundances or ratios in one plot.  The defalult
            is 1.

        '''
        max_num = max(cyclist)
        for i in range(len(cyclist)):
            self.iso_abund(cyclist[i],stable,amass_range,mass_range,ylim,ref,\
                 decayed=decayed,show=False,color_plot=color_plot,grid=False,\
                 point_set=1,include_title=include_title)
            if title !=None:
                pl.title(title)
            else:
                name='IsoAbund'
            number_str=_padding_model_number(cyclist[i],max_num)
            if not pdf:
                pl.savefig(name+number_str+'.png', dpi=200)
            else:
                pl.savefig(name+number_str+'.pdf', dpi=200)
            pl.clf()

        return None


    def elem_abund(self, cycle, Z_range=None, mass_range=None,
                   ylim=[0,0], ref=-1, decayed=False, bracket=False,
                   solar_file='./iniab2.0E-02GN93.ppn', ref_ZZ=1, fancy=False):
        '''
        plot the abundance of all elements

        WARNING: Pavel and Falk spend an hour 150205 looking into this method
                 and we convinced our case that at least for the mulit-zone
                 i-process case that we looked at (mppnp-hif example) this method
                 does not give the right plots. We think that elemental_abund does
                 work, though. But neither does provide decayed elemental plots yet.


        Parameters
        ----------
        cycle : string or integer
            The cycle of interest.
        Z_range : list, optional
            A 1x2 array containing the lower and upper proton number
            range. If None plot entire available proton number range.
            The default is None.
        mass_range : list, optional
            A 1x2 array containing the lower and upper mass range. If
            None it will plot over the entire mass range.
            The default is None.
        ylim : list, optional
            A 1x2 array containing the lower and upper Y limits. If
            [0,0] then ylim will be determined automatically.
            The default is [0,0].
        ref : integer, optional
            reference cycle.  If it is not -1, this method will plot
            the abundances of cycle divided by either the solar abundances
            (ref=0) or divided by the cycle of the same
            instance given in the ref variable (ref>0). If any abundence in
            the reference cycle is zero, it will replace it with 1e-99.
            The default is -1.
        decayed : boolean, optional
            If True plot decayed distributions, else plot live
            distribution.
            The default is False.
        bracket : boolean, optional
            If True, plot bracket notation [X/ZZ], which requires
            knowledge of the solar abundance distribution from an
            iniab file and the user to input which ZZ
            The default is False
        solar_file : string, optional
            Path to the file containing the solar abundance distribution.
            The default is './iniab2.0E-02GN93.ppn'
        ref_ZZ : integer, optional
            Proton number of denominator in bracket notation plot.
            The default is 1.
        fancy : boolean, optional
            whether or not to try and use Times New Roman for the font.
            The default is False
        shape : string
            linestyle string
        '''

        if fancy:
            fsize=18
            try:
                params = {'axes.labelsize':  fsize,
                        'text.fontsize':   fsize,
                        'legend.fontsize': fsize,
                        'xtick.labelsize': fsize,
                        'ytick.labelsize': fsize,
                        'text.usetex': False,
                        'font.family': 'Times New Roman',
                        'figure.facecolor': 'white',
                        'ytick.minor.pad': 8,
                        'ytick.major.pad': 8,
                        'xtick.minor.pad': 8,
                        'xtick.major.pad': 8,
                        'figure.subplot.bottom' : 0.15,
                        'lines.markersize': 6}
                pl.rcParams.update(params)
                elemfont = {'family' : 'sans-serif'}
            except:
                pass

        plotType=self._classTest()
        if plotType=='se':
            if decayed:
                raise IOError("Decayed not supported yet! Sorry!")

            def get_av_elem(cycle):
                mass=self.se.get(cycle,'mass')
                print(mass)
                if mass_range is not None:
                    idx1=abs(mass-mass_range[0]).argmin()
                    idx2=abs(mass-mass_range[1]).argmin()
                    if idx1>idx2:
                        mass=mass[idx2:idx1]
                        ea=self.get_elemental_abunds(cycle,index=[idx2,idx1])
                    else:
                        mass=mass[idx1:idx2]
                        ea=self.get_elemental_abunds(cycle,index=[idx1,idx2])
                else:
                    ea=self.get_elemental_abunds(cycle)
                # mass-average the abundances:
                dm=diff(np.insert(mass,0,0.))
                totmass=sum(dm)
                mea=[array(ea[i])*dm[i] for i in range(len(mass))]
                eamsum=zeros(len(ea[0]))
                for i in range(len(eamsum)):
                    eamsum[i]=sum([zone[i] for zone in mea])
                ea=old_div(eamsum,totmass)
                return np.maximum(ea,1.e-99)

            Z=array(self.se.Z)
            Zuq=array(list(set(Z))) # unique list of Z
            Zuq.sort()
            names=self.se.isos[1:] # abolish the neutron

            ea=get_av_elem(cycle)
            if ref>0:
                ea=old_div(ea,get_av_elem(ref))
            elif ref==0 or bracket:
                # initialise solar data
                # If this import statment goes at the top of the file
                # then it will cause a circular import loop witch will
                # cause all modules that try to import data_plot to crash.
                from . import utils as u
                u.solar(solar_file,1)
                selem=u.solar_elem_abund
                selem=np.maximum(selem,1.e-99)
                Zsol=array(u.z_sol)
                Zuqsol=array(list(set(Zsol)))
                Zuqsol.sort()
                # get rid of MPPNP data not available for the Sun:
                rm=[]
                for ZZ in Zuq:
                    if ZZ not in Zuqsol:
                        rm.append(ZZ)

                names=delete(names,[where(Zuq==rmz)[0] for rmz in rm])
                ea=delete(ea,[where(Zuq==rmz)[0] for rmz in rm])
                Zuq=delete(Zuq,[where(Zuq==rmz)[0] for rmz in rm])

                # get rid of solar data not available from MPPNP:
                rm=[]
                for ZZ in Zuqsol:
                    if ZZ not in Zuq:
                        rm.append(ZZ)

                selem=delete(selem,[where(Zuqsol==rmz)[0] for rmz in rm])

                if bracket:
                    xrefsol=selem[where(Zuqsol==ref_ZZ)[0]]
                    xrefsol=np.maximum(xrefsol,1.e-99)
                    xref=ea[where(Zuq==ref_ZZ)[0]]
                    xref=np.maximum(xref,1.e-99)
                    selem=old_div(selem,xrefsol)
                    ea=old_div(ea,xref)

                ea=old_div(ea,selem)

        else:
            raise IOError("elem_abund currently unsupported for this data type:")
            print(plotType)

        pl.plot(Zuq,np.log10(ea),'ro')
        pl.plot(Zuq,np.log10(ea),'r-')
        for i in range(len(Zuq)):
            loc=(-1)**(i%2)*.5
            if fancy:
                pl.text(Zuq[i],np.log10(ea[i])+loc,names[i],fontsize=8,horizontalalignment='center',verticalalignment='center',clip_on=True,fontdict=elemfont)
            else:
                pl.text(Zuq[i],np.log10(ea[i])+loc,names[i],fontsize=8,horizontalalignment='center',verticalalignment='center',clip_on=True)

        if ref>0:
            pl.ylabel('$\log\,X('+str(cycle)+')/X('+str(ref)+')$')
        elif ref==0:
            pl.ylabel('$\log\,X/X_\odot$')
        elif bracket:
            dname=names[where(Zuq==ref_ZZ)[0]][0]
            pl.ylabel('$[X/{\\rm '+dname+'}]$')
        else:
            pl.ylabel('$\log\,X$')
        if ylim!=[0,0]:
            pl.ylim((ylim[0],ylim[1]))
        if Z_range is not None:
            pl.xlim((Z_range[0],Z_range[1]))
        pl.xlabel('$Z$')
        pl.show()


    def iso_abund(self, cycle, stable=False, amass_range=None,
                  mass_range=None, ylim=[0,0], ref=-1, show=True,
                  log_logic=True, decayed=False, color_plot=True,
                  grid=False, point_set=1, include_title=False,
                  data_provided=False,thedata=None, verbose=True,
                  mov=False,drawfig=None,drawax=None,show_names=True,
                  label=None,colour=None,elemaburtn=False,mypoint=None):
        '''
        plot the abundance of all the chemical species

        Parameters
        ----------
        cycle : string, integer or list
            The cycle of interest.  If it is a list of cycles, this
            method will do a plot for each cycle and save them to a
            file.
        stable : boolean, optional
            A boolean of whether to filter out the unstables.  The
            defaults is False.
        amass_range : list, optional
            A 1x2 array containing the lower and upper atomic mass
            range.  If None plot entire available atomic mass range.
            The default is None.
        mass_range : list, optional
            A 1x2 array containing the lower and upper mass range.  If
            this is an instance of abu_vector this will only plot
            isotopes that have an atominc mass within this range.  This
            will throw an error if this range does not make sense ie
            [45,2].  If None, it will plot over the entire range.  The
            defaults is None.
        ylim : list, optional
            A 1x2 array containing the lower and upper Y limits.  If
            it is [0,0], then ylim will be determined automatically.
            The default is [0,0].
        ref : integer or list, optional
            reference cycle.  If it is not -1, this method will plot
            the abundences of cycle devided by the cycle of the same
            instance given in the ref variable.  If ref is a list it
            will be interpreted to have two elements:
            ref=['dir/of/ref/run',cycle] which uses a refernece cycle
            from another run.  If any abundence in the reference cycle
            is zero, it will replace it with 1e-99.  The default is -1.
        show : boolean, optional
            Boolean of if the plot should be displayed.  The default is
            True.
        log_logic : boolean, optional
            Plot abundances in log scale or linear.  The default is
            True.
        decayed : boolean, optional
            If True plot decayed distributions, else plot life
            distribution.  The default is False.
        color_plot : boolean, optional
            Color dots and lines [True/False].  The default is True.
        grid : boolean, optional
            print grid.  The default is False.
        point_set : integer, optional
            Set to 0, 1 or 2 to select one of three point sets, useful
            for multiple abundances or ratios in one plot.  The defalult
            is 1.
        include_title : boolean, optional
            Include a title with the plot.  The default is False.
        drawfig, drawax, mov : optional, not necessary for user to set these variables
            The figure and axes containers to be drawn on, and whether or not a movie is
            being made (only True when se.movie is called, which sets mov to True
            automatically
        elemaburtn : boolean, private
            If true, iso_abund() returns after writing self.***_iso_to_plot for
            use with other plotting routines.
        mypoint : string, optional
            fix the marker style of all the points in this plot to one type, given
            as a string. If None, multiple point styles are used as per point_set.
            The default is None
        '''

        plotType=self._classTest()
        if str(cycle.__class__)=="<type 'list'>":
            self.iso_abundMulti(cycle, stable,amass_range,mass_range,ylim,ref,
                 decayed,include_title,color_plot=color_plot,grid=False,point_set=point_set)
            return

        if mass_range!=None and mass_range[0]>mass_range[1]:
            print('Please input a proper mass range')
            print('Returning None')
            return None
        if amass_range!=None and amass_range[0]>amass_range[1]:
            print('Please input a proper Atomic mass range')
            print('Returning None')
            return None
        if plotType=='se':
            if decayed:
                print('Decay option not yet implemented for mppnp - but it is easy do! Consider investing the time!')
                return None


            # get things as arrays
            if not data_provided:
                cycle=self.se.findCycle(cycle)
                a_iso_to_plot   = array(self.se.A)
                abunds          = self.get(cycle,'iso_massf')
                isotope_to_plot = array(self.se.isotopes)
                z_iso_to_plot   = array(self.se.Z)
                isomers_to_plot = array(self.se.isomeric_states)
                if ref >-1:
                    ref=self.se.findCycle(ref)
                    abundsRef=self.se.get(ref,'iso_massf')
                masses = self.se.get(cycle,'mass')
            else:
                cycle=cycle # why so serious?
                a_iso_to_plot   = array(self.se.A)
                abunds          = thedata[0]
                isotope_to_plot = array(self.se.isotopes)
                z_iso_to_plot   = array(self.se.Z)
                isomers_to_plot = array(self.se.isomeric_states)

                if ref >-1:
                    raise IOError("No. It's not ready yet.")
                    #ref=self.se.findCycle(ref)
                    #abundsRef=self.se.get(ref,'iso_massf')
                masses = thedata[1]

            if mass_range == None:
                if verbose:
                    print('Using default mass range')
                mass_range = [min(masses),max(masses)]
            masses.sort()
            mass_range.sort()

            if amass_range == None:
                amass_range=[int(min(a_iso_to_plot)),int(max(a_iso_to_plot))]

            # remove neutrons - this could move in the non- se/PPN specific part below
            if 0 in z_iso_to_plot:
                ind_neut        = where(z_iso_to_plot==0)[0][0]
                a_iso_to_plot   = delete(a_iso_to_plot,ind_neut)
                z_iso_to_plot   = delete(z_iso_to_plot,ind_neut)
                isomers_to_plot = delete(isomers_to_plot,ind_neut)
                isotope_to_plot = delete(isotope_to_plot,ind_neut)
                abunds = delete(abunds,ind_neut,1)
                if ref >-1:
                    abundsRef = delete(abundsRef,ind_neut,1)

            # extract amass_range
            acon=(a_iso_to_plot>=amass_range[0]) & (a_iso_to_plot<=amass_range[1])
            isomers_to_plot = isomers_to_plot[acon]
            isotope_to_plot = isotope_to_plot[acon]
            z_iso_to_plot   = z_iso_to_plot[acon]
            abunds          = abunds.T[acon].T
            if ref >-1:
                abundsRef = abundsRef.T[acon].T
            a_iso_to_plot   = a_iso_to_plot[acon]
            el_iso_to_plot  = array([x.split('-')[0] for x in isotope_to_plot.tolist()])
            # apply mass range
            if mass_range == None:
                if verbose:
                    print('Using default mass range')
                mass_range = [min(masses),max(masses)]
            mass_range.sort()
            aabs = []
            if ref >-1:
                cyc  = [cycle,ref]
                abus = [abunds,abundsRef]
            else:
                cyc  = [cycle]
                abus = [abunds]
            for cc,aa in zip(cyc,abus):
                if not data_provided:
                    masses = self.se.get(cc,'mass')
                else:
                    masses=masses # why so serious?
                masses.sort()
                dmass  = masses[1:] - masses[:-1]    # I should check the grid definition
                dmass  = append(dmass,0.)
                mcon   = (masses>=mass_range[0]) & (masses<=mass_range[1])
                dmass  = dmass[mcon]
                aa = aa[mcon]

                # average over mass range:
                aa = (aa.T*dmass).T.sum(0)
                aa = old_div(aa, (mass_range[1] - mass_range[0]))
                # abunds has now length of isotope_to_plot
                aabs.append(aa)
            if ref >-1:
                abunds = old_div(aabs[0],(aabs[1]+1.e-99))
            else:
                abunds = aabs[0]

            self.a_iso_to_plot=a_iso_to_plot
            self.isotope_to_plot=isotope_to_plot
            self.z_iso_to_plot=z_iso_to_plot
            self.el_iso_to_plot=el_iso_to_plot
            self.abunds=abunds
            self.isomers_to_plot=isomers_to_plot

            if elemaburtn: return
#                        self.isotopes = self.se.isotopes

        elif plotType=='PPN':
            print("This method adds the following variables to the instance:")
            print("a_iso_to_plot      mass number of plotted range of species")
            print("isotope_to_plot    corresponding list of isotopes")
            print("z_iso_to_plot      corresponding charge numbers")
            print("el_iso_to_plot     corresponding element names")
            print("abunds             corresponding abundances")
            print("isom               isomers and their abundance")

            self.get(cycle,decayed=decayed)
            if ref is not -1:
                if type(ref) is list: # reference cycle from other run
                    import ppn
                    pp=ppn.abu_vector(ref[0])
                    abunds_pp=pp.get(ref[1],decayed=decayed)
                    self.abunds=old_div(self.abunds,pp.abunds)
                else:
                    abunds=self.abunds
                    self.get(ref,decayed=decayed)
                    self.abunds=old_div(abunds,(self.abunds+1.e-99))
            if amass_range == None:
                amass_range=[min(self.a_iso_to_plot),max(self.a_iso_to_plot)]

            aa=ma.masked_outside(self.a_iso_to_plot,amass_range[0],amass_range[1])
            isotope_to_plot=ma.array(self.isotope_to_plot,mask=aa.mask).compressed()
            z_iso_to_plot=ma.array(self.z_iso_to_plot,mask=aa.mask).compressed()
            el_iso_to_plot=ma.array(self.el_iso_to_plot,mask=aa.mask).compressed()
            abunds=ma.array(self.abunds,mask=aa.mask).compressed()
            a_iso_to_plot=aa.compressed()
            isomers_to_plot=[]
            for i in range(len(self.isom)):
                if int(self.isom[i][0].split('-')[1])>100:
                    isomers_to_plot.append(self.isom[i])

            self.a_iso_to_plot=a_iso_to_plot
            self.isotope_to_plot=isotope_to_plot
            self.z_iso_to_plot=z_iso_to_plot
            self.el_iso_to_plot=el_iso_to_plot
            self.abunds=abunds
            self.isomers_to_plot=isomers_to_plot
        else:
            print('This method, iso_abund, is not supported by this class')
            print('Returning None')
            return None

        if verbose:
            print('Using the following conditions:')
            if plotType=='se':
                print('\tmass_range:', mass_range[0], mass_range[1])
            print('\tAtomic mass_range:', amass_range[0], amass_range[1])
            print('\tcycle:           ',cycle)
            print('\tplot only stable:',stable)
            print('\tplot decayed:    ',decayed)


        if stable: # remove unstables:
            # For the element that belongs to the isotope at index 5 in isotope_to_plot
            # (C-12) the following gives the mass numbers of stable elements:
            # self.stable_el[self.stable_names.index(el_iso_to_plot[5])][1:]
            ind_delete=[]
            for i in range(len(isotope_to_plot)):
                if a_iso_to_plot[i] not in self.stable_el[self.stable_names.index(el_iso_to_plot[i])][1:]:
                    ind_delete.append(i)
            a_iso_to_plot   = delete(a_iso_to_plot,  ind_delete)
            z_iso_to_plot   = delete(z_iso_to_plot,  ind_delete)
            isomers_to_plot = delete(isomers_to_plot,ind_delete)
            isotope_to_plot = delete(isotope_to_plot,ind_delete)
            el_iso_to_plot  = delete(el_iso_to_plot, ind_delete)
            abunds          = delete(abunds,         ind_delete)

#        el_list=[] # list of elements in el_iso_to_plot
#
#        for el in self.elements_names:
#            if el in el_iso_to_plot:
#                el_list.append(el)

        # SJONES implicit loop:
        el_list = [el for el in self.elements_names if el in el_iso_to_plot]

        abund_plot = []  # extract for each element an abundance and associated
        mass_num   = []  # mass number array, sorted by mass number

        for el in el_list:
            numbers = a_iso_to_plot[(el_iso_to_plot==el)]
            abund_plot.append(abunds[(el_iso_to_plot==el)][argsort(numbers)])
            mass_num.append(sort(numbers))
        # now plot:
        plot_type = ['-','--','-.',':','-']
        pl_index = 0
        if mypoint is None:
            points = [['o','^','p','h','*'],['x','+','D','>','s'],['H','v','<','*','3']]
        else:
            points = [ [mypoint]*5 , [mypoint]*5 , [mypoint]*5]
        if color_plot:
            colors = ['g','r','c','m','k']
        elif colour is not None:
            colors = [colour]*5
        else:
            colors = ['k','k','k','k','k']

        ylim1 =  1.e99
        ylim2 = -1.e99

        # initialise movie-related things:
        if mov:
            artists=[]
            ax=drawax
            fig=drawfig
        elif drawax is not None:
            ax=drawax
        else:
            ax=pl.axes()

        if drawfig is not None:
            fig=drawfig

        for j in range(len(abund_plot)):        #Loop through the elements of interest
#            for l in xrange(len(abund_plot[j])):
#                if abund_plot[j][l] == 0:
#                    abund_plot[j][l] = 1e-99

            abund_plot[j] = np.maximum(abund_plot[j],1.e-99) # SJONES instead of looping

#            a_dum=zeros(len(abund_plot[j]))   # this I (FH) have to do because for some
            if log_logic == False:            # reason log10(abu_abund[j]) does not work
                a_dum = abund_plot[j]         # although abu_abund[j] is a numpy array?!?
            else:
#                for ii in range(len(abund_plot[j])):
#                    a_dum[ii]=log10(abund_plot[j][ii])
                a_dum=np.log10(abund_plot[j]) # SJONES this seems to work fine for me
            if type(colors[0]) is str:
                this_label=str(colors[pl_index]+points[point_set][pl_index]+\
                           plot_type[pl_index])
            else:
                this_label=None
            if mov:
                artist1,=ax.plot(mass_num[j],a_dum,this_label,markersize=6,
                                 markeredgecolor='None')
            else:
                if this_label is not None:
                    if label is not None and j==0:
                        pl.plot(mass_num[j],a_dum,this_label,markersize=6,
                                label=label,markeredgecolor='None')
                        pl.legend(loc='best').draw_frame(False)
                    else:
                        pl.plot(mass_num[j],a_dum,this_label,markersize=6,
                                markeredgecolor='None')
                else:
                    if label is not None and j==0:
                        pl.plot(mass_num[j],a_dum,
                                color=colors[pl_index],
                                marker=points[point_set][pl_index],
                                linestyle=plot_type[pl_index],
                                markersize=6,label=label,
                                markeredgecolor='None')
                        pl.legend(loc='best').draw_frame(False)
                    else:
                        pl.plot(mass_num[j],a_dum,
                                color=colors[pl_index],
                                marker=points[point_set][pl_index],
                                linestyle=plot_type[pl_index],
                                markersize=6,markeredgecolor='None')

            abu_max = max(a_dum)
            max_index=where(a_dum==abu_max)[0][0]
            coordinates=[mass_num[j][max_index],abu_max]
            if mov:
                artist2=ax.text(coordinates[0]+0.1,1.05*coordinates[1],el_list[j],clip_on=True)
            else:
                if show_names:
#                    pl.text(coordinates[0]+0.1,1.05*coordinates[1],el_list[j],clip_on=True)
                    pl.text(coordinates[0],np.log10(2.2*10.**coordinates[1]),
                            el_list[j],clip_on=True,
                            horizontalalignment='center')

            pl_index+=1
            if pl_index > 4:
                pl_index = 0
            ylim1=min(ylim1,min(a_dum))
            ylim2=max(ylim2,max(a_dum))

            if mov:
                artists.extend([artist1,artist2])

        # now trimming the ylims
        if log_logic:
            dylim=0.05*(ylim2-ylim1)
            ylim1 = ylim1 -dylim
            ylim2 = ylim2 +dylim
            if ref is not -1:
                ylim2 = min(ylim2,4)
                ylim1 = max(ylim1,-4)
            else:
                ylim2 = min(ylim2,0.2)
                ylim1 = max(ylim1,-13)
        else:
            ylim1 = ylim1 *0.8
            ylim2 = ylim2 *1.1
        if include_title:
            if plotType=='se':
                if ref == -1:
                    title = str('Range %4.2f' %mass_range[0]) + str('-%4.2f' %mass_range[1]) +\
                        str(' for cycle %d' %int(cycle))
                else:
                    title = str('Range %4.2f' %mass_range[0]) + \
                        str('-%4.2f' %mass_range[1]) + str(' for cycle %d' %int(cycle))+\
                        str(' relative to cycle %d'  %int(ref))
            else:
                if ref == -1:
                    title = str('Cycle %d' %int(cycle))
                else:
                    title = str('Cycle %d' %int(cycle))+\
                            str(' relative to cycle %d'  %int(ref))
            print("including title: ...")
            if mov:
                artist1,=ax.title(title)
                artists.append(artist1)
            else:
                pl.title(title)
        if ylim[0] == 0 and ylim[1] == 0:
            pl.ylim(ylim1,ylim2)
        else:
            pl.ylim(ylim[0],ylim[1])
        pl.xlim([amass_range[0]-.5,amass_range[1]+.5])
        pl.xlabel('mass number (A)',fontsize=14)
        if ref is not -1:
            if log_logic:
                pl.ylabel(r'log abundance ratio',fontsize=14)
            else:
                pl.ylabel(r'abundance ratio',fontsize=14)
        else:
            if log_logic:
                pl.ylabel(r'log mass fraction ',fontsize=14)
            else:
                pl.ylabel(r'mass fraction',fontsize=14)


        if amass_range != None:
            minimum_mass = amass_range[0]
            maximum_mass = amass_range[1]

        elif mass_range != None:
            minimum_mass = mass_range[0]
            maximum_mass = mass_range[1]

        else:
            minimum_mass = 0
            maximum_mass = 200

        if log_logic == False:
            if mov:
                artist1,=ax.plot([amass_range[0]-.5,amass_range[1]+.5],[1,1],'k-')
                artists.append(artist1)
            else:
                pl.plot([amass_range[0]-.5,amass_range[1]+.5],[1,1],'k-')
        else:
            if mov:
                artist1,=ax.plot([amass_range[0]-.5,amass_range[1]+.5],[0,0],'k-')
                artists.append(artist1)
            else:
                pl.plot([amass_range[0]-.5,amass_range[1]+.5],[0,0],'k-')

        labelsx=[]
        if (maximum_mass-minimum_mass) > 100:
            delta_labelsx = 10
        else:
            delta_labelsx = 5

        iii = amass_range[0]%delta_labelsx
        if iii == 0:
            labelsx.append(str(amass_range[0]))
        else:
            labelsx.append(' ')
        iii = iii+1
        kkk = 0
        for label1 in range(amass_range[1]-amass_range[0]):
            if iii == 5:
                kkk = kkk+1
                labelsx.append(str((iii*kkk)+amass_range[0]-(amass_range[0]%5)))
                iii = 0
                iii = iii+1
            else:
                labelsx.append(' ')
                iii = iii+1

        if delta_labelsx == 5:
            xticks = arange(amass_range[0],amass_range[1],1)
            pl.xticks(xticks,labelsx)
        else:
            pl.xticks()

        # SJONES moved the pl.grid and pl.show to the very end
        if grid:
            pl.grid()
        if show:
            pl.show()

##!!FOR!!###### print 'LEN LABELS= ', len(labelsx)
##DEBUGGING####
####!!!######## for bbb in range (len(labelsx)):
###############     print labelsx[bbb]
        if mov:
            return artists

    def elemental_abund(self,cycle,zrange=[1,15],ylim=[-12,0],title_items=None,
                        ref=0,solar_filename='',show_names=True,label='',
                        colour='',**kwargs):
        '''
        Plot the decayed elemental abundance distribution (PPN).
        Plot the elemental abundance distribution (nugridse).
        (FH, 06/2014; SJ 07/2014)

        Parameters
        ----------
        cycle : string, integer or list
            The cycle of interest.  If it is a list of cycles, this
            method will do a plot for each cycle and save them to a
            file.
        zrange : list, optional
            A 1x2 array containing the lower and upper atomic number
            limit
        ylim : list, optional
            A 1x2 array containing the lower and upper Y limits.  If
            it is [0,0], then ylim will be determined automatically.
            The default is [0,0].
        title_items : list, optional
            A list of cycle attributes that will be added to the title.
            For possible cycle attributes see self.cattrs.
        ref : integer, optional
            ref = 0: plot abundaces as mass fraction
            ref = 1: plot abundances relative to solar
                for this option, an additional keyword argument 'solar_filename'
                must be passed with the path to the solar abundance dist. file.
            ref = 2: in progress
        label : string, optional
            The label for the abundance distribution
            The default is '' (i.e. do not show a label)
        show_names : boolean, optional
            Whether or not to show the element names on the figure.
        colour : string, optional
            In case you want to dictate line colours. Takes cymkrgb single-character colours
            or any other colour string accepted by matplotlib.
            The default is '' (automatic colour selection)
        kwargs : additional keyword arguments
            These arguments are equivalent to those of iso_abund, e.g.
            mass_range. Routines from iso_abund are called, to perform
            averages and get elemental abundances in the correct form.

        Output
        ------
        z_el : array
            proton number of elements being returned
        el_to_plot : array
            elemental abundances (as you asked for them, could be ref to something else)

        '''
        plotType=self._classTest()
        if plotType=='PPN':
            self.get(cycle,decayed=True)
            z_el=unique(self.z_iso_to_plot)
            zmin_ind=min(where(z_el>=zrange[0])[0])
            zmax_ind=max(where(z_el<=zrange[1])[0])
            # extract some elemental quantities:
            a_el=[]; el_name=[]; el_abu=[]; el_abu_hash={}
            for z in z_el[zmin_ind:zmax_ind]:
                el=self.el_iso_to_plot[where(self.z_iso_to_plot==z)[0].tolist()[0]]
                X_el=self.abunds[where(self.el_iso_to_plot==el)[0].tolist()].sum()
                a_el.append(self.a_iso_to_plot[where(self.z_iso_to_plot==z)[0].tolist()[0]])
                el_abu.append(X_el)
                el_name.append(el)
                el_abu_hash[el]=X_el
            # plot an elemental abundance distribution with labels:
            pl.plot(z_el[zmin_ind:zmax_ind],np.log10(el_abu),**kwargs)
            j=0        # add labels
            for z in z_el[zmin_ind:zmax_ind]:
               pl.text(z+0.15,log10(el_abu[j])+0.05,el_name[j])
               j += 1
            if title_items is not None:
               pl.title(self._do_title_string(title_items,cycle))
            pl.ylim(ylim[0],ylim[1])
            pl.xlabel('Z')
            pl.ylabel('log mass fraction')
#           savefig('elemental'+str(i)+'.png')
        elif plotType=='se':
            # get self.***_iso_to_plot by calling iso_abund function, which writes them
            self.iso_abund(cycle,elemaburtn=True,**kwargs)
            z_el=unique(self.se.Z)
            zmin_ind=min(where(z_el>=zrange[0])[0])
            zmax_ind=max(where(z_el<=zrange[1])[0])
            # extract some elemental quantities:
            a_el=[]; el_name=[]; el_abu=[]; el_abu_hash={}
            for z in z_el[zmin_ind:zmax_ind]:
                el=self.el_iso_to_plot[where(self.se.Z==z)[0].tolist()[0]]
                X_el=self.abunds[where(self.el_iso_to_plot==el)[0].tolist()].sum()
                a_el.append(self.a_iso_to_plot[where(self.z_iso_to_plot==z)[0].tolist()[0]])
                el_abu.append(X_el)
                el_name.append(el)
                el_abu_hash[el]=X_el
            # plot an elemental abundance distribution with labels:
            if ref==0:
                el_abu_plot=el_abu
                ylab='log mass fraction'
            elif ref==1:
                from . import utils
                if solar_filename=='':
                    raise IOError('You chose to plot relative to the solar abundance dist. However, you did not supply the solar abundance file!')
                else:
                    nuutils.solar(solar_filename,1)
                    menow = where(unique(nuutils.z_sol)==44.)[0][0]
                    print(1, menow, nuutils.solar_elem_abund[menow])
                    el_abu_sun=np.array(nuutils.solar_elem_abund)
                    print(2, el_abu_sun)
                    print(3, el_abu_sun[42])
                    el_abu_plot=np.zeros(len(el_abu))
                    for zs in z_el[zmin_ind:zmax_ind]:
                        zelidx=where(z_el[zmin_ind:zmax_ind]==zs)[0]
                        zsolidx=zs-1
                        if el_abu_sun[zsolidx] > 0. :
                            el_abu_plot[zelidx]=old_div(el_abu[zelidx],el_abu_sun[zsolidx])
                        else:
                            el_abu_plot[zelidx]=-1

                    ylab='log X/X$_\odot$'
            else:
                raise IOError('Your choice of ref is not available yet. Please use another.')
            if label != '':
                if colour!='':
                    print("Plotting without color and label:")
                    pl.plot(z_el[zmin_ind:zmax_ind],np.log10(el_abu_plot),
                            'o-',label=label,color=colour,markeredgecolor='None')
                else:
                    pl.plot(z_el[zmin_ind:zmax_ind],np.log10(el_abu_plot)
                            ,'o-',label=label,markeredgecolor='None')
            else:
                if colour!='':
                    pl.plot(z_el[zmin_ind:zmax_ind],np.log10(el_abu_plot),
                            'o-',color=colour,markeredgecolor='None')
                else:
                    pl.plot(z_el[zmin_ind:zmax_ind],np.log10(el_abu_plot),
                            'o-',markeredgecolor='None')
            if show_names:
                j=0        # add labels
                for z in z_el[zmin_ind:zmax_ind]:
#                    pl.text(z+0.15,log10(el_abu_plot[j])+0.05,el_name[j])
                    if el_abu_plot[j] > 0.:
                        pl.text(z,log10(el_abu_plot[j])+0.5,el_name[j],
                            horizontalalignment='center')
                    j += 1
            if title_items is not None:
                pl.title(self._do_title_string(title_items,cycle))
            pl.ylim(ylim[0],ylim[1])
            pl.xlabel('Z')
            pl.ylabel(ylab)
            if label != '':
                pl.legend(loc='best').draw_frame(False)
            return z_el[zmin_ind:zmax_ind],el_abu_plot
        else:
            print('This method is not supported for '+plotType)
            return

    def _do_title_string(self,title_items,cycle):
        '''
        Create title string

        Private method that creates a title string for a cycle plot
        out of a list of title_items that are cycle attributes and can
        be obtained with self.get

        Parameters
        ----------
        title_items : list
            A list of cycle attributes.
        cycle : scalar
            The cycle for which the title string should be created.

        Returns
        -------
        title_string: string
            Title string that can be used to decorate plot.

        '''

        title_string=[]
        form_str='%4.1F'

        for item in title_items:
            num=self.get(item,fname=cycle)
            if num > 999:
                num=log10(num)
                prefix='log '
            else:
                prefix=''
            title_string.append(prefix+item+'='+form_str%num)
        tt=''
        for thing in title_string:
            tt = tt+thing+", "
        return tt.rstrip(', ')

    def plotprofMulti(self, ini, end, delta, what_specie, xlim1, xlim2,
                      ylim1, ylim2, symbol=None):

        '''
        create a movie with mass fractions vs mass coordinate between
        xlim1 and xlim2, ylim1 and ylim2. Only works with instances of
        se.

        Parameters
        ----------
        ini : integer
            Initial model i.e. cycle.
        end : integer
            Final model i.e. cycle.
        delta : integer
            Sparsity factor of the frames.
        what_specie : list
            Array with species in the plot.
        xlim1, xlim2 : integer or float
            Mass coordinate range.
        ylim1, ylim2 : integer or float
            Mass fraction coordinate range.
        symbol : list, optional
            Array indicating which symbol you want to use. Must be of
            the same len of what_specie array.  The default is None.

        '''
        plotType=self._classTest()
        if plotType=='se':
            for i in range(ini,end+1,delta):
                step = int(i)
                #print step
                if symbol==None:
                    symbol_dummy = '-'
                    for j in range(len(what_specie)):
                        self.plot_prof_1(step,what_specie[j],xlim1,xlim2,ylim1,ylim2,symbol_dummy)
                else:
                    for j in range(len(what_specie)):
                        symbol_dummy = symbol[j]
                        self.plot_prof_1(step,what_specie[j],xlim1,xlim2,ylim1,ylim2,symbol_dummy)

                #
                filename = str('%03d' % step)+'_test.png'
                pl.savefig(filename, dpi=400)
                print('wrote file ', filename)
                #
                pl.clf()

        else:
            print('This method is not supported for '+str(self.__class__))
            return

    def movie(self, cycles, plotstyle='',movname='',fps=12,**kwargs):
        from matplotlib import animation
        '''
            Make an interactive movie in the matplotlib window for a number of
            different plot types:

            Plot types
            ----------
            'iso_abund' : abundance distribution a la se.iso_abund()
            'abu_chart' : abundance chart a la se.abu_chart()
            'plot'      : plot any number of y_items against an x_item

            Parameters
            ----------
            cycles : list
                Which cycles do you want to plot as movie frames?
            plotstyle : string
                What type of plot should the movie show? Currently supported is
                'iso_abund', 'abu_chart' and 'plot'
            movname : string, optional
                Name of movie (+ extension, e.g. .mp4 or .avi) if the movie is
                to be saved
                The default is ''
            args : list
                Arguments to should be passed to the plotting function. These are
                the arguments of the respective methods that make the frames. See
                the docstrings of those functions for details

            'plot' Parameters
            -----------------
            'xlims'    : tuple, optional
            'ylims'    : tuple, optional
            'xlabel'   : string, optional
            'ylabel'   : string, optional
            'legend'   : boolean, optional
                The default is False
            'loc'      : string or integer, optional
                Set the legend location if legend is True.
                The default is 'best'
            'interval' : frame interval in ms

            FAQ:
            ----
            If ffmpeg is not installed on OSX (and you don't want to wait for port to do it) check out
            these binaries:
            http://stackoverflow.com/questions/18833731/how-to-set-ffmpeg-for-matplotlib-in-mac-os-x

            '''

        modelself=self
        supported_styles=['iso_abund','abu_chart','plot']
        class mov(object):

            def __init__(self,cyc,style,movname,fps,**kwargs):
                self.fig = None
                self.ax = None
                self.ani = None
                self.cyc = cyc
                self.movname=movname
                self.fps = fps
                self.style=style
                if self.style in supported_styles:
                    animateFunc=draw_frame
                else:
                    raise IOError("this type of movie is not available yet! Sorry!")
                if self.style=='plot':
                    self.y_ditems=kwargs['y_items']
                    self.data=kwargs['data']
                self._init_animation(animateFunc)

            def _init_animation(self, animateFunc):
                if self.style=='plot':
                    fsize=14
                    params = {'axes.labelsize':  fsize,
                        'text.fontsize':   fsize,
                        'legend.fontsize': fsize*0.8,
                        'xtick.labelsize': fsize,
                        'ytick.labelsize': fsize,
                        'text.usetex': False,
                        'figure.facecolor': 'white',
                        'ytick.minor.pad': 8,
                        'ytick.major.pad': 8,
                        'xtick.minor.pad': 8,
                        'xtick.major.pad': 8,
                        'figure.subplot.bottom' : 0.15,
                        'lines.markersize': 8}
                    matplotlib.rcParams.update(params)

                    self.fig, self.ax = pl.subplots()
                    tmp=[]
                    for i in range(len(self.y_ditems)):
                        tmp.append(self.data[0][0])
                        tmp.append(self.data[0][i+1])
                    self.lines = self.ax.plot(*tmp)
                    if 'ylims' in kwargs:
                        pl.ylim(kwargs['ylims'])
                    if 'xlims' in kwargs:
                        pl.xlim(kwargs['xlims'])
                    if 'xlabel' in kwargs:
                        pl.xlabel(kwargs['xlabel'])
                    else:
                        pl.xlabel(kwargs['x_item'])
                    if 'ylabel' in kwargs:
                        pl.ylabel(kwargs['ylabel'])
                    else:
                        if type(y_items) is str:
                            lab=y_items
                        elif type(y_items) is list and len(y_items) == 1:
                            lab=y_items[0]
                        else:
                            lab=''
                            for el in kwargs['y_items']:
                                lab+=el+', '
                            lab=lab[:-2]
                        pl.ylabel(lab)
                    if 'legend' in kwargs and kwargs['legend']:
                        if 'loc' in kwargs:
                            pl.legend([line for line in self.lines], self.y_ditems,
                                      loc=kwargs['loc']).draw_frame(False)
                        else:
                            pl.legend([line for line in self.lines], self.y_ditems,
                                      loc='best').draw_frame(False)

                self._animation(animateFunc)

            def _animation(self, animateFunc):
                if plotstyle=='plot' and 'interval' in kwargs:
                    self.ani = animation.FuncAnimation(self.fig, animateFunc, arange(0, len(self.cyc)), interval=kwargs['interval'], blit=False, fargs=[self])
                elif plotstyle=='iso_abund':
                    self.fig, self.ax = pl.subplots()
                    ims=[]
                    for i in arange(0,len(self.cyc)):
                        im=draw_frame(i,self)
                        ims.append(im)
                    self.ani = animation.ArtistAnimation(self.fig,ims,interval=50,
                                                             blit=False)
                    self.fig.canvas.draw()
                elif plotstyle=='abu_chart':
                    self.fig=pl.figure()
                    axx = 0.10
                    axy = 0.10
                    axw = 0.85
                    axh = 0.8
                    self.ax=pl.axes([axx,axy,axw,axh])
                    ims=[]
                    for i in arange(0,len(self.cyc)):
                        im=draw_frame(i,self)
                        # draw_frame here returns the patch for the abundance squares
                        # im[0] as well as the artists im[1], so that the colorbar
                        # can be plotted only once (on the first plot)
                        ims.append(im[1])
                        if i==0:
                            cb=pl.colorbar(im[0])
                            cb.set_label('log$_{10}$(X)',fontsize='x-large')
                    self.ani = animation.ArtistAnimation(self.fig,ims,interval=50,
                                                         blit=False)
                    self.fig.canvas.draw()

                if self.movname is not '':
                    print('\n generating animation: '+self.movname)
                    self.ani.save(self.movname,fps=self.fps)
                    print('animation '+self.movname+' saved with '+str(self.fps)+' frames per second')

        plotType=self._classTest()

        if plotType=='se':
            if plotstyle == 'iso_abund':
                data = self.se.get(cycles,['iso_massf','mass'])
                def draw_frame(i,self=None):
                    artists=modelself.iso_abund(self.cyc[i],stable=True,show=False,
                                             data_provided=True,thedata=data[i],
                                             verbose=False,drawfig=self.fig,drawax=self.ax,
                                             mov=True,**kwargs)
                    return artists
            if plotstyle == 'abu_chart':
                data = self.se.get(cycles,['iso_massf','mass'])
                def draw_frame(i,self=None):
                    artists=modelself.abu_chart(self.cyc[i],show=False,data_provided=True,
                                        thedata=data[i],lbound=(-12, -6),drawfig=self.fig,
                                        drawax=self.ax,mov=True,**kwargs)
                    return artists
            if plotstyle=='plot':
                if 'x_item' not in kwargs or 'y_items' not in kwargs:
                    raise IOError("Please specify both x_item and y_items")
                x_item = kwargs['x_item']
                y_items = kwargs['y_items']
                tx, ty = type(x_item), type(y_items)
                if tx is list and ty is list:
                    data=self.se.get(cycles,x_item+y_items)
                elif tx is str and ty is list:
                    data=self.se.get(cycles,[x_item]+y_items)
                elif tx is str and ty is str:
                    data=self.se.get(cycles,[x_item]+[y_items])
                def draw_frame(i, self=None):
#                    pl.title("cycle: " + self.cyc[i])
                    for j in range(len(self.lines)):
                        if 'logy' in kwargs and kwargs['logy']:
                            self.lines[j].set_data(self.data[i][0],
                                                   np.log10(self.data[i][j+1]))
                        else:
                            self.lines[j].set_data(self.data[i][0],
                                                   self.data[i][j+1])
                    return  self.lines

        if plotstyle=='plot':
            return mov(cycles,plotstyle,movname,fps,data=data,**kwargs).ani
        else:
            return mov(cycles,plotstyle,movname,fps).ani

    # From mesa_profile
    def plot_prof_1(self, species, keystring, xlim1, xlim2, ylim1,
                    ylim2, symbol=None, show=False):
        '''
        Plot one species for cycle between xlim1 and xlim2 Only works
        with instances of se and mesa _profile.

        Parameters
        ----------
        species : list
            Which species to plot.
        keystring : string or integer
            Label that appears in the plot or in the case of se, a
            cycle.
        xlim1, xlim2 : integer or float
            Mass coordinate range.
        ylim1, ylim2 : integer or float
            Mass fraction coordinate range.
        symbol : string, optional
            Which symbol you want to use.  If None symbol is set to '-'.
            The default is None.
        show : boolean, optional
            Show the ploted graph.  The default is False.

        '''
        plotType=self._classTest()
        if plotType=='se':
            #tot_mass=self.se.get(keystring,'total_mass')
            tot_mass=self.se.get('mini')
            age=self.se.get(keystring,'age')
            mass=self.se.get(keystring,'mass')
            Xspecies=self.se.get(keystring,'iso_massf',species)

            mod=keystring
        elif plotType=='mesa_profile':
            tot_mass=self.header_attr['star_mass']
            age=self.header_attr['star_age']
            mass=self.get('mass')
            mod=self.header_attr['model_number']
            Xspecies=self.get(species)
        else:
            print('This method is not supported for '+str(self.__class__))
            return

        if symbol == None:
            symbol = '-'

        x,y=self._logarithm(Xspecies,mass,True,False,10)
        #print x
        pl.plot(y,x,symbol,label=str(species))
        pl.xlim(xlim1,xlim2)
        pl.ylim(ylim1,ylim2)
        pl.legend()

        pl.xlabel('$Mass$ $coordinate$', fontsize=20)
        pl.ylabel('$X_{i}$', fontsize=20)
        #pl.title('Mass='+str(tot_mass)+', Time='+str(age)+' years, cycle='+str(mod))
        pl.title('Mass='+str(tot_mass)+', cycle='+str(mod))
        if show:
            pl.show()

    def density_profile(self,ixaxis='mass',ifig=None,colour=None,label=None,fname=None):
        '''
        Plot density as a function of either mass coordiate or radius.

        Parameters
        ----------
        ixaxis : string
            'mass' or 'radius'
            The default value is 'mass'
        ifig : integer or string
            The figure label
            The default value is None
        colour : string
            What colour the line should be
            The default value is None
        label : string
            Label for the line
            The default value is None
        fname : integer
            What cycle to plot from (if SE output)
            The default value is None
        '''

        pT=self._classTest()

        # Class-specific things:
        if pT is 'mesa_profile':
            x = self.get(ixaxis)
            if ixaxis is 'radius':
                x = x*ast.rsun_cm
            y = self.get('logRho')
        elif pT is 'se':
            if fname is None:
                raise IOError("Please provide the cycle number fname")
            x = self.se.get(fname,ixaxis)
            y = np.log10(self.se.get(fname,'rho'))
        else:
            raise IOError("Sorry. the density_profile method is not available \
                          for this class")

        # Plot-specific things:
        if ixaxis is 'radius':
            x = np.log10(x)
            xlab='$\log_{10}(r\,/\,{\\rm cm})$'
        else:
            xlab='${\\rm Mass}\,/\,M_\odot$'

        if ifig is not None:
            pl.figure(ifig)
        if label is not None:
            if colour is not None:
                pl.plot(x,y,color=colour,label=label)
            else:
                pl.plot(x,y,label=label)
            pl.legend(loc='best').draw_frame(False)
        else:
            if colour is not None:
                pl.plot(x,y,color=colour)
            else:
                pl.plot(x,y)

        pl.xlabel(xlab)
        pl.ylabel('$\log_{10}(\\rho\,/\,{\\rm g\,cm}^{-3})$')

    def abu_profile(self,ixaxis='mass',isos=None,ifig=None,fname=None,logy=False,
                    colourblind=False):
        '''
        Plot common abundances as a function of either mass coordiate or radius.

        Parameters
        ----------
        ixaxis : string, optional
            'mass',  'logradius' or 'radius'
            The default value is 'mass'
        isos : list, optional
            list of isos to plot, i.e. ['h1','he4','c12'] for MESA or
            ['H-1','He-4','C-12'] for SE output. If None, the code decides
            itself what to plot.
            The default is None.
        ifig : integer or string, optional
            The figure label
            The default value is None
        fname : integer, optional
            What cycle to plot from (if SE output)
            The default value is None
        logy : boolean, optional
            Should y-axis be logarithmic?
            The default value is False
        colourblind : boolean, optional
            do you want to use the colourblind colour palette from the NuGrid
            nuutils module?
        '''

        pT=self._classTest()
        # Class-specific things:
        if pT is 'mesa_profile':
            x = self.get(ixaxis)
            if ixaxis is 'radius':
                x = x*ast.rsun_cm
            if isos is None:
                isos=['h1','he4','c12','c13','n14','o16','ne20','ne22','mg24','mg25',
                      'al26','si28','si30','s32','s34','cl35','ar36','ar38','cr52',
                      'cr56','fe56','ni56']
            risos=[i for i in isos if i in self.cols]
            abunds = [self.get(riso) for riso in risos]
            names=risos
        elif pT is 'se':
            if fname is None:
                raise IOError("Please provide the cycle number fname")
            x = self.se.get(fname,ixaxis)
            if isos is None:
                isos=['H-1','He-4','C-12','C-13','N-14','O-16','Ne-20','Ne-22','Mg-24','Mg-25',
                      'Sl-26','Si-28','Si-30','S-32','S-34','Cl-35','Ar-36','Ar-38','Cr-52',
                      'Cr-56','Fe-56','Ni-56']
            risos=[i for i in isos if i in self.se.isotopes]
            abunds = self.se.get(fname,'iso_massf',risos)
            names=risos
        else:
            raise IOError("Sorry. the density_profile method is not available \
                          for this class")

        # Plot-specific things:
        if ixaxis is 'logradius':
            x = np.log10(x)
            xlab='$\log_{10}(r\,/\,{\\rm cm})$'
        elif ixaxis is 'radius':
            x = old_div(x, 1.e8)
            xlab = 'r / Mm'
        else:
            xlab='${\\rm Mass}\,/\,M_\odot$'

        if ifig is not None:
            pl.figure(ifig)
        from . import utils as u
        cb = u.colourblind
        lscb = u.linestylecb # colourblind linestyle function
        for i in range(len(risos)):
            if logy:
                y = np.log10(abunds if len(risos) < 2 else abunds[i])
            else:
                y = abunds if len(risos) < 2 else abunds[i]
            if colourblind:
                pl.plot(x,y,ls=lscb(i)[0],marker=lscb(i)[1],
                        color=lscb(i)[2],markevery=u.linestyle(i)[1]*20,
                        label=names[i],mec='None')
            else:
                pl.plot(x,y,u.linestyle(i)[0],markevery=u.linestyle(i)[1]*20,
                        label=names[i],mec='None')

        pl.legend(loc='best').draw_frame(False)
        pl.xlabel(xlab)
        pl.ylabel('$\log(X)$')

    def ye_profile(self,ixaxis='mass',ifig=None,colour=None,label=None,fname=None):
        '''
            Plot electron fraction Y_e as a function of either mass coordiate or radius.

            Parameters
            ----------
            ixaxis : string
            'mass' or 'radius'
            The default value is 'mass'
            ifig : integer or string
            The figure label
            The default value is None
            colour : string
            What colour the line should be
            The default value is None
            label : string
            Label for the line
            The default value is None
            fname : integer
            What cycle to plot from (if SE output)
            The default value is None
            '''

        pT=self._classTest()

        # Class-specific things:
        if pT is 'mesa_profile':
            x = self.get(ixaxis)
            if ixaxis is 'radius':
                x = x*ast.rsun_cm
            y = self.get('ye')
        elif pt is 'se':
            if fname is None:
                raise IOError("Please provide the cycle number fname")
            raise IOError("Sorry, Ye profile is not yet available for\
                          the SE data class")
        #            x = self.se.get(fname,ixaxis)
        #            y = None # not implemented yet...
        else:
            raise IOError("Sorry. the density_profile method is not available \
                          for this class")

        # Plot-specific things:
        if ixaxis is 'radius':
            x = np.log10(x)
            xlab='$\log_{10}(r\,/\,{\\rm cm})$'
        else:
            xlab='${\\rm Mass}\,/\,M_\odot$'

        if ifig is not None:
            pl.figure(ifig)
        if label is not None:
            if colour is not None:
                pl.plot(x,y,color=colour,label=label)
            else:
                pl.plot(x,y,label=label)
            pl.legend(loc='best').draw_frame(False)
        else:
            if colour is not None:
                pl.plot(x,y,color=colour)
            else:
                pl.plot(x,y)

        pl.xlabel(xlab)
        pl.ylabel('$Y_{\\rm e}$')


    # From mesa.star_log

def flux_chart(file_name, plotaxis, plot_type, which_flux=None,
               I_am_the_target=None, prange=None):
    '''
    Plots a chart with fluxes

    Parameters
    ----------
    file_name : string
        Name of the file of fluxes we are looking at.
    plotaxis : list
        [xmin, xmax, ymin, ymax], where on x axis there is neutron
        number and on y axis there is Z.
    plot_types : integer
        Set to 0 for standard flux plot.  Set to 1 if fluxes focused
        on one specie.
    which_flux : integer, optional
        Set to 0 for nucleosynthesis flux plot.  Set to 1 is for energy
        flux plot.  Seting to None is the same a 0.  The default is
        None.
    I_am_the_target : list, optional
        A 2xArray used only if plot_type=1, and is given by [neutron
        number, proton number].  The default is None.
    prange : integer, optional
        The range of fluxes to be considered.  If prange is None, then
        8 fluxes are consdered.  The default is None.

    Notes
    -----
    This script is terribly slow and needs to be improved.  For now I
    put here in data_plot:

    [1]: import data_plot

    [2]: data_plot.flux_chart('file_name', [xmin, xmax, ymin, ymax],
                              int, which_flux, I_am_the_target, prange)

    The pdf is created, but an error bumped up and the gui is empty.
    To avoid this, I had to set 'text.usetex': False.  See below.  Also,
    for the same reason no label in x axys is written using
    'text.usetex': True.

    Note also that the GUI works really slow with this plot.  So, we
    need to optimize from the graphic point of view.  This need to be
    included in ppn.py I think, and set in multi option too, in case
    we want to read more flux files at the same time.

    Finally, you need to have stable.dat to read in to make it work ...

    '''

    import numpy as np
    import matplotlib.pyplot as plt
    #from matplotlib.mpl import colors,cm # deppreciated in mpl ver 1.3
                                          # use line below instead
    from matplotlib import colors,cm
    from matplotlib.patches import Rectangle, Arrow
    from matplotlib.collections import PatchCollection
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea
    import sys

    print_max_flux_in_plot =  True


    f = open(file_name)
    lines = f.readline()
    lines = f.readlines()
    f.close()

    # starting point of arrow
    coord_x_1 = []
    coord_y_1 = []
    # ending point of arrow (option 1)
    coord_x_2 = []
    coord_y_2 = []
    # ending point of arrow (option 2)
    coord_x_3 = []
    coord_y_3 = []
    # fluxes
    flux_read = []
    flux_log10 = []

    if which_flux == None or which_flux == 0:
        print('chart for nucleosynthesis fluxes [dYi/dt]')
        line_to_read = 9
    elif which_flux == 1:
        print('chart for energy fluxes')
        line_to_read = 10
    elif which_flux == 2:
        print('chart for timescales')
        line_to_read = 11
    elif which_flux > 2:
        print("you have only option 0, 1 or 2, not larger than 2")

    single_line = []
    for i in range(len(lines)):
        single_line.append(lines[i].split())
        coord_y_1.append(float(single_line[i][1]))
        coord_x_1.append(float(single_line[i][2])-coord_y_1[i])
        coord_y_2.append(float(single_line[i][5]))
        coord_x_2.append(float(single_line[i][6])-coord_y_2[i])
        coord_y_3.append(float(single_line[i][7]))
        coord_x_3.append(float(single_line[i][8])-coord_y_3[i])
        try:
            flux_read.append(float(single_line[i][line_to_read]))
        except ValueError: # this is done to avoid format issues like 3.13725-181...
            flux_read.append(1.0E-99)
        flux_log10.append(np.log10(flux_read[i]+1.0e-99))

    print('file read!')

    # I need to select smaller sample, with only fluxes inside plotaxis.
    coord_y_1_small=[]
    coord_x_1_small=[]
    coord_y_2_small=[]
    coord_x_2_small=[]
    coord_y_3_small=[]
    coord_x_3_small=[]
    flux_log10_small = []
    for i in range(len(flux_log10)):
        I_am_in = 0
        if coord_y_1[i] > plotaxis[2] and coord_y_1[i] < plotaxis[3] and coord_x_1[i] > plotaxis[0] and coord_x_1[i] < plotaxis[1]:
            I_am_in = 1
            coord_y_1_small.append(coord_y_1[i])
            coord_x_1_small.append(coord_x_1[i])
            coord_y_2_small.append(coord_y_2[i])
            coord_x_2_small.append(coord_x_2[i])
            coord_y_3_small.append(coord_y_3[i])
            coord_x_3_small.append(coord_x_3[i])
            flux_log10_small.append(flux_log10[i])
        if coord_y_3[i] > plotaxis[2] and coord_y_3[i] < plotaxis[3] and coord_x_3[i] > plotaxis[0] and coord_x_3[i] < plotaxis[1] and I_am_in == 0:
            I_am_in = 1
            coord_y_1_small.append(coord_y_1[i])
            coord_x_1_small.append(coord_x_1[i])
            coord_y_2_small.append(coord_y_2[i])
            coord_x_2_small.append(coord_x_2[i])
            coord_y_3_small.append(coord_y_3[i])
            coord_x_3_small.append(coord_x_3[i])
            flux_log10_small.append(flux_log10[i])


    if plot_type == 1:
        print('I_am_the_target=',I_am_the_target)
        #I_am_the_target = [56.-26.,26.]
    # here below need for plotting
    # plotaxis = [xmin,xmax,ymin,ymax]
    #plotaxis=[1,20,1,20]
    #plotaxis=[0,0,0,0]

    # elemental labels off/on [0/1]
    ilabel = 1

    # label for isotopic masses off/on [0/1]
    imlabel = 1

    # turn lines for magic numbers off/on [0/1]
    imagic = 0

    # flow is plotted over "prange" dex. If flow < maxflow-prange it is not plotted
    if prange == None:
        print('plot range given by default')
        prange = 8.

    #############################################

    # we should scale prange on plot_axis range, not on max_flux!
    max_flux = max(flux_log10)
    ind_max_flux = flux_log10.index(max_flux)
    max_flux_small = max(flux_log10_small)
    min_flux = min(flux_log10)
    ind_min_flux = flux_log10.index(min_flux)
    min_flux_small = min(flux_log10_small)

    #nzmax = int(max(max(coord_y_1),max(coord_y_2),max(coord_y_3)))+1
    #nnmax = int(max(max(coord_x_1),max(coord_x_2),max(coord_x_3)))+1
    nzmax = int(max(max(coord_y_1_small),max(coord_y_2_small),max(coord_y_3_small)))+1
    nnmax = int(max(max(coord_x_1_small),max(coord_x_2_small),max(coord_x_3_small)))+1


    nzycheck = np.zeros([nnmax,nzmax,3])
    #coord_x_out = np.zeros(len(coord_x_2))
    #coord_y_out = np.zeros(len(coord_y_2))
    #for i in range(len(flux_log10)):
    #       nzycheck[coord_x_1[i],coord_y_1[i],0] = 1
    #       nzycheck[coord_x_1[i],coord_y_1[i],1] = flux_log10[i]
    #       if coord_x_2[i] >= coord_x_3[i]:
    #                coord_x_out[i] = coord_x_2[i]
    #                coord_y_out[i] = coord_y_2[i]
    #               nzycheck[coord_x_out[i],coord_y_out[i],0] = 1
    #               nzycheck[coord_x_out[i],coord_y_out[i],1] = flux_log10[i]
    #        elif coord_x_2[i] < coord_x_3[i]:
    #                coord_x_out[i] = coord_x_3[i]
    #                coord_y_out[i] = coord_y_3[i]
    #               nzycheck[coord_x_out[i],coord_y_out[i],0] = 1
    #               nzycheck[coord_x_out[i],coord_y_out[i],1] = flux_log10[i]
    #       if flux_log10[i]>max_flux-prange:
    #               nzycheck[coord_x_1[i],coord_y_1[i],2] = 1
    #               nzycheck[coord_x_out[i],coord_y_out[i],2] = 1
    coord_x_out = np.zeros(len(coord_x_2_small))
    coord_y_out = np.zeros(len(coord_y_2_small))
    for i in range(len(flux_log10_small)):
        nzycheck[coord_x_1_small[i],coord_y_1_small[i],0] = 1
        nzycheck[coord_x_1_small[i],coord_y_1_small[i],1] = flux_log10_small[i]
        if coord_x_2_small[i] >= coord_x_3_small[i]:
            coord_x_out[i] = coord_x_2_small[i]
            coord_y_out[i] = coord_y_2_small[i]
            nzycheck[coord_x_out[i],coord_y_out[i],0] = 1
            nzycheck[coord_x_out[i],coord_y_out[i],1] = flux_log10_small[i]
        elif coord_x_2_small[i] < coord_x_3_small[i]:
            coord_x_out[i] = coord_x_3_small[i]
            coord_y_out[i] = coord_y_3_small[i]
            nzycheck[coord_x_out[i],coord_y_out[i],0] = 1
            nzycheck[coord_x_out[i],coord_y_out[i],1] = flux_log10_small[i]
        if which_flux == None or which_flux < 2 and flux_log10_small[i]>max_flux_small-prange:
            nzycheck[coord_x_1_small[i],coord_y_1_small[i],2] = 1
            nzycheck[coord_x_out[i],coord_y_out[i],2] = 1
        elif which_flux == 2 and flux_log10_small[i]<min_flux_small+prange:
            nzycheck[coord_x_1_small[i],coord_y_1_small[i],2] = 1
            nzycheck[coord_x_out[i],coord_y_out[i],2] = 1

    #######################################################################
    # elemental names: elname(i) is the name of element with Z=i
    elname= ('none','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe',
    'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb',
    'Te', 'I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os',
    'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu')


    #### create plot

    ## define axis and plot style (colormap, size, fontsize etc.)
    if plotaxis==[0,0,0,0]:
        xdim=10
        ydim=6
    else:
        dx = plotaxis[1]-plotaxis[0]
        dy = plotaxis[3]-plotaxis[2]
        ydim = 6
        xdim = ydim*dx/dy

    format = 'pdf'
    # note that I had to set 'text.usetex': False, to avoid Exception in Tkinter callback.
    # and to make the GUI work properly. Why? some missing package?
    params = {'axes.labelsize':  15,
      'text.fontsize':   15,
      'legend.fontsize': 15,
      'xtick.labelsize': 15,
      'ytick.labelsize': 15,
      'text.usetex': False}
    plt.rcParams.update(params)
    fig=plt.figure(figsize=(xdim,ydim),dpi=100)
    axx = 0.10
    axy = 0.10
    axw = 0.85
    axh = 0.8
    ax=plt.axes([axx,axy,axw,axh])

    # color map choice for abundances
    cmapa = cm.jet
    # color map choice for arrows
    if which_flux == None or which_flux < 2:
        cmapr = cm.autumn
    elif  which_flux == 2:
        cmapr = cm.autumn_r
    # if a value is below the lower limit its set to white
    cmapa.set_under(color='w')
    cmapr.set_under(color='w')
    # set value range for abundance colors (log10(Y))
    norma = colors.Normalize(vmin=-20,vmax=0)
    # set x- and y-axis scale aspect ratio to 1
    ax.set_aspect('equal')
    #print time,temp and density on top
    #temp = '%8.3e' %ff['temp']
    #time = '%8.3e' %ff['time']
    #dens = '%8.3e' %ff['dens']

    #box1 = TextArea("t : " + time + " s~~/~~T$_{9}$ : " + temp + "~~/~~$\\rho_{b}$ : " \
    #      + dens + ' g/cm$^{3}$', textprops=dict(color="k"))
    #anchored_box = AnchoredOffsetbox(loc=3,
    #        child=box1, pad=0.,
    #        frameon=False,
    #        bbox_to_anchor=(0., 1.02),
    #        bbox_transform=ax.transAxes,
    #        borderpad=0.,
    #        )
    #ax.add_artist(anchored_box)

    # Add black frames for stable isotopes
    f = open('stable.dat')

    head = f.readline()
    stable = []

    for line in f.readlines():
        tmp = line.split()
        zz = int(tmp[2])
        nn = int(tmp[3])
        xy = nn-0.5,zz-0.5
        rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=3.)
        rect.set_zorder(2)
        ax.add_patch(rect)

    apatches = []
    acolor = []
    m = old_div(0.8,prange)#0.8/prange
    if which_flux == None or which_flux < 2:
        vmax=np.ceil(max(flux_log10_small))
        vmin=max(flux_log10_small)-prange
        b=-vmin*m+0.1
    elif which_flux == 2:
        vmin=min(flux_log10_small)
        vmax=np.ceil(min(flux_log10_small)+prange)
        b=vmax*m+0.1
    if which_flux == None or which_flux < 3:
        normr = colors.Normalize(vmin=vmin,vmax=vmax)
        print('vmin and vmax =',vmin,vmax)
    ymax=0.
    xmax=0.

    for i in range(len(flux_log10_small)):
        x = coord_x_1_small[i]
        y = coord_y_1_small[i]
        dx = coord_x_out[i]-coord_x_1_small[i]
        dy = coord_y_out[i]-coord_y_1_small[i]
        if plot_type == 0:
            if which_flux == None or which_flux < 2:
                if flux_log10_small[i]>=vmin:
                    arrowwidth = flux_log10_small[i]*m+b
                    arrow = Arrow(x,y,dx,dy, width=arrowwidth)
                    if xmax<x:
                        xmax=x
                    if ymax<y:
                        ymax=y
                    acol = flux_log10_small[i]
                    apatches.append(arrow)
                    acolor.append(acol)
            elif which_flux == 2:
                if flux_log10_small[i]<=vmax:
                    arrowwidth = -flux_log10_small[i]*m+b
                    arrow = Arrow(x,y,dx,dy, width=arrowwidth)
                    if xmax<x:
                        xmax=x
                    if ymax<y:
                        ymax=y
                    acol = flux_log10_small[i]
                    apatches.append(arrow)
                    acolor.append(acol)
        elif plot_type == 1 and which_flux != 2:
            if x==I_am_the_target[0] and y==I_am_the_target[1] and flux_log10_small[i]>=vmin:
                arrowwidth = flux_log10_small[i]*m+b
                arrow = Arrow(x,y,dx,dy, width=arrowwidth)
                if xmax<x:
                    xmax=x
                if ymax<y:
                    ymax=y
                acol = flux_log10_small[i]
                apatches.append(arrow)
                acolor.append(acol)
            if x+dx==I_am_the_target[0] and y+dy==I_am_the_target[1] and flux_log10_small[i]>=vmin:
                arrowwidth = flux_log10_small[i]*m+b
                arrow = Arrow(x,y,dx,dy, width=arrowwidth)
                if xmax<x:
                    xmax=x
                if ymax<y:
                    ymax=y
                acol = flux_log10_small[i]
                apatches.append(arrow)
                acolor.append(acol)
        elif plot_type == 1 and which_flux == 2:
            if x==I_am_the_target[0] and y==I_am_the_target[1] and flux_log10_small[i]<=vmax:
                arrowwidth = -flux_log10_small[i]*m+b
                arrow = Arrow(x,y,dx,dy, width=arrowwidth)
                if xmax<x:
                    xmax=x
                if ymax<y:
                    ymax=y
                acol = flux_log10_small[i]
                apatches.append(arrow)
                acolor.append(acol)
            if x+dx==I_am_the_target[0] and y+dy==I_am_the_target[1] and flux_log10_small[i]<=vmax:
                arrowwidth = -flux_log10_small[i]*m+b
                arrow = Arrow(x,y,dx,dy, width=arrowwidth)
                if xmax<x:
                    xmax=x
                if ymax<y:
                    ymax=y
                acol = flux_log10_small[i]
                apatches.append(arrow)
                acolor.append(acol)

    #apatches = []
    #acolor = []
    #m = 0.8/prange
    #vmax=np.ceil(max(flux_log10))
    #vmin=max(flux_log10)-prange
    #b=-vmin*m+0.1
    #normr = colors.Normalize(vmin=vmin,vmax=vmax)
    #ymax=0.
    #xmax=0.

    #for i in range(len(flux_log10)):
    #       x = coord_x_1[i]
    #       y = coord_y_1[i]
    #       dx = coord_x_out[i]-coord_x_1[i]
    #       dy = coord_y_out[i]-coord_y_1[i]
    #       if plot_type == 0:
    #               if flux_log10[i]>=vmin:
    #                       arrowwidth = flux_log10[i]*m+b
    #                       arrow = Arrow(x,y,dx,dy, width=arrowwidth)
    #                       if xmax<x:
    #                               xmax=x
    #                       if ymax<y:
    #                               ymax=y
    #                       acol = flux_log10[i]
    #                       apatches.append(arrow)
    #                       acolor.append(acol)
    #       elif plot_type == 1:
    #               if x==I_am_the_target[0] and y==I_am_the_target[1] and flux_log10[i]>=vmin:
    #                       arrowwidth = flux_log10[i]*m+b
    #                       arrow = Arrow(x,y,dx,dy, width=arrowwidth)
    #                       if xmax<x:
    #                               xmax=x
    #                       if ymax<y:
    #                               ymax=y
    #                       acol = flux_log10[i]
    #                       apatches.append(arrow)
    #                       acolor.append(acol)
    #               if x+dx==I_am_the_target[0] and y+dy==I_am_the_target[1] and flux_log10[i]>=vmin:
    #                       arrowwidth = flux_log10[i]*m+b
    #                       arrow = Arrow(x,y,dx,dy, width=arrowwidth)
    #                       if xmax<x:
    #                               xmax=x
    #                       if ymax<y:
    #                               ymax=y
    #                       acol = flux_log10[i]
    #                       apatches.append(arrow)
    #                       acolor.append(acol)
    #


        xy = x-0.5,y-0.5
        rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=1.)
        rect.set_zorder(2)
        ax.add_patch(rect)
        xy = x+dx-0.5,y+dy-0.5
        rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=1.)
        rect.set_zorder(2)
        ax.add_patch(rect)


    a = PatchCollection(apatches, cmap=cmapr, norm=normr)
    a.set_array(np.array(acolor))
    a.set_zorder(3)
    ax.add_collection(a)
    cb = plt.colorbar(a)

    # colorbar label
    if which_flux == None or which_flux == 0:
        cb.set_label('log$_{10}$(f)')
    elif which_flux ==1:
        cb.set_label('log$_{10}$(E)')
    elif which_flux ==2:
        cb.set_label('log$_{10}$(timescale)')


    # plot file name
    graphname = 'flow-chart.'+format

    # decide which array to take for label positions
    iarr = 2

    # plot element labels
    for z in range(nzmax):
        try:
            nmin = min(np.argwhere(nzycheck[:,z,iarr-2]))[0]-1
            ax.text(nmin,z,elname[z],horizontalalignment='center',verticalalignment='center',fontsize='medium',clip_on=True)
        except ValueError:
            continue

    # plot mass numbers
    if imlabel==1:
        for z in range(nzmax):
            for n in range(nnmax):
                a = z+n
                if nzycheck[n,z,iarr-2]==1:
                    ax.text(n,z,a,horizontalalignment='center',verticalalignment='center',fontsize='small',clip_on=True)

    # plot lines at magic numbers
    if imagic==1:
        ixymagic=[2, 8, 20, 28, 50, 82, 126]
        nmagic = len(ixymagic)
        for magic in ixymagic:
            if magic<=nzmax:
                try:
                    xnmin = min(np.argwhere(nzycheck[:,magic,iarr-2]))[0]
                    xnmax = max(np.argwhere(nzycheck[:,magic,iarr-2]))[0]
                    line = ax.plot([xnmin,xnmax],[magic,magic],lw=3.,color='r',ls='-')
                except ValueError:
                    dummy=0
            if magic<=nnmax:
                try:
                    yzmin = min(np.argwhere(nzycheck[magic,:,iarr-2]))[0]
                    yzmax = max(np.argwhere(nzycheck[magic,:,iarr-2]))[0]
                    line = ax.plot([magic,magic],[yzmin,yzmax],lw=3.,color='r',ls='-')
                except ValueError:
                    dummy=0

    # set axis limits
    if plotaxis==[0,0,0,0]:
        ax.axis([-0.5,xmax+0.5,-0.5,ymax+0.5])
    else:
        ax.axis(plotaxis)

    # set x- and y-axis label
    ax.set_xlabel('neutron number')
    ax.set_ylabel('proton number')
    if which_flux == None or which_flux == 0:
        max_flux_label="max flux = "+str('{0:.4f}'.format(max_flux))
    elif which_flux == 1:
        max_flux_label="max energy flux = "+str('{0:.4f}'.format(max_flux))
    elif which_flux == 2:
        min_flux_label="min timescale [s] = "+str('{0:.4f}'.format(min_flux))
    if print_max_flux_in_plot:
        if which_flux == None or which_flux < 2:
            ax.text(plotaxis[1]-1.8,plotaxis[2]+0.1,max_flux_label,fontsize=10.)
        elif which_flux == 2:
            ax.text(plotaxis[1]-1.8,plotaxis[2]+0.1,min_flux_label,fontsize=10.)



    fig.savefig(graphname)
    print(graphname,'is done')
    if which_flux == None or which_flux < 2:
        print(max_flux_label,'for reaction =',ind_max_flux+1)
    elif which_flux == 2:
        print(min_flux_label,'for reaction =',ind_min_flux+1)

    plt.show()
