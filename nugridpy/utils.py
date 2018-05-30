'''
utils.py

Utility class for holding extra methods from mesa.py, nuh5p.py

'''
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from builtins import zip
from builtins import str
from builtins import range
from builtins import list
from past.utils import old_div
from builtins import object

import numpy as np
import scipy as sc
from . import ascii_table as att

from scipy import optimize
import matplotlib.pyplot as pl
import os
from decimal import *
class data_fitting(object):
    '''
    Wrapper for the scipy method optimize.leastsq

    Parameters
    ----------
    func : string or function, optional
        If func is a string, then it must be one of the two strings
        'linear' or 'powerlaw'. If func is not a string it must be
        a custom function.  The default is "linear".
    coef : tuple or list
        A guess for the list of coeffiecients.  For 'powerlaw' coef
        must have three enries.  For 'linear' coef must have two
        enries.  If you provide your own function, then provide as
        many coef entries as your function needs.

    Examples
    --------
    Typically you have some data y(x) and you want to
    fit some function to this data.

    in order to test this create some sample data y(x)

    >>> x=arange(0,100)
    >>> y=(((random_sample(100)+0.5)*50.)-50)+x # this is the data to fit

    This class provides two simple fit functions: linear and power law.

    Feel free to add more, or provide your custom function as an
    argument.  For example, in order to provide an exponential fit
    function, first define the function.

    >>> def ff(coef,x):
    >>>     return coef[0]*sc.exp(coef[1]*x)+coef[2]

    and then call the data_fitting instance.

    >>> f=utils.data_fitting(ff,coef=(1,0,0))

    Once you have initialized this class, the instance provides two
    methods; fit and plot to check the fit.

    >>> f.fit(x,y)
    >>> f.plot()
    >>> plot(f.x,(y-f.func(f.fcoef,f.x))/(0.5*(y[-1]-y[0])),'o')
    The fit coefficients are stored in self.fcoef.

    (FH)

    '''

    def __init__(self, func='linear', coef=(1, 1)):
        '''
        Parameters
        ----------
        func : string or function, optional
            If func is a string, then it must be one of the two strings
            'linear' or 'powerlaw'. If func is not a string it must be
            a custom function.  The default is "linear".
        coef : list, optional
            A guess for the list of coeffiecients.  For 'powerlaw' coef
            must have three enries.  For 'linear' coef must have two
            enries.  If you provide your own function, then provide as
            many coef entries as your function needs.  The default is
            (1, 1).

        '''
        if func is 'linear':
            print("Information: 'linear' fit needs coef list with 2 entries")
            print(" -> will use default: coef = "+str(coef))
            if len(coef) is not 2:
                print("Warning: you want a linear fit but you have not")
                print("         provided a guess for coef with the")
                print("         right length (2).")
                print(" -> I will continue and assume coef=(1,1)")
                coef = (1,1)
            def ff(coef,x):
                return coef[0]*x + coef[1]
            self.__name__ = func
        elif func is 'powerlaw':
            print("Information: 'powerlaw' fit needs coef list with 3 entries")
            print(" -> will use default: coef = "+str(coef))
            if len(coef) is not 3:
                print("Warning: you want a power law fit but you have")
                print("         not provided a guess for coef with the")
                print("         right length (3).")
                print(" -> I will continue and assume coef=(1,1,1)")
                coef = (1,1,1)
            def ff(coef,x):
                return coef[0]*x**coef[1] + coef[2]
            self.__name__ = func
        else:
            print("Information: You provide a fit function yourself. I trust")
            print("             you have provided a matching guess for the ")
            print("             coefficient list!")
            ff = func
            self.__name__ = func.__name__


        # we want to determine the coefficients that
        # fit the power law to the data
        # this is done by finding the minimum to a
        # residual function:
        # func(params) = ydata - f(xdata, params)
        # therefore we define a residual function
        def fres(coef,y,ff,x):
            return y-ff(coef,x)

        self.residual = fres
        self.coef     = coef
        self.func     = ff

    def fit(self, x, y, dcoef='none'):
        '''
        performs the fit

        x, y : list
            Matching data arrays that define a numerical function y(x),
            this is the data to be fitted.
        dcoef : list or string
            You can provide a different guess for the coefficients, or
            provide the string 'none' to use the inital guess.  The
            default is 'none'.

        Returns
        -------
        ierr
            Values between 1 and 4 signal success.

        Notes
        -----
        self.fcoef, contains the fitted coefficients.

        '''
        self.x = x
        self.y = y

        if dcoef is not 'none':
            coef = dcoef
        else:
            coef = self.coef

        fcoef=optimize.leastsq(self.residual,coef,args=(y,self.func,x))
        self.fcoef = fcoef[0].tolist()
        return fcoef[1]

    def plot(self, ifig=1, data_label='data', fit_label='fit',
             data_shape='o', fit_shape='-'):
        '''
        plot the data and the fitted function.

        Parameters
        ----------
        ifig : integer
            Figure window number.  The default is 1.
        data_label : string
            Legend for data.  The default is 'data'.
        fit_label : string
            Legend for fit.  If fit_lable is 'fit', then substitute fit
            function type self.func_name.  The default is 'fit'.
        data_shape : character
            Shape for data.  The default is 'o'.
        fit_shape : character
            Shape for fit.  The default is '-'.

        '''

        if len(self.coef) is not len(self.fcoef):
            print("Warning: the fitted coefficient list is not same")
            print("         length as guessed list - still I will try ...")

        pl.figure(ifig)
        pl.plot(self.x,self.y,data_shape,label=data_label)
        if fit_label is 'fit':
            fit_label=self.__name__
        pl.plot(self.x,self.func(self.fcoef,self.x),fit_shape,label=fit_label)
        pl.legend()

class constants(object):
    mass_sun=1.9891e+33
    mass_sun_unit='g'
    one_year=31558149.984
    avogadro=6.02214179e23
    avogadro_unit='mol^-1'

class Utils(object):
    '''
    This private class contains utilities that are used by methods,
    mostly in the ppn and mppnp classes.  Users whould normally not use
    these methods directly.  Things go here when it can be imagined that
    they may be used not in immediate conjunction with plotting.
    Otherwise they would go into the superclass data_plot.

    '''

    #elements_names is marked for deletion (FH, Oct2011) and
    #should be replaced with self.stable_names from self._stable_names()
    elements_names = ['Neutron','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al',
        'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',
        'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc',
        'Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te', 'I','Xe','Cs','Ba','La','Ce',
        'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',
        'W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
        'Pa','U','Np','Pu','Am','Cm','Bk','Cf']

    stable_el = [['Neutron','999'],['H',1, 2],['He', 3, 4],['Li', 6, 7],['Be', 9],
                 ['B', 10, 11],['C', 12, 13],['N', 14, 15],['O', 16, 17, 18],['F', 19],
                 ['Ne', 20, 21, 22],['Na', 23],['Mg', 24, 25, 26],['Al', 27],['Si', 28, 29, 30],
                 ['P', 31],['S', 32, 33, 34, 36],['Cl', 35, 37],['Ar', 36, 38, 40],['K', 39, 40, 41],
                 ['Ca', 40, 42, 43, 44, 46, 48],['Sc', 45],['Ti', 46, 47, 48, 49, 50],['V', 50, 51],
                 ['Cr', 50, 52, 53, 54],['Mn', 55],['Fe', 54, 56, 57, 58],['Co', 59],
                 ['Ni', 58, 60, 61, 62, 64],['Cu', 63, 65],['Zn', 64, 66, 67, 68, 70],['Ga', 69, 71],
                 ['Ge', 70, 72, 73, 74, 76],['As', 75],['Se', 74, 76, 77, 78, 80, 82],['Br', 79, 81],
                 ['Kr', 78, 80, 82, 83, 84, 86],['Rb', 85, 87],['Sr', 84, 86, 87, 88],['Y', 89],
                 ['Zr', 90, 91, 92, 94, 96],['Nb', 93],['Mo', 92, 94, 95, 96, 97, 98, 100],
                 ['Tc',999],['Ru', 96, 98, 99, 100, 101, 102, 104],['Rh', 103],
                 ['Pd', 102, 104, 105, 106, 108, 110],['Ag', 107, 109],
                 ['Cd', 106, 108, 110, 111, 112, 113, 114, 116],['In', 113, 115],
                 ['Sn', 112, 114, 115, 116, 117, 118, 119, 120, 122, 124],['Sb', 121, 123],
                 ['Te', 120, 122, 123, 124, 125, 126, 128, 130],['I', 127],
                 ['Xe', 124, 126, 128, 129, 130, 131, 132, 134, 136],['Cs', 133],
                 ['Ba', 130, 132, 134, 135, 136, 137, 138],['La', 138, 139],['Ce', 136, 138, 140, 142],
                 ['Pr', 141],['Nd', 142, 143, 144, 145, 146, 148, 150],['Pm',999],
                 ['Sm', 144, 147, 148, 149, 150, 152, 154],['Eu', 151, 153],
                 ['Gd', 152, 154, 155, 156, 157, 158, 160],['Tb', 159],
                 ['Dy', 156, 158, 160, 161, 162, 163, 164],['Ho', 165],
                 ['Er', 162, 164, 166, 167, 168, 170],['Tm', 169],['Yb', 168, 170, 171, 172, 173, 174, 176],
                 ['Lu', 175, 176],['Hf', 174, 176, 177, 178, 179, 180],['Ta', 180, 181],
                 ['W', 180, 182, 183, 184, 186],['Re', 185, 187],['Os', 184, 186, 187, 188, 189, 190, 192],
                 ['Ir', 191, 193],['Pt', 190, 192, 194, 195, 196, 198],['Au', 197],
                 ['Hg', 196, 198, 199, 200, 201, 202, 204],['Tl', 203, 205],['Pb', 204, 206, 207, 208],
                 ['Bi', 209],['Th', 232],['U',235,238],['Po',999],['At',999],['Rn',999],['Fr',999],['Ra',999],['Ac',999],['Pa',999],['Np',999],['Pu',999],['Am',999],['Cm',999],['Bk',999],['Cf',999],]


    def _stable_names(self):
        '''
        This private method extracts the element names from stable_el.
        Note that stable_names is a misnomer as stable_el also contains
        unstable element names with a number 999 for the *stable* mass
        numbers. (?!??)

        '''
        stable_names=[]
        for i in range(len(self.stable_el)):
            stable_names.append(self.stable_el[i][0])
        self.stable_names=stable_names

    def _process_abundance_vector(self, a, z, isomers, yps):
        '''
        This private method takes as input one vector definition and
        processes it, including sorting by charge number and
        mass number. It returns the processed input variables
        plus an element and isotope vector and a list of
        isomers.

        '''
        def cmp_to_key(mycmp):
            'Convert a cmp= function into a key= function'
            class K(object):
                def __init__(self, obj, *args):
                    self.obj = obj
                def __lt__(self, other):
                    return mycmp(self.obj, other.obj) < 0
                def __gt__(self, other):
                    return mycmp(self.obj, other.obj) > 0
                def __eq__(self, other):
                    return mycmp(self.obj, other.obj) == 0
                def __le__(self, other):
                    return mycmp(self.obj, other.obj) <= 0  
                def __ge__(self, other):
                    return mycmp(self.obj, other.obj) >= 0
                def __ne__(self, other):
                    return mycmp(self.obj, other.obj) != 0
            return K
        
        tmp=[]
        isom=[]
        for i in range(len(a)):
            if z[i]!=0 and isomers[i]==1: #if its not 'NEUt and not an isomer'
                tmp.append([self.stable_names[int(z[i])]+'-'+str(int(a[i])),yps[i],z[i],a[i]])
            elif isomers[i]!=1: #if it is an isomer
                if yps[i]==0:
                    isom.append([self.stable_names[int(z[i])]+'-'+str(int(a[i]))+'-'+str(int(isomers[i]-1)),1e-99])
                else:
                    isom.append([self.stable_names[int(z[i])]+'-'+str(int(a[i]))+'-'+str(int(isomers[i]-1)),yps[i]])
        tmp.sort(key = cmp_to_key(self.compar))
        tmp.sort(key = cmp_to_key(self.comparator))
        abunds=[]
        isotope_to_plot=[]
        z_iso_to_plot=[]
        a_iso_to_plot=[]
        el_iso_to_plot=[]
        for i in range(len(tmp)):
            isotope_to_plot.append(tmp[i][0])
            abunds.append(tmp[i][1])
            z_iso_to_plot.append(int(tmp[i][2]))
            a_iso_to_plot.append(int(tmp[i][3]))
            el_iso_to_plot.append(self.stable_names[int(tmp[i][2])])

        return a_iso_to_plot,z_iso_to_plot,abunds,isotope_to_plot,el_iso_to_plot,isom
    

    
    def compar(self, x, y):
        '''
        simple comparator method

        '''

        indX=0
        indY=0

        a= int(x[0].split('-')[1])

        b= int(y[0].split('-')[1])


        if a>b:
            return 1
        if a==b:
            return 0
        if a<b:
            return -1

    def comparator(self, x, y):
        '''
        simple comparator method

        '''

        indX=0
        indY=0
        for i in range(len(self.stable_names)):
            if self.stable_names[i] == x[0].split('-')[0]:
                indX=i
            if self.stable_names[i] == y[0].split('-')[0]:
                indY=i

        if indX>indY:
            return 1
        if indX==indY:
            return 0
        if indX<indY:
            return -1

    def _read_isotopedatabase(self, ffname='isotopedatabase.txt'):
        '''
        This private method reads the isotopedatabase.txt file in sldir
        run dictory and returns z, a, elements, the cutoff mass for each
        species that delineate beta+ and beta- decay and the logical in
        the last column. Also provides charge_from_element dictionary
        according to isotopedatabase.txt.

        '''
        name=self.sldir+ffname
        z_db, a_db, el_db, stable_a_db,logic_db=\
            np.loadtxt(name,unpack=True,dtype='str')
        z_db=np.array(z_db,dtype='int')
        a_db=np.array(a_db,dtype='int')
        stable_a_db=np.array(stable_a_db,dtype='int')

        # charge number for element name from dictionary in isotopedatabase.txt
        charge_from_element_name={}
        for name in self.stable_names:
            if name=='Neutron' or name=='Neut' or name=='NEUT' or name=='N-1':
                name='nn'
            try:
                zz=z_db[np.where(el_db==name)][0]
                charge_from_element_name[name]=zz
            except IndexError:
                print(name+" does not exist in this run")
        return z_db, a_db, el_db, stable_a_db,logic_db,charge_from_element_name

    def decay_indexpointer(self):
        '''
        This private method provides decay indexpointers which allow to
        instantaneously decay an abundance vector. These are attributes
        are.

        decay_idp : list
            points in the iso_to_plot (i.e. the undecayed abundance
            vector index space) to the decay target.
        idp_to_stables_in_isostoplot : list
            points to the stable isotopes in the undecayed abundance
            vector index space.

        Notes
        -----
        For an application example see ppn.py-abu_vector-_getcycle.

        '''
        a_iso_to_plot   =self.a_iso_to_plot
        isotope_to_plot =self.isotope_to_plot
        z_iso_to_plot   =self.z_iso_to_plot
        el_iso_to_plot  =self.el_iso_to_plot
        abunds          =self.abunds
        isom            =self.isom

        z_db, a_db, el_db, stable_a_db,logic_db,charge_from_element_name=\
            self._read_isotopedatabase()
        # find out which species  beta+ and which beta- decay:
        beta=np.sign(stable_a_db-a_db) # if a species is unstable and if beta < 0 => beta- decay
                                    # else beta > 0 => beta+ decay

        # now we need an index array on the scale of the abundance
        # distribution to be plotted that points to itself for stable species,
        # and to the stable element to which it decays in case of an unstable
        # species
        decay_index_pointer=np.zeros(len(isotope_to_plot), dtype='int')-1
        idp_to_stables_in_isostoplot=[]
        for i in range(len(isotope_to_plot)):
            element_name=isotope_to_plot[i].split('-')[0]
            try:
                stable_a=stable_a_db[np.where(el_db==element_name)][0] # 4th column for that element in isotopedatabase.txt
            except IndexError:
                print("Can't find element "+element_name+" in isotopedatabase.txt")
            if a_iso_to_plot[i] <= 209 and stable_a <=209:  # Bi209 is last stable element
                stable_mass_numbers=self.stable_el[self.stable_names.index(element_name)][1:]
                iso_db_index_range_el=np.where(el_db==element_name)
                beta_for_this_species=\
                    int(beta[iso_db_index_range_el][np.where(a_db[iso_db_index_range_el]==a_iso_to_plot[i])])
                if beta_for_this_species == 0:  # if there are no stable species for an element (Tc,Pm) the cutoff specifies
                    beta_for_this_species = -1  # the lowest mass beta- isotope
                if a_iso_to_plot[i] in stable_mass_numbers:
                    # print isotope_to_plot[i]+" is stable"
                    decay_index_pointer[i]=i
                    idp_to_stables_in_isostoplot.append(i)
                elif  a_iso_to_plot[i]==8: # Be8 -> He4
                    decay_index_pointer[i]=isotope_to_plot.index('He-4')
                else: # beta decay
                    found_decay_target=False
                    i_search=-1*beta_for_this_species
                    while not found_decay_target:
                        try:
                            try_target_el=self.stable_names[charge_from_element_name[element_name]+i_search]
                        except TypeError:
                            print("Maybe information about species "+isotope_to_plot[i]+" is not available in isotopedatabase.txt")
                            decay_index_pointer[i]=-1
                            break
                        # print try_target_el
                        try:
                            stable_mass_numbers=self.stable_el[self.stable_names.index(try_target_el)][1:]
                        except ValueError:
                            print("Can not find decay target for "+isotope_to_plot[i])
                        if a_iso_to_plot[i] in stable_mass_numbers:
                            ind_range=np.where(np.array(el_iso_to_plot)==try_target_el)[0]
                            if a_iso_to_plot[i] in np.array(a_iso_to_plot)[ind_range]:
                                this_ind=\
                                    ind_range[np.where(np.array(a_iso_to_plot)[ind_range]==a_iso_to_plot[i])[0]]
                                # print isotope_to_plot[i]+" is unstable and decays to "+isotope_to_plot[this_ind]
                                decay_index_pointer[i]=this_ind
                            else:
                                print("It seems unstable species "+isotope_to_plot[i]+" wants to decay to " \
                                    +try_target_el+"-"+str(a_iso_to_plot[i])+", however this species is not in this run." \
                                    +" This points to an inconsistency in the network build. Here we will ignore the abundance of " \
                                    +isotope_to_plot[i]+'.')
                                decay_index_pointer[i]=-1
                            found_decay_target=True
                        else:
                            i_search += -1*beta_for_this_species
        if self.debug:
            print("Decay rules:")
            for i in range(len(isotope_to_plot)):
                if decay_index_pointer[i]>= 0:
                    print(isotope_to_plot[i]+" -> "+isotope_to_plot[decay_index_pointer[i]])
        ind_tmp=idp_to_stables_in_isostoplot
        #ind_tmp=utils.strictly_monotonic(decay_index_pointer)  # this would do the same, but the method above is more straight forward

        self.decay_idp=decay_index_pointer
        self.idp_to_stables_in_isostoplot=ind_tmp

    def is_stable(self,species):
        '''
        This routine accepts input formatted like 'He-3' and checks with
        stable_el list if occurs in there.  If it does, the routine
        returns True, otherwise False.

        Notes
        -----
        this method is designed to work with an se instance from
        nugridse.py. In order to make it work with ppn.py some
        additional work is required.

        FH, April 20, 2013.

        '''
        element_name_of_iso            = species.split('-')[0]
        try:
            a_of_iso               = int(species.split('-')[1])
        except ValueError: # if the species name contains in addition to the
                           # mass number some letters, e.g. for isomere, then
                           # we assume it is unstable. This is not correct but
                           # related to the fact that in nugridse.py we do not
                           # identify species properly by the three numbers A, Z
                           # and isomeric_state. We should do that!!!!!!
            a_of_iso               = 999
        idp_of_element_in_stable_names = self.stable_names.index(element_name_of_iso)
        if a_of_iso in self.stable_el[idp_of_element_in_stable_names][1:]:
            return True
        else:
            return False

class iniabu(Utils):
    '''
    This class in the utils package reads an abundance distribution file
    of the types iniab.dat or iso_massf. It then provides you with methods to change
    some abundances, modify, normalise and eventually write out the
    final distribution in a format that can be used as an initial
    abundance file for ppn. This class also contains a method to write
    initial abundance files for a MESA run, for a given MESA netowrk.

    '''
    # clean variables that we will use in this class


    filename = ''
    
    
    def __init__(self,filename,filetype):
        if filetype == 'iniabu.ppn':
            f0=open(filename)
            sol=f0.readlines()
            f0.close

        # Now read in the whole file and create a hashed array:
            names=[]
            z=[]
            yps=np.zeros(len(sol))
            mass_number=np.zeros(len(sol))
            for i in range(len(sol)):
                z.append(int(sol[i][1:3]))
                names.extend([sol[i].split("         ")[0][4:]])
                yps[i]=float(sol[i].split("         ")[1])
               
                try:
                    mass_number[i]=int(names[i][2:5])
                except ValueError:
                    print("WARNING:")
                    print("This initial abundance file uses an element name that does")
                    print("not contain the mass number in the 3rd to 5th position.")
                    print("It is assumed that this is the proton and we will change")
                    print("the name to 'h   1' to be consistent with the notation used")
                    print("in iniab.dat files")
                    names[i]='h   1'
                mass_number[i]=int(names[i][2:5])
            # now zip them together:
            hash_abu={}
            hash_index={}
            for a,b in zip(names,yps):
                hash_abu[a] = b

            for i in range(len(names)):
                hash_index[names[i]] = i

            self.z=z
            self.abu=yps
            self.a=mass_number
            self.names=names
            self.habu=hash_abu
            self.hindex=hash_index



        if filetype == 'iso_massf':  #this has only been tested with write and set_and_normalize functions within this class.
                                     # If you have an issue contact Ondrea.
            
            f0=open(filename)
            ppn_out=f0.readlines()
            f0.close

        # Now read in the whole file and create a hashed array:
            names1=[]
            z1=[]
            yps=np.zeros(len(ppn_out))
            mass_number=np.zeros(len(ppn_out))
            yps1=[]

            isomers_to_remove = ['AL*26', 'KR*85', 'CD*15', 'LU*76', 'TAg80']            
            for i in range(7,len(ppn_out)): #skip the header and NEUT
                element_name = ppn_out[i].split("         ")[0][37:42]
                if element_name in isomers_to_remove:
                    print('Skipped %s' % element_name)
                    continue

                z1.append(float(ppn_out[i][9:12]))
                names1.extend([element_name])
                yps1.append(ppn_out[i].split("         ")[0][24:35])
             

            #convert data
            yps=np.zeros(len(yps1))
            z=np.zeros(len(z1),dtype=int)
        
            for j in range(len(z1)):
                z[j]= int(z1[j])
                yps[j]= float(yps1[j])

            names=[k.lower() for k in names1]
            names[0] = 'NEUT'
            names[1] = 'h   1'

            hash_abu={}
            hash_index={}
            for a,b in zip(names,yps):
                hash_abu[a] = b

            for i in range(len(names)):
                hash_index[names[i]] = i
            

          #  print(names)
            self.z=z
            self.abu=yps
            self.a=mass_number
            self.names=names
            self.habu=hash_abu
            self.hindex=hash_index




           
    def write(self, outfile='initial_abundance.dat',
              header_string='initial abundances for a PPN run'):
        '''
        Write initial abundance file (intended for use with ppn)

        Parameters
        ----------
        outfile : string
            Name of output file.  The default is
            'initial_abundance.dat'.
        header_string : string
            A string with header line.  The default is
            'initial abundances for a PPN run'.

        '''

        dcols=['Z', 'species','mass fraction']
        data=[self.z,self.names,self.abu]
        hd=[header_string]
        att.write(outfile,hd,dcols,data)

    def write_mesa(self, mesa_isos_file='isos.txt',
                   add_excess_iso='fe56', outfile='xa_iniabu.dat',
                   header_string='initial abundances for a MESA run',
                   header_char='!'):
        '''
        Write initial abundance file, returns written abundances and
        mesa names.

        Parameters
        ----------
        mesa_isos_file : string, optional
            List with isos copied from mesa network definition file in
            mesa/data/net_data/nets.  The default is 'isos.txt'.
        add_excess_iso : string, optional
            Add 1.-sum(isos in mesa net) to this isotope.  The defualt
            is 'fe56'.
        outfile : string, optional
            name of output file.  The default file is 'xa_iniabu.dat'.
        header_string : string, optional
            Srting with header line.  The default is
            'initial abundances for a MESA run'.
        header_char : character, optional
            The default is '!'.

        Examples
        --------

        >>> from NuGridPy import utils
        >>> !ls ~/PPN/forum.astro.keele.ac.uk/frames/mppnp/USEEPP/   # find ppn initial abundance file
        >>> !cat ~/mesa/data/net_data/nets/agb.net                   # find isos needed in mesa net
        >>> !cat > isos.txt                                          # paste needed isos into file
        >>> help(utils.iniabu)                                       # check documentation of method
        >>> x=utils.iniabu('path_to_here/forum.astro.keele.ac.uk/frames/mppnp/USEEPP/iniab2.0E-02GN93.ppn')
        >>> x.write_mesa?
        >>> mnames,mabus = x.write_mesa(add_excess_iso='ne22',
        ...                header_string='mppnp/USEEPP/iniab2.0E-02GN93.ppn for mesa/agb.net',
        ...                outfile='xa_2.0E-02GN93.mesa')

        '''

        
        f=open('isos.txt')
        a=f.readlines()
        isos=[]
        for i in range(len(a)):
            isos.append(a[i].strip().rstrip(','))

        mesa_names=[]
        abus=[]
        for i in range(len(self.z)):
            b=self.names[i].split()
            a=''
            a=a.join(b)
            if a in isos:
                mesa_names.append(a)
                abus.append(self.abu[i])
            # mesa_names.append(elements_names[int(x.z[i])].lower()+str(int(x.a[i])))

        for i in range(len(isos)):
            if isos[i] not in mesa_names:
                mesa_names.append(isos[i])
                abus.append(0.0)

        excess=1.-np.sum(np.array(abus))
        abus=np.array(abus)
        abus[mesa_names.index(add_excess_iso)] += excess

        dcols=['','']
        data=[mesa_names,abus]
        hd=[header_string]
        att.write(outfile,hd,dcols,data,header_char=header_char)
        return mesa_names,abus
        
    def set_and_normalize(self,species_hash):
        '''
        species_hash is a hash array in which you provide abundances
        referenced by species names that you want to set to some
        particular value; all other species are then normalised so that
        the total sum is 1.

        Examples
        --------

        You can set up the argument array for this method for example
        in the following way.

        >>> sp={}
        >>> sp['he  4']=0.2
        >>> sp['h   1']=0.5

        '''
        
        sum_before = sum(self.abu)
        for i in range(len(species_hash)):
            sum_before -= self.abu[self.hindex[list(species_hash.keys())[i]]]
        print("sum_before = "+str(sum_before))
        normalization_factor=old_div(1.0-sum(species_hash.values()),sum_before)
        print("normalizing the rest witih factor "+str(normalization_factor))
        self.abu *= normalization_factor 
        for i in range(len(species_hash)):
            self.abu[self.hindex[list(species_hash.keys())[i]]]=list(species_hash.values())[i]
                  
        for l in range(len(self.abu)):
            if self.abu[l] <= 1e-99:   #otherwise we might write e-100 which will be read as e-10 by ppn
                self.abu[l] = 1.0e-99
        for name in self.habu:
            self.habu[name]=self.abu[self.hindex[name]]
    


    
    def isoratio_init(self,isos):
        '''
        This file returns the isotopic ratio of two isotopes specified
        as iso1 and iso2. The isotopes are given as, e.g.,
        ['Fe',56,'Fe',58] or ['Fe-56','Fe-58'] (for compatibility)
        -> list.

        '''
        if len(isos) == 2:
            dumb = []
            dumb = isos[0].split('-')
            dumb.append(isos[1].split('-')[0])
            dumb.append(isos[1].split('-')[1])
            isos = dumb
        ssratio = old_div(self.habu[isos[0].ljust(2).lower() + str(int(isos[1])).rjust(3)], self.habu[isos[2].ljust(2).lower() + str(int(isos[3])).rjust(3)])
        return ssratio

    def iso_abundance(self,isos):
        '''
        This routine returns the abundance of a specific isotope.
        Isotope given as, e.g., 'Si-28' or as list
        ['Si-28','Si-29','Si-30']

        '''
        if type(isos) == list:
            dumb = []
            for it in range(len(isos)):
                dumb.append(isos[it].split('-'))
            ssratio = []
            isos = dumb
            for it in range(len(isos)):
                ssratio.append(self.habu[isos[it][0].ljust(2).lower() + str(int(isos[it][1])).rjust(3)])
        else:
            isos = isos.split('-')
            ssratio = self.habu[isos[0].ljust(2).lower() + str(int(isos[1])).rjust(3)]
        return ssratio


def trajectory_SgConst(Sg=0.1, delta_logt_dex=-0.01):
    '''
    setup trajectories for constant radiation entropy.

    S_gamma/R where the radiation constant R = N_A*k
    (Dave Arnett, Supernova book, p. 212)
    This relates rho and T but the time scale for this
    is independent.

    Parameters
    ----------
    Sg : float
        S_gamma/R, values between 0.1 and 10. reflect conditions in
        massive stars.  The default is 0.1.
    delta_logt_dex : float
        Sets interval between time steps in dex of logtimerev.  The
        default is -0.01.

    '''

    # reverse logarithmic time
    logtimerev=np.arange(5.,-6.,delta_logt_dex)
    logrho=np.linspace(0,8.5,len(logtimerev))
    logT = (old_div(1.,3.))*(logrho + 21.9161 + np.log10(Sg))

    #rho_6=10**logrho/(0.1213*1.e6)
    #T9=rho_6**(1./3.)
    #logT_T3=np.log10(T9*1.e9)


    pl.close(3);pl.figure(3);pl.plot(logrho,logT,label='$S/\mathrm{N_Ak}='+str(Sg)+'$')
    pl.legend(loc=2);pl.xlabel('$\log \\rho$'); pl.ylabel('$\log T$')
    pl.close(5);pl.figure(5);pl.plot(logtimerev, logrho)
    pl.xlabel('$\log (t_\mathrm{final}-t)$'); pl.ylabel('$\log \\rho$')
    pl.xlim(8,-6)

    pl.close(6);pl.figure(6);pl.plot(logtimerev)
    pl.ylabel('$\log (t_\mathrm{final}-t)$'); pl.xlabel('cycle')

    # [t] logtimerev yrs
    # [rho] cgs
    # [T]   K

    T9=old_div(10**logT,1.e9)
    data=[logtimerev,T9,logrho]
    att.writeTraj(filename='trajectory.input', data=data, ageunit=2, tunit=1, rhounit=1, idNum=1)

    # data: A list of 1D data vectors with time, T and rho                                                       # ageunit: If 1 ageunit = SEC, If 0 ageunit = YRS. If 2 agunit = logtimerev in yrs. Default is 0             # logtimerev is log of time until end


def _xlimrev(self):
    ''' reverse xrange'''
    xmax,xmin=pyl.xlim()
    pyl.xlim(xmin,xmax)

def close_wins(win_min,win_max):
    '''
    close all windows in a certain window number range

    win_min/max  minumum and maximum window number to close
    '''

    for i in range(win_min,win_max+1):
        close(i)

def _xlimrev():
    ''' reverse xrange'''
    xmax,xmin=pl.xlim()
    pl.xlim(xmin,xmax)

def species_list(what_list):
    '''
    provide default lists of elements to plot.

    what_list : string
        String name of species lists provided.

        If what_list is "CNONe", then C, N, O and some other light
        elements.

        If what_list is "s-process", then s-process indicators.

    '''
    if what_list is "CNONe":
        list_to_print = ['H-1','He-4','C-12','N-14','O-16','Ne-20']
    elif what_list is "sprocess":
        list_to_print = ['Fe-56','Ge-70','Zn-70','Se-76','Kr-80','Kr-82','Kr-86','Sr-88','Ba-138','Pb-208']
    elif what_list is "burn_stages":
        list_to_print = ['H-1','He-4','C-12','O-16','Ne-20','Si-28']
    elif what_list is "list_marco_1":
        list_to_print = ['C-12','O-16','Ne-20','Ne-22','Na-23','Fe-54','Fe-56','Zn-70','Ge-70','Se-76','Kr-80','Kr-82','Sr-88','Y-89','Zr-96','Te-124','Xe-130','Xe-134','Ba-138']

    return list_to_print

def linestyle(i,a=5,b=3):
    '''
    provide one out of 25 unique combinations of style, color and mark

    use in combination with markevery=a+mod(i,b) to add spaced points,
    here a would be the base spacing that would depend on the data
    density, modulated with the number of lines to be plotted (b)

    Parameters
    ----------
    i : integer
        Number of linestyle combination - there are many....
    a : integer
        Spacing of marks.  The default is 5.
    b : integer
        Modulation in case of plotting many nearby lines.  The default
        is 3.

    Examples
    --------

    >>> plot(x,sin(x),linestyle(7)[0], markevery=linestyle(7)[1])


    (c) 2014 FH
    '''

    lines=['-','--','-.',':']
    points=['v','^','<','>','1','2','3','4','s','p','*','h','H','+','x','D','d','o']
    colors=['b','g','r','c','m','k']
    ls_string = colors[sc.mod(i,6)]+lines[sc.mod(i,4)]+points[sc.mod(i,18)]
    mark_i    = a+sc.mod(i,b)
    return ls_string,int(mark_i)

def colourblind(i):
    '''
        colour pallete from http://tableaufriction.blogspot.ro/
        allegedly suitable for colour-blind folk

        SJ
    '''

    rawRGBs = [(162,200,236),
               (255,128,14),
               (171,171,171),
               (95,158,209),
               (89,89,89),
               (0,107,164),
               (255,188,121),
               (207,207,207),
               (200,82,0),
               (137,137,137)]

    scaledRGBs = []
    for r in rawRGBs:
        scaledRGBs.append((old_div(r[0],255.),old_div(r[1],255.),old_div(r[2],255.)))

    idx = sc.mod(i,len(scaledRGBs))
    return scaledRGBs[idx]

def colourblind2(i):
    '''
        another colour pallete from http://www.sron.nl/~pault/
        allegedly suitable for colour-blind folk

        SJ
    '''

    hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',
               '#CC6677', '#882255', '#AA4499']
    idx = sc.mod(i,len(hexcols))
    return hexcols[idx]

def linestylecb(i,a=5,b=3):
    '''
        version of linestyle function with colourblind colour scheme

        returns linetyle, marker, color (see example)

        Parameters
        ----------
        i : integer
        Number of linestyle combination - there are many....
        a : integer
        Spacing of marks.  The default is 5.
        b : integer
        Modulation in case of plotting many nearby lines.  The default
        is 3.

        Examples
        --------

        >>> plot(x,sin(x),ls=linestyle(7)[0], marker=linestyle(7)[1], \
                 color=linestyle(7)[2],markevery=linestyle(7)[3])


        (c) 2014 FH
        '''

    lines=['-','--','-.',':']
    points=['v','^','<','>','1','2','3','4','s','p','*','h','H','+','x','D','d','o']
    colors=['b','g','r','c','m','k']
    col=colourblind(i)
    style=lines[sc.mod(i,4)]
    point=points[sc.mod(i,18)]
    mark_i    = a+sc.mod(i,b)
    return style,point,col,mark_i


def symbol_list(what_list):
    '''
    provide default symbol lists

    Parameters
    ----------
    what_list : string
        String name of symbol lists provided; "list1", "list2",
        "lines1" or "lines2".

    '''
    if what_list is "list1":
        symbol=['ro','bo','ko','go','mo'\
                ,'r-','b-','k-','g-','m-','r--','b--','k--'\
                ,'g--','r1']
        #symbol=['r+','ro','r-']
    elif what_list is "list2":
        symbol=['r-','b--','g-.','k:','md','.','o','v','^','<','>','1','2',\
                '3','4','s','p','*','h','H','+']
    elif what_list is "lines1":
        symbol=['b--','k--','r--','c--','m--','g--','b-','k-','r-','c-','m-','g-','b.','b-.','k-.','r-.','c-.','m-.','g-.','b:','k:','r:','c:','m:','g:']
    elif what_list is "lines2":
        symbol=['g:','r-.','k-','b--','k-.','b+','r:','b-','c--','m--','g--','r-','c-','m-','g-','k-.','c-.','m-.','g-.','k:','r:','c:','m:','b-.','b:']
    return symbol

def make_list(default_symbol_list, len_list_to_print):
    '''
    provide the list of symbols to use according for the list of
    species/arrays to plot.

    Parameters
    ----------
    default_symbol_list : list
        Symbols that the user choose to use.
    len_list_to_print : integer
        len of list of species/arrays to print.

    '''

    symbol_used = []
    for i in range(len_list_to_print):
        symbol_used.append(default_symbol_list[sc.mod(i,len(default_symbol_list))])

    return symbol_used

def strictly_monotonic(bb):
    '''
    bb is an index array which may have numerous double or triple
    occurrences of indices, such as for example the decay_index_pointer.
    This method removes all entries <= -, then all dublicates and
    finally returns a sorted list of indices.

    '''
    cc=bb[np.where(bb>=0)]
    cc.sort()
    dc=cc[1:]-cc[:-1] # subsequent equal entries have 0 in db
    dc=np.insert(dc,0,1) # the first element is always unique (the second occurence is the dublicate)
    dc_mask=np.ma.masked_equal(dc,0)
    return np.ma.array(cc,mask=dc_mask.mask).compressed()

def solar(filename_solar, solar_factor):
    '''
    read solar abundances from filename_solar.

    Parameters
    ----------
    filename_solar : string
        The file name.
    solar_factor : float
        The correction factor to apply, in case filename_solar is not
        solar, but some file used to get initial abundances at
        metallicity lower than solar. However, notice that this is
        really rude, since alpha-enahncements and things like that are
        not properly considered.  Only H and He4 are not multiplied. So,
        for publications PLEASE use proper filename_solar at...solar,
        and use solar_factor = 1. Marco

    '''

    f0=open(filename_solar)
    sol=f0.readlines()
    f0.close
    sol[0].split("         ")

    # Now read in the whole file and create a hashed array:
    global names_sol
    names_sol=[]
    global z_sol
    z_sol=[]
    yps=np.zeros(len(sol))
    mass_number=np.zeros(len(sol))
    for i in range(len(sol)):
        z_sol.append(int(sol[i][1:3]))
        names_sol.extend([sol[i].split("         ")[0][4:]])
        yps[i]=float(sol[i].split("         ")[1]) * solar_factor
        try:
            mass_number[i]=int(names_sol[i][2:5])
        except ValueError:
            print("WARNING:")
            print("This initial abundance file uses an element name that does")
            print("not contain the mass number in the 3rd to 5th position.")
            print("It is assumed that this is the proton and we will change")
            print("the name to 'h   1' to be consistent with the notation used in")
            print("iniab.dat files")
            names_sol[i]='h   1'
            mass_number[i]=int(names_sol[i][2:5])
        if mass_number[i] == 1 or mass_number[i] == 4:
            yps[i] = old_div(yps[i],solar_factor)
    #  convert 'h   1' in prot, not needed any more??
    #names_sol[0] = 'prot '


    # now zip them together:
    global solar_abundance
    solar_abundance={}
    for a,b in zip(names_sol,yps):
        solar_abundance[a] = b



    z_bismuth = 83
    global solar_elem_abund
    solar_elem_abund = np.zeros(z_bismuth)


    for i in range(z_bismuth):
        dummy = 0.
        for j in range(len(solar_abundance)):
            if z_sol[j] == i+1:
                dummy = dummy + float(solar_abundance[names_sol[j]])
        solar_elem_abund[i] = dummy


def convert_specie_naming_from_h5_to_ppn(isotope_names):
    '''
    read isotopes names from h5 files, and convert them
    according to standard scheme used inside ppn and mppnp.  Also
    Z and A are recalculated, for these species. Isomers are
    excluded for now, since there were recent changes in isomers
    name. As soon as the isomers names are settled, than Z and A
    provided here will be obsolete, and can be changed by usual Z
    and A.

    '''

    spe_rude1 = []
    spe_rude2 = []
    spe_rude3 = []
    for i in range(len(isotope_names)):
        spe_rude1.append(isotope_names[i].split('-')[0])
        spe_rude2.append(isotope_names[i].split('-')[1])
    # spe_rude1 is elem name and spe_rude2 is mass number.
    #print spe_rude1,spe_rude2
    k = 0
    for i in range(len(spe_rude1)):
        try:
            if int(spe_rude2[i]) < 10:
                spe_rude3.append(str(spe_rude1[i][0:2])+str('  ')+str(spe_rude2[i][0:3]))
            elif int(spe_rude2[i]) >= 10 and int(spe_rude2[i]) < 100 :
                spe_rude3.append(str(spe_rude1[i][0:2])+str(' ')+str(spe_rude2[i][0:3]))
            elif int(spe_rude2[i]) >= 100 :
                spe_rude3.append(str(spe_rude1[i][0:2])+str(spe_rude2[i][0:3]))
        except ValueError:
            k = k+1
            None

    global spe
    spe = []
    global n_array
    n_array = []
    for i in range(len(spe_rude3)):
        if len(str(spe_rude1[i])) == 1:
            spe.append(str(spe_rude3[i][0:1])+str(' ')+str(spe_rude3[i][1:4]))
        else:
            spe.append(spe_rude3[i])
        n_array.append(i)
    if spe[0]=='Ne  1':
        spe[0] = 'N   1'

    # spe_rude is the isotope name, in agreement with what we use in ppn, etc.
    # need to do this to can use other functions without changing them drastically.



    # here I skip isomers...
    global amass_int
    amass_int=np.zeros(len(spe_rude2))
    for i in range(len(spe_rude2)-k):
        amass_int[i]=int(spe_rude2[i])
        #print amass_int


    # here I have to create an array for the atomic numbers.
    # I need to this when I calculate and plot element abundances

    global znum_int
    znum_int=np.zeros(len(spe))

    for i in range(len(spe)):
        znum_int[i] = Utils.elements_names.index(str(spe[i][0:2]).strip())
        # changed by alex
        # if str(spe[i][0:2]) == 'H ':
        #     znum_int[i] = 1
        # elif str(spe[i][0:2]) == 'He':
        #     znum_int[i] = 2
        # elif str(spe[i][0:2]) == 'Li':
        #     znum_int[i] = 3
        # elif str(spe[i][0:2]) == 'Be':
        #     znum_int[i] = 4
        # elif str(spe[i][0:2]) == 'B ':
        #     znum_int[i] = 5
        # elif str(spe[i][0:2]) == 'C ':
        #     znum_int[i] = 6
        # elif str(spe[i][0:2]) == 'N ':
        #     znum_int[i] = 7
        # elif str(spe[i][0:2]) == 'O ':
        #     znum_int[i] = 8
        # elif str(spe[i][0:2]) == 'F ':
        #     znum_int[i] = 9
        # elif str(spe[i][0:2]) == 'Ne':
        #     znum_int[i] = 10
        # elif str(spe[i][0:2]) == 'Na':
        #     znum_int[i] = 11
        # elif str(spe[i][0:2]) == 'Mg':
        #     znum_int[i] = 12
        # elif str(spe[i][0:2]) == 'Al':
        #     znum_int[i] = 13
        # elif str(spe[i][0:2]) == 'Si':
        #     znum_int[i] = 14
        # elif str(spe[i][0:2]) == 'P ':
        #     znum_int[i] = 15
        # elif str(spe[i][0:2]) == 'S ':
        #     znum_int[i] = 16
        # elif str(spe[i][0:2]) == 'Cl':
        #     znum_int[i] = 17
        # elif str(spe[i][0:2]) == 'Ar':
        #     znum_int[i] = 18
        # elif str(spe[i][0:2]) == 'K ':
        #     znum_int[i] = 19
        # elif str(spe[i][0:2]) == 'Ca':
        #     znum_int[i] = 20
        # elif str(spe[i][0:2]) == 'Sc':
        #     znum_int[i] = 21
        # elif str(spe[i][0:2]) == 'Ti':
        #     znum_int[i] = 22
        # elif str(spe[i][0:2]) == 'V ':
        #     znum_int[i] = 23
        # elif str(spe[i][0:2]) == 'Cr':
        #     znum_int[i] = 24
        # elif str(spe[i][0:2]) == 'Mn':
        #     znum_int[i] = 25
        # elif str(spe[i][0:2]) == 'Fe':
        #     znum_int[i] = 26
        # elif str(spe[i][0:2]) == 'Co':
        #     znum_int[i] = 27
        # elif str(spe[i][0:2]) == 'Ni':
        #     znum_int[i] = 28
        # elif str(spe[i][0:2]) == 'Cu':
        #     znum_int[i] = 29
        # elif str(spe[i][0:2]) == 'Zn':
        #     znum_int[i] = 30
        # elif str(spe[i][0:2]) == 'Ga':
        #     znum_int[i] = 31
        # elif str(spe[i][0:2]) == 'Ge':
        #     znum_int[i] = 32
        # elif str(spe[i][0:2]) == 'As':
        #     znum_int[i] = 33
        # elif str(spe[i][0:2]) == 'Se':
        #     znum_int[i] = 34
        # elif str(spe[i][0:2]) == 'Br':
        #     znum_int[i] = 35
        # elif str(spe[i][0:2]) == 'Kr':
        #     znum_int[i] = 36
        # elif str(spe[i][0:2]) == 'Rb':
        #     znum_int[i] = 37
        # elif str(spe[i][0:2]) == 'Sr':
        #     znum_int[i] = 38
        # elif str(spe[i][0:2]) == 'Y ':
        #     znum_int[i] = 39
        # elif str(spe[i][0:2]) == 'Zr':
        #     znum_int[i] = 40
        # elif str(spe[i][0:2]) == 'Nb':
        #     znum_int[i] = 41
        # elif str(spe[i][0:2]) == 'Mo':
        #     znum_int[i] = 42
        # elif str(spe[i][0:2]) == 'Tc':
        #     znum_int[i] = 43
        # elif str(spe[i][0:2]) == 'Ru':
        #     znum_int[i] = 44
        # elif str(spe[i][0:2]) == 'Rh':
        #     znum_int[i] = 45
        # elif str(spe[i][0:2]) == 'Pd':
        #     znum_int[i] = 46
        # elif str(spe[i][0:2]) == 'Ag':
        #     znum_int[i] = 47
        # elif str(spe[i][0:2]) == 'Cd':
        #     znum_int[i] = 48
        # elif str(spe[i][0:2]) == 'In':
        #     znum_int[i] = 49
        # elif str(spe[i][0:2]) == 'Sn':
        #     znum_int[i] = 50
        # elif str(spe[i][0:2]) == 'Sb':
        #     znum_int[i] = 51
        # elif str(spe[i][0:2]) == 'Te':
        #     znum_int[i] = 52
        # elif str(spe[i][0:2]) == 'I ':
        #     znum_int[i] = 53
        # elif str(spe[i][0:2]) == 'Xe':
        #     znum_int[i] = 54
        # elif str(spe[i][0:2]) == 'Cs':
        #     znum_int[i] = 55
        # elif str(spe[i][0:2]) == 'Ba':
        #     znum_int[i] = 56
        # elif str(spe[i][0:2]) == 'La':
        #     znum_int[i] = 57
        # elif str(spe[i][0:2]) == 'Ce':
        #     znum_int[i] = 58
        # elif str(spe[i][0:2]) == 'Pr':
        #     znum_int[i] = 59
        # elif str(spe[i][0:2]) == 'Nd':
        #     znum_int[i] = 60
        # elif str(spe[i][0:2]) == 'Pm':
        #     znum_int[i] = 61
        # elif str(spe[i][0:2]) == 'Sm':
        #     znum_int[i] = 62
        # elif str(spe[i][0:2]) == 'Eu':
        #     znum_int[i] = 63
        # elif str(spe[i][0:2]) == 'Gd':
        #     znum_int[i] = 64
        # elif str(spe[i][0:2]) == 'Tb':
        #     znum_int[i] = 65
        # elif str(spe[i][0:2]) == 'Dy':
        #     znum_int[i] = 66
        # elif str(spe[i][0:2]) == 'Ho':
        #     znum_int[i] = 67
        # elif str(spe[i][0:2]) == 'Er':
        #     znum_int[i] = 68
        # elif str(spe[i][0:2]) == 'Tm':
        #     znum_int[i] = 69
        # elif str(spe[i][0:2]) == 'Yb':
        #     znum_int[i] = 70
        # elif str(spe[i][0:2]) == 'Lu':
        #     znum_int[i] = 71
        # elif str(spe[i][0:2]) == 'Hf':
        #     znum_int[i] = 72
        # elif str(spe[i][0:2]) == 'Ta':
        #     znum_int[i] = 73
        # elif str(spe[i][0:2]) == 'W ':
        #     znum_int[i] = 74
        # elif str(spe[i][0:2]) == 'Re':
        #     znum_int[i] = 75
        # elif str(spe[i][0:2]) == 'Os':
        #     znum_int[i] = 76
        # elif str(spe[i][0:2]) == 'Ir':
        #     znum_int[i] = 77
        # elif str(spe[i][0:2]) == 'Pt':
        #     znum_int[i] = 78
        # elif str(spe[i][0:2]) == 'Au':
        #     znum_int[i] = 79
        # elif str(spe[i][0:2]) == 'Hg':
        #     znum_int[i] = 80
        # elif str(spe[i][0:2]) == 'Tl':
        #     znum_int[i] = 81
        # elif str(spe[i][0:2]) == 'Pb':
        #     znum_int[i] = 82
        # elif str(spe[i][0:2]) == 'Bi':
        #     znum_int[i] = 83
        # elif str(spe[i][0:2]) == 'Po':
        #     znum_int[i] = 84
        # elif str(spe[i][0:2]) == 'At':
        #     znum_int[i] = 85
        # elif str(spe[i][0:2]) == 'Rn':
        #     znum_int[i] = 86
        # elif str(spe[i][0:2]) == 'Fr':
        #     znum_int[i] = 87
        # elif str(spe[i][0:2]) == 'Ra':
        #     znum_int[i] = 88
        # elif str(spe[i][0:2]) == 'Ac':
        #     znum_int[i] = 89
        # elif str(spe[i][0:2]) == 'Th':
        #     znum_int[i] = 90
        # elif str(spe[i][0:2]) == 'Pa':
        #     znum_int[i] = 91
        # elif str(spe[i][0:2]) == 'U ':
        #     znum_int[i] = 92
        # elif str(spe[i][0:2]) == 'Np':
        #     znum_int[i] = 93
        # elif str(spe[i][0:2]) == 'Pu':
        #     znum_int[i] = 94
        # elif str(spe[i][0:2]) == 'Am':
        #     znum_int[i] = 95
        # elif str(spe[i][0:2]) == 'Cm':
        #     znum_int[i] = 96
        # elif str(spe[i][0:2]) == 'Bk':
        #     znum_int[i] = 97
        # elif str(spe[i][0:2]) == 'Cf':
        #     znum_int[i] = 98

    if spe[0] == 'N   1':
        znum_int[0] = 0

    # here the index to connect name and atomic numbers.
    global index_atomic_number
    index_atomic_number = {}
    for a,b in zip(spe,znum_int):
        index_atomic_number[a]=b


def   define_zip_index_for_species(names_ppn_world,
                                   number_names_ppn_world):
    ''' This just give back cl, that is the original index as it is read from files from a data file.'''

    #connect the specie number in the list, with the specie name
    global cl
    cl={}
    for a,b in zip(names_ppn_world,number_names_ppn_world):
        cl[a] = b



def element_abund_marco(i_decay, stable_isotope_list,
                        stable_isotope_identifier,
                        mass_fractions_array_not_decayed,
                        mass_fractions_array_decayed):
    '''
    Given an array of isotopic abundances not decayed and a similar
    array of isotopic abundances not decayed, here elements abundances,
    and production factors for elements are calculated

    '''


    # this way is done in a really simple way. May be done better for sure, in a couple of loops.
    # I keep this, since I have only to copy over old script. Falk will probably redo it.

    #import numpy as np
    #from NuGridPy import utils as u

    global elem_abund
    elem_abund = np.zeros(z_bismuth)
    global elem_abund_decayed
    elem_abund_decayed = np.zeros(z_bismuth)
    global elem_prod_fac
    elem_prod_fac = np.zeros(z_bismuth)
    global elem_prod_fac_decayed
    elem_prod_fac_decayed = np.zeros(z_bismuth)


    # notice that elem_abund include all contribution, both from stables and unstables in
    # that moment.
    for i in range(z_bismuth):
        dummy = 0.
        for j in range(len(spe)):
            if znum_int[j] == i+1 and stable_isotope_identifier[j] > 0.5:
                dummy = dummy + float(mass_fractions_array_not_decayed[j])
        elem_abund[i] = dummy


    for i in range(z_bismuth):
        if index_stable[i] == 1:
            elem_prod_fac[i] = float(old_div(elem_abund[i],solar_elem_abund[i]))
        elif index_stable[i] == 0:
            elem_prod_fac[i] = 0.


    if i_decay == 2:
        for i in range(z_bismuth):
            dummy = 0.
            for j in range(len(mass_fractions_array_decayed)):
                if znum_int[cl[stable_isotope_list[j].capitalize()]] == i+1:
                #print znum_int[cl[stable[j].capitalize()]],cl[stable[j].capitalize()],stable[j]
                    dummy = dummy + float(mass_fractions_array_decayed[j])
            elem_abund_decayed[i] = dummy


        for i in range(z_bismuth):
            if index_stable[i] == 1:
                elem_prod_fac_decayed[i] = float(old_div(elem_abund_decayed[i],solar_elem_abund[i]))
            elif index_stable[i] == 0:
                elem_prod_fac_decayed[i] = 0.


def stable_specie():
    ''' provide the list of stable species, and decay path feeding stables '''


    #import numpy as np


    stable_raw=[]
    stable_raw = ['H   1', 'H   2',\
    'HE  3', 'HE  4',\
    'LI  6', 'LI  7',\
    'BE  9',\
    'B  10', 'B  11',\
    'C  12', 'C  13',\
    'N  14', 'N  15',\
    'O  16', 'O  17', 'O  18',\
    'F  19',\
    'NE 20', 'NE 21', 'NE 22',\
    'NA 23',\
    'MG 24', 'MG 25', 'MG 26',\
    'AL 27',\
    'SI 28', 'SI 29', 'SI 30',\
    'P  31',\
    'S  32', 'S  33', 'S  34', 'S  36',\
    'CL 35', 'CL 37',\
    'AR 36', 'AR 38', 'AR 40',\
    'K  39', 'K  40', 'K  41',\
    'CA 40', 'CA 42', 'CA 43', 'CA 44', 'CA 46', 'CA 48',\
    'SC 45',\
    'TI 46', 'TI 47', 'TI 48', 'TI 49', 'TI 50',\
    'V  50', 'V  51',\
    'CR 50', 'CR 52', 'CR 53', 'CR 54',\
    'MN 55',\
    'FE 54', 'FE 56', 'FE 57', 'FE 58',\
    'CO 59',\
    'NI 58', 'NI 60', 'NI 61', 'NI 62', 'NI 64',\
    'CU 63', 'CU 65',\
    'ZN 64', 'ZN 66', 'ZN 67', 'ZN 68', 'ZN 70',\
    'GA 69', 'GA 71',\
    'GE 70', 'GE 72', 'GE 73', 'GE 74', 'GE 76',\
    'AS 75',\
    'SE 74', 'SE 76', 'SE 77', 'SE 78', 'SE 80', 'SE 82',\
    'BR 79', 'BR 81',\
    'KR 78', 'KR 80', 'KR 82', 'KR 83', 'KR 84', 'KR 86',\
    'RB 85', 'RB 87',\
    'SR 84', 'SR 86', 'SR 87', 'SR 88',\
    'Y  89',\
    'ZR 90', 'ZR 91', 'ZR 92', 'ZR 94', 'ZR 96',\
    'NB 93',\
    'MO 92', 'MO 94', 'MO 95', 'MO 96', 'MO 97', 'MO 98', 'MO100',\
    'RU 96', 'RU 98', 'RU 99', 'RU100', 'RU101', 'RU102', 'RU104',\
    'RH103',\
    'PD102', 'PD104', 'PD105', 'PD106', 'PD108', 'PD110',\
    'AG107', 'AG109',\
    'CD106', 'CD108', 'CD110', 'CD111', 'CD112', 'CD113', 'CD114', 'CD116',\
    'IN113',  'IN115',\
    'SN112', 'SN114', 'SN115', 'SN116', 'SN117', 'SN118', 'SN119', 'SN120', 'SN122', 'SN124',\
    'SB121', 'SB123',\
    'TE120', 'TE122', 'TE123', 'TE124', 'TE125', 'TE126', 'TE128', 'TE130',\
    'I 127',\
    'XE124', 'XE126', 'XE128', 'XE129', 'XE130', 'XE131', 'XE132', 'XE134', 'XE136',\
    'CS133',\
    'BA130', 'BA132', 'BA134', 'BA135', 'BA136', 'BA137', 'BA138',\
    'LA138', 'LA139',\
    'CE136', 'CE138', 'CE140', 'CE142',\
    'PR141',\
    'ND142', 'ND143', 'ND144', 'ND145', 'ND146', 'ND148', 'ND150',\
    'SM144', 'SM147', 'SM148', 'SM149', 'SM150', 'SM152', 'SM154',\
    'EU151', 'EU153',\
    'GD152', 'GD154', 'GD155', 'GD156', 'GD157', 'GD158', 'GD160',\
    'TB159',\
    'DY156', 'DY158', 'DY160', 'DY161', 'DY162', 'DY163', 'DY164',\
    'HO165',\
    'ER162', 'ER164', 'ER166', 'ER167', 'ER168', 'ER170',\
    'TM169',\
    'YB168', 'YB170', 'YB171', 'YB172', 'YB173', 'YB174', 'YB176',\
    'LU175', 'LU176',\
    'HF174', 'HF176', 'HF177', 'HF178', 'HF179', 'HF180',\
    'TA180', 'TA181',\
    'W 180', 'W 182', 'W 183', 'W 184', 'W 186',\
    'RE185', 'RE187',\
    'OS184', 'OS186', 'OS187', 'OS188', 'OS189', 'OS190', 'OS192',\
    'IR191', 'IR193',\
    'PT190', 'PT192', 'PT194', 'PT195', 'PT196', 'PT198',\
    'AU197',\
    'HG196', 'HG198', 'HG199', 'HG200', 'HG201', 'HG202', 'HG204',\
    'TL203', 'TL205',\
    'PB204', 'PB206', 'PB207', 'PB208',\
    'BI209',\
    'TH232',\
    'U 235','U 238']

    jj=-1
    global count_size_stable
    count_size_stable=[]
    global stable
    stable=[]
    global jdum
    jdum=np.zeros(len(stable_raw))
    global jjdum
    jjdum=np.zeros(len(spe))
    for i in range(len(stable_raw)):
        dum_str = stable_raw[i]
        for j in range(len(spe)):
            if stable_raw[i].capitalize() == spe[j]:
                stable.append(stable_raw[i])
                jdum[i]=1
                jjdum[j]=1
                jj=jj+1
                count_size_stable.append(int(jj))
    #print stable
    # back_ind is an index to go back, to use the order of stable
    # useful for example for decayed yields.
    global back_ind
    back_ind={}
    for a,b in zip(stable,count_size_stable):
        back_ind[a]=b
    #print 'in stable:',back_ind['SE 74']
    # definition of decay paths
    global decay_raw
    decay_raw=[]
    decay_raw=[['H   1'],\
    ['H   2'],\
    ['HE  3'],\
    ['HE  4','B   8'],\
    ['LI  6'],\
    ['LI  7','BE  7'],\
    ['BE  9'],\
    ['B  10','BE 10'],\
    ['B  11','C  11','BE 11'],\
    ['C  12'],\
    ['C  13','N  13','O  13'],\
    ['N  14','C  14','O  14'],\
    ['N  15','C  15','O  15','F  15'],\
    ['O  16'],\
    ['O  17','F  17'],\
    ['O  18','F  18','NE 18'],\
    ['F  19','O  19','NE 19'],\
    ['NE 20','F  20','NA 20'],\
    ['NE 21','F  21','NA 21'],\
    ['NE 22','NA 22','MG 22','F  22'],\
    ['NA 23','MG 23','NE 23'],\
    ['MG 24','AL 24','NA 24','NE 24',],\
    ['MG 25','NA 25','AL 25'],\
    ['MG 26','SI 26','AL 26','NA 26','AL*26'],\
    ['AL 27','SI 27','MG 27'],\
    ['SI 28','AL 28','MG 28'],\
    ['SI 29','P  29','AL 29','MG 29'],\
    ['SI 30','S  30','P  30','AL 30','MG 30'],\
    ['P  31','S  31','SI 31'],\
    ['S  32','P  32','SI 32'],\
    ['S  33','CL 33','P  33','SI 33'],\
    ['S  34','CL 34','P  34'],\
    ['S  36','P  36'],\
    ['CL 35','S  35','AR 35'],\
    ['CL 37','S  37','P  37','AR 37','K  37'],\
    ['AR 36','CL 36'],\
    ['AR 38','CL 38','S  38','P  38','K  38'],\
    ['AR 40','CL 40','S  40'],\
    ['K  39','AR 39','CL 39','S  39','CA 39'],\
    ['K  40'],\
    ['K  41','AR 41','CL 41','S  41','CA 41','SC 41'],\
    ['CA 40','SC 40'],\
    ['CA 42','K  42','AR 42','CL 42','S  42','SC 42','TI 42'],\
    ['CA 43','K  43','AR 43','CL 43','SC 43','TI 43','V  43'],\
    ['CA 44','K  44','AR 44','CL 44','SC 44','TI 44'],\
    ['CA 46','K  46','AR 46'],\
    ['CA 48','K  48','AR 48'],\
    ['SC 45','CA 45','K  45','AR 45','CL 45','TI 45','V  45'],\
    ['TI 46','SC 46','V  46','CR 46'],\
    ['TI 47','SC 47','CA 47','K  47','AR 47','V  47','CR 47'],\
    ['TI 48','SC 48','V  48','CR 48'],\
    ['TI 49','SC 49','CA 49','K  49','V  49','CR 49','MN 49'],\
    ['TI 50','SC 50','CA 50','K  50'],\
    ['V  50'],\
    ['V  51','CR 51','TI 51','SC 51','CA 51','MN 51'],\
    ['CR 50','MN 50'],\
    ['CR 52','MN 52','FE 52','V  52','TI 52','SC 52','CA 52'],\
    ['CR 53','MN 53','FE 53','V  53','TI 53','SC 53'],\
    ['CR 54','MN 54','V  54','TI 54','SC 54'],\
    ['MN 55','FE 55','CR 55','V  55','TI 55','CO 55'],\
    ['FE 54','CO 54'],\
    ['FE 56','NI 56','CO 56','MN 56','CR 56'],\
    ['FE 57','NI 57','CO 57','MN 57','CR 57'],\
    ['FE 58','CO 58','MN 58','CR 58'],\
    ['CO 59','FE 59','MN 59','CR 59','NI 59','CU 59'],\
    ['NI 58','CU 58'],\
    ['NI 60','CO 60','FE 60','MN 60','CR 60','CU 60','ZN 60'],\
    ['NI 61','CO 61','FE 61','MN 61','CU 61','ZN 61'],\
    ['NI 62','CO 62','FE 62','CU 62','ZN 62'],\
    ['NI 64','CO 64','FE 64','CU 64'],\
    ['CU 63','NI 63','CO 63','FE 63','MN 63','ZN 63','GA 63'],\
    ['CU 65','NI 65','CO 65','FE 65','ZN 65','GA 65','GE 65'],\
    ['ZN 64','CU 64','GA 64','GE 64'],\
    ['ZN 66','CU 66','NI 66','CO 66','FE 66','GA 66','GE 66'],\
    ['ZN 67','CU 67','NI 67','CO 67','FE 67','GA 67','GE 67','AS 77'],\
    ['ZN 68','NI 68','CO 68','GA 68','GE 68','CU 68','AS 68','SE 68'],\
    ['ZN 70','CU 70','NI 70','CO 70'],\
    ['GA 69','ZN 69','CU 69','NI 69','GE 69','AS 69','SE 69'],\
    ['GA 71','ZN 71','CU 71','NI 71','GE 71','AS 71','SE 71','BR 71'],\
    ['GE 70','GA 70','AS 70','SE 70','BR 70'],\
    ['GE 72','GA 72','ZN 72','CU 72','NI 72','AS 72','SE 72','BR 72','KR 72'],\
    ['GE 73','GA 73','ZN 73','CU 73','NI 73','AS 73','SE 73','BR 73','KR 73'],\
    ['GE 74','GA 74','ZN 74','CU 74','NI 74','AS 74'],\
    ['GE 76','GA 76','ZN 76','CU 76'],\
    ['AS 75','GE 75','GA 75','ZN 75','CU 75','SE 75','BR 75','KR 75','RB 75'],\
    ['SE 74','AS 74','BR 74','KR 74'],\
    ['SE 76','AS 76','BR 76','KR 76','RB 76','SR 76'],\
    ['SE 77','AS 77','GE 77','BR 77','GA 77','ZN 77','KR 77','RB 77','SR 77'],\
    ['SE 78','AS 78','GE 78','GA 78','ZN 78','BR 78'],\
    ['SE 80','AS 80','GE 80','GA 80','ZN 80'],\
    ['SE 82','AS 82','GE 82','GA 82'],\
    ['BR 79','SE 79','AS 79','GE 79','GA 79','ZN 79','KR 79','RB 79','SR 79','Y  79'],\
    ['BR 81','SE 81','KR 81','AS 81','GE 81','GA 81','RB 81','SR 81','Y  81','ZR 81'],\
    ['KR 78','RB 78','SR 78','Y  78'],\
    ['KR 80','BR 80','RB 80','SR 80','ZR 80'],\
    ['KR 82','BR 82','RB 82','SR 82','Y  82','ZR 82'],\
    ['KR 83','BR 83','SE 83','AS 83','GE 83','RB 83','SR 83','Y  83','ZR 83','NB 83'],\
    ['KR 84','BR 84','SE 84','AS 84','GE 84','RB 84'],\
    ['KR 86','BR 86','SE 86','AS 86'],\
    ['RB 85','KR 85','SR 85','KR*85','BR 85','SE 85','AS 85','Y  85','ZR 85','NB 85','MO 85'],\
    ['RB 87','KR 87','BR 87','SE 87','AS 87'],\
    ['SR 84','Y  84','ZR 84','NB 84','MO 84'],\
    ['SR 86','RB 86','Y  86','ZR 86','NB 86','MO 86'],\
    ['SR 87','Y  87','ZR 87','NB 87','MO 87','TC 87'],\
    ['SR 88','RB 88','KR 88','BR 88','SE 88','Y  88','ZR 88','NB 88','MO 88','TC 88'],\
    ['Y  89','SR 89','RB 89','KR 89','BR 89','ZR 89','NB 89','MO 89','TC 89','RU 89'],\
    ['ZR 90','Y  90','SR 90','RB 90','KR 90','BR 90','NB 90','MO 90','TC 90','RU 90'],\
    ['ZR 91','Y  91','SR 91','RB 91','KR 91','BR 91','SE 91','NB 91','MO 91','TC 91','RU 91','RH 91'],\
    ['ZR 92','Y  92','SR 92','RB 92','KR 92','BR 92','NB 92'],\
    ['ZR 94','Y  94','SR 94','RB 94','KR 94'],\
    ['ZR 96','Y  96','SR 96'],\
    ['NB 93','ZR 93','Y  93','SR 93','RB 93','KR 93','MO 93','TC 93','RU 93','RH 93','PD 93'],\
    ['MO 92','TC 92','RU 92','RH 92','PD 92'],\
    ['MO 94','NB 94','TC 94','RU 94','RH 94','PD 94'],\
    ['MO 95','NB 95','ZR 95','Y  95','SR 95','TC 95','RU 95','RH 95','PD 95','AG 95'],\
    ['MO 96','NB 96','TC 96'],\
    ['MO 97','NB 97','ZR 97','Y  97','SR 97','TC 97','RU 97','RH 97','PD 97','AG 97','CD 97'],\
    ['MO 98','NB 98','ZR 98','Y  98','SR 98'],\
    ['MO100','NB100','ZR100','Y 100'],\
    ['RU 96','RH 96','PD 96','AG 96'],\
    ['RU 98','TC 98','RH 98','PD 98','AG 98','CD 98'],\
    ['RU 99','TC 99','MO 99','NB 99','ZR 99','Y  99','RH 99','PD 99','AG 99','CD 99','IN 99'],\
    ['RU100','TC100','RH100','PD100','AG100','CD100','IN100','SN100'],\
    ['RU101','TC101','MO101','NB101','ZR101','Y 101','RH101','PD101','AG101','CD101','IN101','SN101'],\
    ['RU102','MO102','TC102','NB102','ZR102','Y 102','RH102'],\
    ['RU104','TC104','MO104','NB104'],\
    ['RH103','RU103','TC103','MO103','NB103','ZR103','Y 103','PD103','AG103','CD103','IN103','SN103'],\
    ['PD102','AG102','CD102','IN102','SN102'],\
    ['PD104','RH104','AG104','CD104','IN104','SN104','SB104'],\
    ['PD105','RH105','RU105','TC105','MO105','NB105','ZR105','AG105','CD105','IN105','SN105','SB105'],\
    ['PD106','RH106','RU106','TC106','MO106','NB106','AG106'],\
    ['PD108','RH108','RU108','TC108','MO108','NB108'],\
    ['PD110','RH110','RU110','TC110','MO110','NB110'],\
    ['AG107','PD107','RH107','RU107','TC107','MO107','CD107','IN107','SN107','SB107'],\
    ['AG109','PD109','RH109','RU109','TC109','MO109','NB109','CD109','IN109','SN109','SB109','TE109'],\
    ['CD106','IN106','SN106','SB106'],\
    ['CD108','AG108','IN108','SN108','SB108'],\
    ['CD110','AG110','IN110','SN110','SB110','TE110'],\
    ['CD111','AG111','PD111','RH111','RU111','TC111','IN111','SN111','SB111','TE111','I 111'],\
    ['CD112','AG112','PD112','RH112','RU112','TC112'],\
    ['CD113','AG113','PD113','RH113','RU113'],\
    ['CD114','AG114','PD114','RH114','RU114'],\
    ['CD116','AG116','PD116','RH116'],\
    ['IN113','SN113','SB113','TE113','I 113','XE113'],\
    ['IN115','CD115','AG115','PD115','RH115','RU115'],\
    ['SN112','IN112','SB112','TE112','I 112','XE112'],\
    ['SN114','IN114','SB114','TE114','I 114','XE114','CS114'],\
    ['SN115','SB115','TE115','I 115','XE115','CS115','BA115'],\
    ['SN116','IN116','SB116','TE116','I 116','XE116','CS116','BA116'],\
    ['SN117','IN117','AG117','PD117','RH117','SB117','CD117','TE117','I 117','XE117','CS117','BA117'],\
    ['SN118','IN118','CD118','AG118','PD118','SB118','TE118','XE118','CS118','BA118'],\
    ['SN119','IN119','CD119','AG119','PD119','SB119','TE119','XE119','CS119','BA119'],\
    ['SN120','SB120','IN120','CD120','AG120','PD120'],\
    ['SN122','IN122','CD122','AG122'],\
    ['SN124','IN124','CD124'],\
    ['SB121','SN121','IN121','CD121','AG121','TE121','I 121','XE121','CS121','XE121','CS121','BA121','LA121','CE121'],
    ['SB123','SN123','IN123','CD123','AG123'],\
    ['TE120','I 120','XE120','CS120','BA120','LA120'],\
    ['TE122','SB122','I 122','XE122','CS122','BA122','LA122'],\
    ['TE123','I 123','XE123','CS123','BA123','LA123','CE123'],\
    ['TE124','SB124','I 124'],\
    ['TE125','SB125','SN125','IN125','CD125','I 125','XE125','CS125','BA125','LA125','CE125','PR125'],\
    ['TE126','SB126','SN126','IN126','CD126'],\
    ['TE128','SB128','SN128','IN128','CD128'],\
    ['TE130','SB130','SN130','IN130'],\
    ['I 127','TE127','SB127','SN127','IN127','CD127','XE127','CS127','BA127','LA127','CE127','PR127','ND127'],\
    ['XE124','CS124','BA124','LA124','CE124','PR124'],\
    ['XE126','I 126','CS126','BA126','LA126','CE126','PR126'],\
    ['XE128','I 128','CS128','BA128','LA128','CE128','PR128'],\
    ['XE129','I 129','TE129','SB129','SN129','IN129','CD129','CS129','BA129','LA129','CE129','PR129','ND129'],\
    ['XE130','I 130','CS130'],\
    ['XE131','I 131','TE131','SB131','SN131','IN131','CS131','BA131','LA131','CE131','PR131','ND131','PM131','SM131'],\
    ['XE132','I 132','TE132','SB132','SN132','IN132','CS132'],
    ['XE134','I 134','TE134','SB134','SN134'],\
    ['XE136','I 136','TE136','SB136'],\
    ['CS133','XE133','I 133','TE133','SB133','SN133','BA133','LA133','CE133','PR133','ND133','PM133','SM133'],\
    ['BA130','LA130','CE130','PR130','ND130','PM130'],\
    ['BA132','LA132','CE132','PR132','ND132','PM132','SM132'],\
    ['BA134','CS134','LA134','CE134','PR134','ND134','PM134','SM134','EU134'],\
    ['BA135','CS135','XE135','I 135','TE135','SB135','SN135','LA135','CE135','PR135','ND135','PM135','SM135','EU135','GD135'],\
    ['BA136','CS136','LA136'],\
    ['BA137','CS137','XE137','I 137','TE137','SB137','LA137','CE137','PR137','ND137','PM137','SM137','EU137','GD137'],\
    ['BA138','CS138','XE138','I 138','TE138'],\
    ['LA138'],\
    ['LA139','BA139','CS139','XE139','I 139','CE139','PR139','ND139','PM139','SM139','EU139','GD139','TB139','DY139'],\
    ['CE136','PR136','ND136','PM136','SM136','EU136'],\
    ['CE138','PR138','ND138','PM138','SM138','EU138','GD138'],\
    ['CE140','LA140','BA140','CS140','XE140','I 140','PR140','ND140','PM140','SM140','EU140','GD140','TB140'],\
    ['CE142','LA142','BA142','CS142','XE142','I 142'],\
    ['PR141','CE141','LA141','BA141','CS141','XE141','I 141','ND141','PM141','SM141','EU141','GD141','TB141','DY141'],\
    ['ND142','PR142','PM142','SM142','EU142','GD142','TB142','DY142','HO142','SM146','EU146','GD146','TB146','DY146','HO146','ER146',\
'GD150','TB150','DY150','HO150','ER150','DY154','HO154','ER154'],\
    ['ND143','PR143','CE143','LA143','BA143','CS143','XE143','PM143','SM143','EU143','GD143','TB143','DY143'],\
    ['ND144','PR144','CE144','LA144','BA144','CS144','XE144','PM144'],\
    ['ND145','PR145','CE145','LA145','BA145','CS145','PM145','SM145','EU145','GD145','TB145','DY145','HO145','ER145'],\
    ['ND146','PR146','CE146','LA146','BA146','CS146','PM146'],\
    ['ND148','PR148','CE148','LA148','BA148'],\
    ['ND150','PR150','CE150','LA150','BA150'],\
    ['SM144','EU144','GD144','TB144','DY144','HO144','GD148','TB148','DY148','HO148','ER148','TM148'],\
    ['SM147','PM147','ND147','PR147','CE147','LA147','BA147','EU147','GD147','TB147','DY147','HO147','ER147'],\
    ['SM148','PM148','EU148'],\
    ['SM149','PM149','ND149','PR149','CE149','LA149','EU149','GD149','TB149','DY149','HO149','ER149','TM149','YB149'],\
    ['SM150','PM150','EU150'],\
    ['SM152','PM152','ND152','PR152','CE152','EU152'],\
    ['SM154','PM154','ND154','PR154'],\
    ['EU151','SM151','PM151','ND151','PR151','CE151','GD151','TB151'],\
    ['EU153','SM153','PM153','GD153','PR153','GD153','TB153','DY153','HO153'],\
    ['GD152','TB152','DY152','HO152'],\
    ['GD154','EU154','TB154'],\
    ['GD155','EU155','SM155','PM155','ND155','TB155','DY155','HO155','ER155','TM155'],\
    ['GD156','EU156','SM156','PM156','ND156','TB156'],\
    ['GD157','EU157','SM157','PM157','TB157','DY157','HO157','ER157','TM157','YB157'],\
    ['GD158','EU158','SM158','PM158','TB158'],\
    ['GD160','EU160','SM160'],\
    ['TB159','GD159','EU159','SM159','PM159','DY159','HO159','ER159','TM159','YB159','LU159'],\
    ['DY156','HO156','ER156','TM156','YB156','HO156','ER156','TM156','YB156'],\
    ['DY158','HO158','ER158','TM158','YB158','LU158'],\
    ['DY160','TB160','HO160','ER160','TM160','YB160','LU160','HF160'],\
    ['DY161','TB161','GD161','EU161','SM161','PM161','ND161','HO161','ER161','TM161','YB161','LU161','HF161','TA161'],\
    ['DY162','TB162','GD162','EU162','SM162','PM162','HO162'],\
    ['DY163','TB163','HO163','GD163','EU163','SM163','PM163','ND163','ER163','TM163','YB163','LU163','HF163','TA163'],\
    ['DY164','TB164','HO164','GD164','EU164','SM164','PM164','ND164','HO164'],\
    ['HO165','DY165','ER165','TB165','GD165','EU165','SM165','PM165','HO165','TM165','YB165','LU165','HF165','TA165','W 165'],\
    ['ER162','TM162','YB162','LU162','HF162','TA162'],\
    ['ER164','TM164','YB164','LU164','HF164','TA164','W 164'],\
    ['ER166','HO166','DY166','TB166','GD166','EU166','SM166','PM166','ND166','TM166','YB166','LU166','HF166','TA166','W 166','RE166'],\
    ['ER167','HO167','DY167','TB167','GD167','EU167','SM167','PM167','ND167','TM167','YB167','LU167','HF167','TA167','W 167','RE167'],\
    ['ER168','HO168','DY168','TB168','GD168','EU168','SM168','PM168','ND168','TM168'],\
    ['ER170','HO170','DY170','TB170','GD170','EU170','SM170','PM170','ND170'],\
    ['TM169','ER169','HO169','DY169','TB169','GD169','EU169','SM169','PM169','ND169','YB169','LU169','HF169','TA169','W 169','RE169'],\
    ['YB168','LU168','HF168','TA168','W 168','RE168'],\
    ['YB170','TM170','LU170','HF170','TA170','W 170','RE170','OS170'],\
    ['YB171','TM171','ER171','HO171','DY171','TB171','GD171','EU171','SM171','PM171','ND171','LU171','HF171','TA171','W 171','RE171','OS171'],\
    ['YB172','TM172','ER172','HO172','DY172','TB172','GD172','EU172','SM172','PM172','ND172','LU172','HF172','TA172','W 172','RE172','OS172','IR172'],\
    ['YB173','TM173','ER173','HO173','DY173','TB173','GD173','EU173','SM173','PM173','ND173','LU173','HF173','TA173','W 173','RE173','OS173','IR173'],\
    ['YB174','TM174','ER174','HO174','DY174','TB174','GD174','EU174','SM174','PM174','ND174','LU174'],\
    ['YB176','TM176','ER176','HO176','DY176','TB176','GD176','EU176','SM176','PM176','ND176'],\
    ['LU175','YB175','TM175','ER175','HO175','DY175','TB175','GD175','EU175','SM175','PM175','ND175','HF175','TA175','W 175','RE175','OS175','IR175'],\
    ['LU176','HF176','TA176','W 176','RE176','OS176','IR176'],\
    ['HF174','TA174','W 174','RE174','OS174','IR174'],\
    ['HF176','TA176','W 176','RE176','OS176','IR176'],\
    ['HF177','LU177','YB177','LU177','YB177','TM177','ER177','HO177','DY177','TB177','GD177','TA177','W 177','RE177','OS177','IR177'],\
    ['HF178','LU178','YB178','LU178','YB178','TM178','ER178','HO178','DY178','TB178','GD178','TA178','W 178','RE178','OS178','IR178'],\
    ['HF179','LU179','YB179','LU179','YB179','TM179','ER179','HO179','DY179','TB179','GD179','TA179','W 179','RE179','OS179','IR179','PT179'],\
    ['HF180','LU180','YB180','LU180','YB180','TM180','ER180','HO180','DY180','TB180','GD180','TA180','W 180','RE180','OS180','IR180','PT180','AU180'],\
    ['TA180'],\
    ['TA181','HF181','LU181','YB181','LU181','YB181','TM181','ER181','HO181','DY181','TB181','GD181','W 181','RE181','OS181','IR181','PT181','AU181'],\
    ['W 180','RE180','OS180','IR180','PT180','AU180'],\
    ['W 182','TA182','HF182','LU182','YB182','LU182','YB182','TM182','ER182','HO182','DY182','TB182','GD182','RE182','OS182','IR182','PT182','AU182'],\
    ['W 183','TA183','HF183','LU183','YB183','LU183','YB183','TM183','ER183','HO183','DY183','TB183','GD183','RE183','OS183','IR183','PT183','AU183'],\
    ['W 184','TA184','HF184','LU184','YB184','LU184','YB184','TM184','ER184','HO184','DY184','TB184','GD184','RE184'],\
    ['W 186','TA186','HF186','LU186','YB186','LU186','YB186','TM186','ER186','HO186','DY186','TB186','GD186'],\
    ['RE185','W 185','TA185','HF185','LU185','YB185','LU185','YB185','TM185','ER185','HO185','DY185','TB185','GD185','OS185','IR185','PT185','AU185','HG185','TL185'],\
    ['RE187','W 187','TA187','HF187','LU187','YB187','LU187','YB187','TM187','ER187','HO187','DY187','TB187','GD187'],\
    ['OS184','IR184','PT184','AU184','HG184','TL184'],\
    ['OS186','RE186','IR186','PT186','AU186','HG186','TL186'],\
    ['OS187','IR187','PT187','AU187','HG187','TL187','PB187'],\
    ['OS188','RE188','W 188','TA188','HF188','LU188','YB188','LU188','YB188','TM188','ER188','HO188','DY188','TB188','GD188','IR188','PT188','AU188','HG188','TL188','PB188'],\
    ['OS189','RE189','W 189','TA189','HF189','LU189','YB189','LU189','YB189','TM189','ER189','HO189','DY189','TB189','GD189','IR189','PT189','AU189','HG189','TL189','PB189'],\
    ['OS190','RE190','W 190','TA190','HF190','LU190','YB190','LU190','YB190','TM190','ER190','HO190','DY190','TB190','GD190','IR190'],\
    ['OS192','RE190','W 192','TA192','HF192','LU192','YB192','LU192','YB192','TM192','ER192','HO192','DY192','TB192','GD192'],\
    ['IR191','OS191','RE191','W 191','TA191','HF191','LU191','YB191','LU191','YB191','TM191','ER191','HO191','DY191','TB191','GD191','PT191','AU191','HG191','TL191','PB191'],\
    ['IR193','OS193','RE193','W 193','TA193','HF193','LU193','YB193','LU193','YB193','TM193','ER193','HO193','DY193','TB193','GD193','PT193','AU193','HG193','TL193','PB193'],\
    ['PT190','AU190','HG190','TL190','PB190'],\
    ['PT192','IR192','AU192','HG192','TL192','PB192'],\
    ['PT194','IR194','OS194','RE194','W 194','TA194','HF194','AU194','HG194','TL194','PB194','BI194'],
    ['PT195','IR195','OS195','RE195','W 195','TA195','HF195','AU195','HG195','TL195','PB195','BI195'],\
    ['PT196','IR196','OS196','RE196','W 196','TA196','HF196','AU196'],\
    ['PT198','IR198','OS198','RE198','W 198','TA198','HF198'],\
    ['AU197','PT197','IR197','OS197','RE197','W 197','TA197','HF197','HG197','TL197','PB197','BI197'],\
    ['HG196','TL196','PB196','BI196'],\
    ['HG198','AU198','TL198','PB198','BI198'],\
    ['HG199','AU199','PT199','IR199','OS199','RE199','W 199','TA199','HF199','TL199','PB199','BI199'],\
    ['HG200','AU200','PT200','IR200','OS200','RE200','W 200','TA200','HF200','TL200','PB200','BI200'],\
    ['HG201','AU201','PT201','IR201','OS201','RE201','W 201','TA201','HF201','TL201','PB201','BI201','PO201'],\
    ['HG202','AU202','PT202','IR202','OS202','RE202','W 202','TA202','HF202','TL202','PB202','BI202','PO202'],\
    ['HG204','AU204','PT204','IR204','OS204','RE204','W 204','TA204','HF204'],\
    ['TL203','HG203','AU203','PT203','IR203','OS203','RE203','W 203','TA203','HF203','PB203','BI203','PO203'],\
    ['TL205','HG205','AU205','PT205','IR205','OS205','RE205','W 205','TA205','HF205','PB205','BI205','PO205'],\
    ['PB204','TL204','BI204','PO204'],\
    ['PB206','TL206','HG206','AU206','PT206','IR206','OS206','RE206','W 206','TA206','HF206','BI206','PO210'],\
    ['PB207','TL207','HG207','AU207','PT207','IR207','OS207','RE207','W 207','TA207','HF207','BI207','PO211','BI211'],\
    ['PB208','TL208','HG208','AU208','PT208','IR208','OS208','RE208','W 208','TA208','HF208','BI208','PO212','BI212'],\
    ['BI209','PB209','TL209','HG209','AU209','PT209','IR209','OS209','RE209','W 209','TA209','HF209'],\
    ['TH232','AC232','RA232','FR232','RN232'],\
    ['U 235','PA235','TH235','AC235','RA235','FR235','RN235'],\
    ['U 238','PA238','TH238','AC238','RA238','FR238','RN238']]
    #print decay_raw

def give_zip_element_z_and_names(element_name):
    ''' create 2 indexes that, given the name of the element/specie, give the atomic number.'''

    #import numpy as np

    global z_bismuth
    z_bismuth = 83
    global z_for_elem
    z_for_elem = []
    global index_stable
    index_stable = []

    i_for_stable = 1
    i_for_unstable = 0
    for i in range(z_bismuth):
        z_for_elem.append(int(i+1))
        # the only elements below bismuth with no stable isotopes are Tc and Pm
        if i+1 == 43 or i+1 == 61:
            index_stable.append(i_for_unstable)
        else:
            index_stable.append(i_for_stable)

        dummy_index = np.zeros(len(element_name))
        for i in range(len(element_name)):
            if element_name[i] == 'Neutron':
                dummy_index[i] = 0
            elif element_name[i] == 'H':
                dummy_index[i] = 1
            elif element_name[i] == 'He':
                dummy_index[i] = 2
            elif element_name[i] == 'Li':
                dummy_index[i] = 3
            elif element_name[i] == 'Be':
                dummy_index[i] = 4
            elif element_name[i] == 'B':
                dummy_index[i] = 5
            elif element_name[i] == 'C':
                dummy_index[i] = 6
            elif element_name[i] == 'N':
                dummy_index[i] = 7
            elif element_name[i] == 'O':
                dummy_index[i] = 8
            elif element_name[i] == 'F':
                dummy_index[i] = 9
            elif element_name[i] == 'Ne':
                dummy_index[i] = 10
            elif element_name[i] == 'Na':
                dummy_index[i] = 11
            elif element_name[i] == 'Mg':
                dummy_index[i] = 12
            elif element_name[i] == 'Al':
                dummy_index[i] = 13
            elif element_name[i] == 'Si':
                dummy_index[i] = 14
            elif element_name[i] == 'P':
                dummy_index[i] = 15
            elif element_name[i] == 'S':
                dummy_index[i] = 16
            elif element_name[i] == 'Cl':
                dummy_index[i] = 17
            elif element_name[i] == 'Ar':
                dummy_index[i] = 18
            elif element_name[i] == 'K':
                dummy_index[i] = 19
            elif element_name[i] == 'Ca':
                dummy_index[i] = 20
            elif element_name[i] == 'Sc':
                dummy_index[i] = 21
            elif element_name[i] == 'Ti':
                dummy_index[i] = 22
            elif element_name[i] == 'V':
                dummy_index[i] = 23
            elif element_name[i] == 'Cr':
                dummy_index[i] = 24
            elif element_name[i] == 'Mn':
                dummy_index[i] = 25
            elif element_name[i] == 'Fe':
                dummy_index[i] = 26
            elif element_name[i] == 'Co':
                dummy_index[i] = 27
            elif element_name[i] == 'Ni':
                dummy_index[i] = 28
            elif element_name[i] == 'Cu':
                dummy_index[i] = 29
            elif element_name[i] == 'Zn':
                dummy_index[i] = 30
            elif element_name[i] == 'Ga':
                dummy_index[i] = 31
            elif element_name[i] == 'Ge':
                dummy_index[i] = 32
            elif element_name[i] == 'As':
                dummy_index[i] = 33
            elif element_name[i] == 'Se':
                dummy_index[i] = 34
            elif element_name[i] == 'Br':
                dummy_index[i] = 35
            elif element_name[i] == 'Kr':
                dummy_index[i] = 36
            elif element_name[i] == 'Rb':
                dummy_index[i] = 37
            elif element_name[i] == 'Sr':
                dummy_index[i] = 38
            elif element_name[i] == 'Y':
                dummy_index[i] = 39
            elif element_name[i] == 'Zr':
                dummy_index[i] = 40
            elif element_name[i] == 'Nb':
                dummy_index[i] = 41
            elif element_name[i] == 'Mo':
                dummy_index[i] = 42
            elif element_name[i] == 'Tc':
                dummy_index[i] = 43
            elif element_name[i] == 'Ru':
                dummy_index[i] = 44
            elif element_name[i] == 'Rh':
                dummy_index[i] = 45
            elif element_name[i] == 'Pd':
                dummy_index[i] = 46
            elif element_name[i] == 'Ag':
                dummy_index[i] = 47
            elif element_name[i] == 'Cd':
                dummy_index[i] = 48
            elif element_name[i] == 'In':
                dummy_index[i] = 49
            elif element_name[i] == 'Sn':
                dummy_index[i] = 50
            elif element_name[i] == 'Sb':
                dummy_index[i] = 51
            elif element_name[i] == 'Te':
                dummy_index[i] = 52
            elif element_name[i] == 'I':
                dummy_index[i] = 53
            elif element_name[i] == 'Xe':
                dummy_index[i] = 54
            elif element_name[i] == 'Cs':
                dummy_index[i] = 55
            elif element_name[i] == 'Ba':
                dummy_index[i] = 56
            elif element_name[i] == 'La':
                dummy_index[i] = 57
            elif element_name[i] == 'Ce':
                dummy_index[i] = 58
            elif element_name[i] == 'Pr':
                dummy_index[i] = 59
            elif element_name[i] == 'Nd':
                dummy_index[i] = 60
            elif element_name[i] == 'Pm':
                dummy_index[i] = 61
            elif element_name[i] == 'Sm':
                dummy_index[i] = 62
            elif element_name[i] == 'Eu':
                dummy_index[i] = 63
            elif element_name[i] == 'Gd':
                dummy_index[i] = 64
            elif element_name[i] == 'Tb':
                dummy_index[i] = 65
            elif element_name[i] == 'Dy':
                dummy_index[i] = 66
            elif element_name[i] == 'Ho':
                dummy_index[i] = 67
            elif element_name[i] == 'Er':
                dummy_index[i] = 68
            elif element_name[i] == 'Tm':
                dummy_index[i] = 69
            elif element_name[i] == 'Yb':
                dummy_index[i] = 70
            elif element_name[i] == 'Lu':
                dummy_index[i] = 71
            elif element_name[i] == 'Hf':
                dummy_index[i] = 72
            elif element_name[i] == 'Ta':
                dummy_index[i] = 73
            elif element_name[i] == 'W':
                dummy_index[i] = 74
            elif element_name[i] == 'Re':
                dummy_index[i] = 75
            elif element_name[i] == 'Os':
                dummy_index[i] = 76
            elif element_name[i] == 'Ir':
                dummy_index[i] = 77
            elif element_name[i] == 'Pt':
                dummy_index[i] = 78
            elif element_name[i] == 'Au':
                dummy_index[i] = 79
            elif element_name[i] == 'Hg':
                dummy_index[i] = 80
            elif element_name[i] == 'Tl':
                dummy_index[i] = 81
            elif element_name[i] == 'Pb':
                dummy_index[i] = 82
            elif element_name[i] == 'Bi':
                dummy_index[i] = 83
            elif element_name[i] == 'Po':
                dummy_index[i] = 84
            elif element_name[i] == 'At':
                dummy_index[i] = 85

        #if spe[0] == 'N   1':
        #       znum_int[0] = 0

        # here the index to connect name and atomic numbers.
        global index_z_for_elements
        index_z_for_elements = {}
        for a,b in zip(element_name,dummy_index):
            index_z_for_elements[a]=b


def get_z_from_el(element):
    '''
    Very simple function that gives the atomic number AS A STRING when given the element symbol.
    Uses predefined a dictionnary.
    Parameter :
    element : string
    For the other way, see get_el_from_z
    '''
    dict_name={'Ru': '44', 'Re': '75', 'Ra': '88', 'Rb': '37', 'Rn': '86', 'Rh': '45', 'Be': '4', 'Ba': '56', 'Bi': '83', 'Br': '35', 'H': '1', 'P': '15', 'Os': '76', 'Hg': '80', 'Ge': '32', 'Gd': '64', 'Ga': '31', 'Pr': '59', 'Pt': '78', 'C': '6', 'Pb': '82', 'Pa': '91', 'Pd': '46', 'Cd': '48', 'Po': '84', 'Pm': '61', 'Ho': '67', 'Hf': '72', 'K': '19', 'He': '2', 'Mg': '12', 'Mo': '42', 'Mn': '25', 'O': '8', 'S': '16', 'W': '74', 'Zn': '30', 'Eu': '63', 'Zr': '40', 'Er': '68', 'Ni': '28', 'Na': '11', 'Nb': '41', 'Nd': '60', 'Ne': '10', 'Fr': '87', 'Fe': '26', 'B': '5', 'F': '9', 'Sr': '38', 'N': '7', 'Kr': '36', 'Si': '14', 'Sn': '50', 'Sm': '62', 'V': '23', 'Sc': '21', 'Sb': '51', 'Se': '34', 'Co': '27', 'Cl': '17', 'Ca': '20', 'Ce': '58', 'Xe': '54', 'Lu': '71', 'Cs': '55', 'Cr': '24', 'Cu': '29', 'La': '57', 'Li': '3', 'Tl': '81', 'Tm': '69', 'Th': '90', 'Ti': '22', 'Te': '52', 'Tb': '65', 'Tc': '43', 'Ta': '73', 'Yb': '70', 'Dy': '66', 'I': '53', 'U': '92', 'Y': '39', 'Ac': '89', 'Ag': '47', 'Ir': '77', 'Al': '13', 'As': '33', 'Ar': '18', 'Au': '79', 'At': '85', 'In': '49'}
    return int(dict_name[element])

def get_el_from_z(z):
    '''
    Very simple Vfunction that gives the atomic number AS A STRING when given the element symbol.
    Uses predefined a dictionnary.
    Parameter :
    z : string or number
    For the other way, see get_z_from_el
    '''
    if(type(z)==float):
        z=int(z)
    if(type(z)==int):
        z=str(z)
    dict_z={'24': 'Cr', '25': 'Mn', '26': 'Fe', '27': 'Co', '20': 'Ca', '21': 'Sc', '22': 'Ti', '23': 'V', '28': 'Ni', '29': 'Cu', '4': 'Be', '8': 'O', '59': 'Pr', '58': 'Ce', '55': 'Cs', '54': 'Xe', '57': 'La', '56': 'Ba', '51': 'Sb', '50': 'Sn', '53': 'I', '52': 'Te', '88': 'Ra', '89': 'Ac', '82': 'Pb', '83': 'Bi', '80': 'Hg', '81': 'Tl', '86': 'Rn', '87': 'Fr', '84': 'Po', '85': 'At', '3': 'Li', '7': 'N', '39': 'Y', '38': 'Sr', '33': 'As', '32': 'Ge', '31': 'Ga', '30': 'Zn', '37': 'Rb', '36': 'Kr', '35': 'Br', '34': 'Se', '60': 'Nd', '61': 'Pm', '62': 'Sm', '63': 'Eu', '64': 'Gd', '65': 'Tb', '66': 'Dy', '67': 'Ho', '68': 'Er', '69': 'Tm', '2': 'He', '6': 'C', '91': 'Pa', '90': 'Th', '92': 'U', '11': 'Na', '10': 'Ne', '13': 'Al', '12': 'Mg', '15': 'P', '14': 'Si', '17': 'Cl', '16': 'S', '19': 'K', '18': 'Ar', '48': 'Cd', '49': 'In', '46': 'Pd', '47': 'Ag', '44': 'Ru', '45': 'Rh', '42': 'Mo', '43': 'Tc', '40': 'Zr', '41': 'Nb', '1': 'H', '5': 'B', '9': 'F', '77': 'Ir', '76': 'Os', '75': 'Re', '74': 'W', '73': 'Ta', '72': 'Hf', '71': 'Lu', '70': 'Yb', '79': 'Au', '78': 'Pt'}
    return dict_z[z]
