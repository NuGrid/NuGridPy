
#
# NuGridpy - Tools for accessing and visualising NuGrid data.
#
# Copyright 2007 - 2014 by the NuGrid Team.
# All rights reserved. See LICENSE.
#

"""
Ascii_table.py: read and write simple ascii tables

By Daniel Alexander Bertolino Conti
Fall 2010
If the user find any bugs or errors, please email fherwig@uvic.ca.

Assumptions for ascii Files:

    Headers are always at the beginning of the file and start with a
    capital H (default, can be reset).  The next line after the header
    lines, is a line of column attribute names.  The column attribute
    names are separated by ' ' by default or whatever the user
    dictates to the class.  Data columns are seperated by spaces and
    each data attribute contains no spaces.  All the data columns are
    of equal length.  Any file name that has 'trajectory' or
    'Trajectory' in it, is assumed to be a trajectory type file.

Assumptions for Trajectory Files:

    The first three lines start with a '#'.  The first of these contains
    'time', 'T' and 'rho', each seperated by a space.  The second is
    '# YRS/SEC; T8K/T9K; CGS/LOG' which are the choices for the header
    attributes AGEUNIT, TUNIT, and RHOUNIT can be.  The third is
    something like 'FORMAT:' The four lines afther those three are the
    header attributes.  Each line has one header attribute, followd by
    a = and then followed after by the value.  After these seven lines
    comes three columns of data, the first is is associated with
    'time', the second with 'T' and the third with rho.

"""
from __future__ import print_function
from __future__ import absolute_import

from builtins import input
from builtins import str
from builtins import range
from numpy import *
from .data_plot import *
import matplotlib.pylab as pyl
import matplotlib.pyplot as pl
import os

class ascii_table(DataPlot):
    '''
    Data structure to read simple data tables and trajectory data
    tables.

    Init method that reads in ascii type files and trajectory type
    files.  By default this method reads ascii type files.  If the user
    wants a trajectory file read, either the file must have 'trajectory'
    in the filename or the user must set the datatype='trajectory'.

    Parameters
    ----------
    filename : string
        The name of the file we are looking at, or writeng to.
    sldir : string, optional
        Standard directory of filename. The default is '.'.
    sep : string, optional
        The seperator that seperates column attributes in filename.  The
        defaults is '  '.
    datatype : string, optional
        What type of ascii table you are reading and or writeing.  The
        only two options currently are 'normal' and 'trajectory'.  The
        default is 'normal'.
    Headers : list, optional
        A list of Header strings or if the file being writen is of type
        trajectory, this is a dictionary of header attributes and their
        associated values. The default is [].
    dcols : list, optional
        A list of data attributes. The default is [].
    data : list, optional
        A list of lists (or of numpy arrays) of columns of data.  The
        default is [].
    header_char : caracter, optional
        The character that indicates a header lines. The default is 'H'.
    read : boolean, optional
        Boolean of weather this is reading or writing a file.  The
        default is True.
    headerlines : list, optional
        Additional list of strings of header data, only used in
        trajectory data Types. The default is [].

    Examples
    --------
    >>> from ascii_table import *
    >>> p=ascii_table('c12cg.dat')
    >>> p.hattrs
    ['1 12  6  1  1  1  0  0  0  1 13  7  0', '55   1.943', 'c12pg']
    >>> a.dcols
    ['upper', 'lower', 'T9', 'ado(gs)', 'ado/CA88', 'tt/gs', 'ado/fit', 'CA88', 'fitted']
    >>> a.get('upper')
    [1.3400000000000001e-24, ... 1590.0]
    >>> a.plot('T9', 'ado/CA88') #plots data

    Also there is a global method called write that will allow a user
    to write ascii Files and trajectory files.

    One Can call this method like

    >>> write('file.txt',dcols,data,headers)

    Where Dcols is a list of data attributes, data is a list of lists
    of data and headers is a list of strings that are each a header
    attribute.

    See Also
    --------
    write()

    '''

    def __init__(self, filename, sldir='.', sep='  ', datatype='normal',
                 headers=[], dcols=[], data=[], header_char='H',
                 read=True, headerlines=[]):
        '''
        Init method that reads in ascii type files and trajectory type
        files.  By default this method reads ascii type files.  If the
        user wants a trajectory file read, either the file must have
        'trajectory' in the filename or the user must set the
        datatype='trajectory'.

        Parameters
        ----------
        filename : string
            The name of the file we are looking at, or writeing to.
        sldir : string, optional
            Standard directory of filename. The default is '.'.
        sep : string, optional
            The seperator that seperates column attributes in filename.
            The defaults is '  '.
        datatype : string, optional
            What type of ascii table you are reading and or writeing.
            The only two options currently are 'normal' and
            'trajectory'.  The default is 'normal'.
        Headers : list, optional
            A list of Header strings or if the file being writen is of
            type trajectory, this is a dictionary of header attributes
            and their associated values. The default is [].
        dcols : list, optional
            A list of data attributes. The default is [].
        data : list, optional
            A list of lists (or of numpy arrays) of columns of data.
            The default is [].
        header_char : caracter, optional
            The character that indicates a header lines.  The default
            is 'H'.
        read : boolean, optional
            Boolean of weather this is reading or writing a file.  The
            default is True.
        headerlines : list, optional
            Additional list of strings of header data, only used in
            trajectory data Types. The default is [].
        '''
        self.sldir=sldir
        self.files = []
        self.header_char=header_char
        self.hattrs=[]
        self.dcols=[]
        self.data ={}
        self.headerLines=[]
        self.files.append(filename)
        self.datatype=datatype
        if 'trajectory' in filename or 'Trajectory' in filename:
            self.datatype='trajectory'

        if read:

            self.hattrs,self.data=self._readFile(sldir,filename,sep)
            self.dcols=list(self.data.keys())
        else:
            a=self.write(filename,headers,dcols,data,headerlines,sldir,sep)
            if a ==None:
                return None
            self.hattrs,self.data=self._readFile(sldir,filename,sep)
        self.dcols=list(self.data.keys())

    def get(self, attri):
        '''
        Method that dynamically determines the type of attribute that is
        passed into this method. Also it then returns that attribute's
        associated data.

        Parameters
        ----------
        attri : string
            The attribute we are looking for.

        '''
        isCol=False
        isHead=False

        if attri in self.dcols:
            isCol=True
        elif attri in self.hattrs:
            isHead=True
        else:
            print("That attribute does not exist in this File")
            print('Returning None')

        if isCol:
            return self.getColData(attri)
        elif isHead:
            return hattrs

    def getColData(self, attri):
        '''
        Method that returns column data

        Parameters
        ----------
        attri : string
            The attribute we are looking for.

        '''
        return self.data[attri]

    def _readFile(self, sldir, fileName, sep):
        '''
        Private method that reads in the header and column data.

        '''

        if sldir.endswith(os.sep):
            fileName = str(sldir)+str(fileName)
        else:
            fileName = str(sldir)+os.sep+str(fileName)


        fileLines=[] #list of lines in the file
        header=[]    #list of Header lines
        dataCols=[]  #Dictionary of data column names
        data=[]      #List of Data lists
        cols=[]      #List of column names

        f=open(fileName,'r')
        fileLines=f.readlines()
        i=0
        if self.datatype != 'trajectory':

            while i<len(fileLines):
                if fileLines[i].startswith(self.header_char):
                    tmp=fileLines[i].lstrip(self.header_char)
                    header.append(tmp.strip())
                else:
                    break
                i+=1

            cols=fileLines[i].split(sep)

            tmp=[]
            tmp1=[]
            for j in range(len(cols)):
                tmp1=cols[j].strip()
                if tmp1 !='':
                    tmp.append(tmp1)
            cols=tmp
            i+=1
        else:
            header={}
            while fileLines[i].startswith('#') or '=' in fileLines[i]:
                if fileLines[i].startswith('#') and cols==[]:
                    cols=fileLines[i].strip('#')
                    cols=cols.strip()
                    cols=cols.split()
                elif fileLines[i].startswith('#'):
                    tmp1=fileLines[i].strip('#')
                    tmp1=tmp1.strip()
                    self.headerLines.append(tmp1)
                elif not fileLines[i].startswith('#'):
                    tmp=fileLines[i].split('=')
                    tmp[0]=tmp[0].strip()
                    tmp[1]=tmp[1].strip()
                    if header=={}:
                        header={str(tmp[0]):str(tmp[1])}
                    else:
                        header[str(tmp[0])]=str(tmp[1])
                i+=1
        while i<len(fileLines):
            tmp=fileLines[i].split()
            for j in range(len(tmp)):
                tmp[j]=tmp[j].strip()
            data.append(tmp)
            i+=1
        tmp=[]
        tmp1=[]
        for j in range(len(data)):
            for k in range(len(data[j])):
                tmp1=data[j][k].strip()
                if tmp1 !='':
                    tmp.append(tmp1)
            data[j]=tmp
            tmp=[]
        tmp=[]

        for j in range(len(cols)):
            for k in range(len(data)):
                tmp.append(float(data[k][j]))
            tmp=array(tmp)

            if j == 0:
                dataCols={cols[j]:tmp}
            else:
                dataCols[cols[j]]=tmp
            tmp=[]

        return header,dataCols

#Global methods
def writeTraj(filename='trajectory.input', data=[], ageunit=0, tunit=0,
              rhounit=0, idNum=0):
    '''
    Method for writeing Trajectory type ascii files files.

    Parameters
    ----------
    filename : string
        The file where this data will be written.
    data : list
        A list of 1D data vectors with time, T and rho.
    ageunit : integer, optional
        If 1 ageunit = SEC, If 0 ageunit = YRS. If 2 agunit =
        logtimerev in yrs.  The default is 0.  logtimerev is log of
        time until end
    tunit : integer, optional
        If 1 TUNIT  = T9K, if 0 TUNIT = T8K.  The default is 0.
    rhounit : integer, optional
        If 1 RHOUNIT = LOG, if 0 RHOUNIT = CGS. The default is 0.
    idNum : optional
        An optional id argument

    '''
    if data==[]:
        print('Please input correct data')
        print('returning None')
        return None
    headers=[]
    if ageunit ==1:
        headers.append('AGEUNIT = SEC')
    elif ageunit==0:
        headers.append('AGEUNIT = YRS')
    elif ageunit==2:
        headers.append('AGEUNIT = logtimerev/yrs')

    if tunit ==1:
        headers.append('TUNIT   = T9K')
    elif tunit==0:
        headers.append('TUNIT   = T8K')

    if rhounit ==1:
        headers.append('RHOUNIT = LOG')
    elif rhounit==0:
        headers.append('RHOUNIT = CGS')

    headers.append('ID = '+str(idNum))

    write(filename,headers,['time','T','rho'],data,['YRS/SEC; T8K/T9K; CGS/LOG',"FORMAT: '(10x,A3)'"],trajectory=True)

def write(filename, headers, dcols, data, headerlines=[],
          header_char='H', sldir='.', sep='  ', trajectory=False,
          download=False):
    '''
    Method for writeing Ascii files.

    Note the attribute name at position i in dcols will be associated
    with the column data at index i in data.  Also the number of data
    columns(in data) must equal the number of data attributes (in dcols)
    Also all the lengths of that columns must all be the same.

    Parameters
    ----------
    filename : string
        The file where this data will be written.
    Headers : list
        A list of Header strings or if the file being written is of
        type trajectory, this is a List of strings that contain header
        attributes and their associated values which are seperated by
        a '='.
    dcols : list
        A list of data attributes.
    data : list
        A list of lists (or of numpy arrays).
    headerlines : list, optional
        Additional list of strings of header data, only used in
        trajectory data Types. The default is [].
    header_char : character, optional
        The character that indicates a header lines. The default is 'H'.
    sldir : string, optional
        Where this fill will be written. The default is '.'.
    sep : string, optional
        What seperates the data column attributes. The default is '  '.
    trajectory : boolean, optional
        Boolean of if we are writeing a trajectory type file. The
        default is False.
    download : boolean, optional
        If using iPython notebook, do you want a download link for
        the file you write?
        The default is False.

    '''
    if sldir.endswith(os.sep):
        filename = str(sldir)+str(filename)
    else:
        filename = str(sldir)+os.sep+str(filename)
    tmp=[] #temp variable
    lines=[]#list of the data lines
    lengthList=[]# list of the longest element (data or column name)
                 # in each column

    if os.path.exists(filename):
        print('Warning this method will overwrite '+ filename)
        print('Would you like to continue? (y)es or (n)no?')
        s = input('--> ')
        if s=='Y' or s=='y' or s=='Yes' or s=='yes':
            print('Yes selected')
            print('Continuing as normal')
        else:
            print('No Selected')
            print('Returning None')
            return None

    if len(data)!=len(dcols):
        print('The number of data columns does not equal the number of Data attributes')
        print('returning none')
        return None
    if trajectory:
        sep=' '
    for i in range(len(headers)):
        if not trajectory:
            tmp.append(header_char+' '+headers[i]+'\n')
        else:
            tmp.append(headers[i]+'\n')
    headers=tmp
    tmp=''

    for i in range(len(data)): #Line length stuff
        length=len(dcols[i])
        for j in range(len(data[i])):
            if len(str(data[i][j]))>length:
                length=len(str(data[i][j]))
        lengthList.append(length)

    tmp=''
    tmp1=''
    if trajectory:
        tmp='#'
    for i in range(len(dcols)):
        tmp1=dcols[i]
        if not trajectory:
            if len(dcols[i]) < lengthList[i]:
                j=lengthList[i]-len(dcols[i])
                for k in range(j):
                    tmp1+=' '
            tmp+=sep+tmp1
        else:
            tmp+=' '+dcols[i]
    tmp+='\n'
    dcols=tmp
    tmp=''


    for i in range(len(data[0])):
        for j in range(len(data)):
            tmp1=str(data[j][i])
            if len(str(data[j][i])) < lengthList[j]:
                l=lengthList[j]-len(str(data[j][i]))
                for k in range(l):
                    tmp1+=' '
            tmp+=sep+tmp1
        lines.append(tmp+'\n')
        tmp=''

    f=open(filename,'w')
    if not trajectory:
        for i in range(len(headers)):
            f.write(headers[i])
        f.write(dcols)
    else:
        f.write(dcols)
        for i in range(len(headerlines)):
            f.write('# '+headerlines[i]+'\n')
        for i in range(len(headers)):
            f.write(headers[i])
    for i in range(len(lines)):
        f.write(lines[i])

    f.close()
    if download:
        from IPython.display import FileLink, FileLinks
        return FileLink(filename)
    else:
        return None

def writeGCE_table(filename,headers,data,dcols=['Isotopes','Yields','Z','A'],header_char='H',sldir='.',sep='&'):
    '''
    Method for writeing data in GCE format in Ascii files.
    Reads either elements or isotopes
    dcols[0] needs to contain either isotopes or elements

    Note the attribute name at position i in dcols will be associated
    with the column data at index i in data.
    Also the number of data columns(in data) must equal the number
    of data attributes (in dcols)
    Also all the lengths of that columns must all be the same.
    Input:
    filename: The file where this data will be written.
    Headers: A list of Header strings or if the file being written
         is of type trajectory, this is a List of strings
         that contain header attributes and their associated
         values which are seperated by a '='.
    dcols: A list of data attributes
    data:  A list of lists (or of numpy arrays).
    header_char  the character that indicates a header lines
    sldir: Where this fill will be written.
    sep: What seperatesa the data column attributes
    trajectory: Boolean of if we are writeing a trajectory type file
    '''

    import re
    from . import utils as u

    #check if input are elements or isotopes
    if not '-' in data[0][0]:
        iso_inp=False
        dcols=dcols+['Z']
    else:
        iso_inp=True
        dcols=dcols+['Z','A']
    #Attach Z and A
    if iso_inp==True:
        data.append([])
        data.append([])
        u.convert_specie_naming_from_h5_to_ppn(data[0])
        Z=u.znum_int
        A=u.amass_int
        for i in range(len(data[0])):
            zz=str(int(Z[i]))
            aa=str(int(A[i]))
            data[1][i]='{:.3E}'.format(data[1][i])+' '
            data[-2].append(zz)
            data[-1].append(aa)


    else:
        #in order to get Z , create fake isotope from element
        fake_iso=[]
        for k in range(len(data[0])):
            fake_iso.append(data[0][k]+'-99')
        #print fake_iso
            data.append([])
            u.convert_specie_naming_from_h5_to_ppn(fake_iso)
            Z=u.znum_int
        for i in range(len(data[0])):
            zz=str(int(Z[i]))
            data[1][i]='{:.3E}'.format(data[1][i])+' '
            data[-1].append(zz)


    if sldir.endswith(os.sep):
        filename = str(sldir)+str(filename)
    else:
        filename = str(sldir)+os.sep+str(filename)
    tmp=[] #temp variable
    lines=[]#list of the data lines
    lengthList=[]# list of the longest element (data or column name)
             # in each column

    if os.path.exists(filename):
        print('This method will add table to existing file '+ filename)

    if len(data)!=len(dcols):
        print('The number of data columns does not equal the number of Data attributes')
        print('returning none')
        return None
    for i in range(len(headers)):
        tmp.append(header_char+' '+headers[i]+'\n')
    headers=tmp
    tmp=''

    for i in range(len(data)): #Line length stuff
        length=len(dcols[i])+1
        for j in range(len(data[i])):
            tmp2=data[i][j]
            if isinstance(data[i][j],float):
                tmp2='{:.3E}'.format(data[i][j])+' '
                data[i][j] = tmp2
            if len(str(tmp2))>length:
                length=len(str(tmp2))
        lengthList.append(length)

    tmp=''
    tmp1=''
    for i in range(len(dcols)):
        tmp1=dcols[i]
        if len(dcols[i]) < lengthList[i]:
            j=lengthList[i]-len(dcols[i])
            for k in range(j):
                tmp1+=' '
        tmp+=sep+tmp1
    tmp+='\n'
    dcols=tmp
    tmp=''
    for i in range(len(data[0])):
        for j in range(len(data)):
            if type(data[j][i]) == str:
                #match = re.match(r"([a-z]+)([0-9]+)",data[j][i], re.I)
                                    #items = match.groups()
                tmp1=data[j][i]#items[0].capitalize()+'-'+items[1]
                if len(str(data[j][i])) < lengthList[j]:
                        l=lengthList[j]-len(tmp1)
                        for k in range(l):
                                tmp1+=' '
                extra=''
            #else:
                            #        tmp1=data[j][i]
                            #        if len(data[j][i]) < lengthList[j]:
                            #                l=lengthList[j]-len(data[j][i]))
                            #                for k in xrange(l):
                            #                        tmp1+=' '


            tmp+=sep+tmp1
        lines.append(tmp+'\n')
        tmp=''

    f=open(filename,'a')
    for i in range(len(headers)):
        f.write(headers[i])
    f.write(dcols)
    for i in range(len(lines)):
        f.write(lines[i])

    f.close()
    return None




def writeGCE_table(filename,headers,data,dcols=['Isotopes','Yields','Z','A'],header_char='H',sldir='.',sep='&'):
    '''
    Method for writeing data in GCE format in Ascii files.
    Reads either elements or isotopes
    dcols[0] needs to contain either isotopes or elements

    Note the attribute name at position i in dcols will be associated
    with the column data at index i in data.
    Also the number of data columns(in data) must equal the number
    of data attributes (in dcols)
    Also all the lengths of that columns must all be the same.
    Input:
    filename: The file where this data will be written.
    Headers: A list of Header strings or if the file being written
         is of type trajectory, this is a List of strings
         that contain header attributes and their associated
         values which are seperated by a '='.
    dcols: A list of data attributes
    data:  A list of lists (or of numpy arrays).
    header_char  the character that indicates a header lines
    sldir: Where this fill will be written.
    sep: What seperatesa the data column attributes
    trajectory: Boolean of if we are writeing a trajectory type file
    '''

    import re
    from . import utils as u

    #check if input are elements or isotopes
    if not '-' in data[0][0]:
        iso_inp=False
        dcols=dcols+['Z']
    else:
        iso_inp=True
        dcols=dcols+['Z','A']
    #Attach Z and A
    if iso_inp==True:
        data.append([])
        data.append([])
        u.convert_specie_naming_from_h5_to_ppn(data[0])
        Z=u.znum_int
        A=u.amass_int
        for i in range(len(data[0])):
            zz=str(int(Z[i]))
            aa=str(int(A[i]))
            data[1][i]='{:.3E}'.format(data[1][i])+' '
            data[-2].append(zz)
            data[-1].append(aa)


    else:
        #in order to get Z , create fake isotope from element
        fake_iso=[]
        for k in range(len(data[0])):
            fake_iso.append(data[0][k]+'-99')
        #print fake_iso
            data.append([])
            u.convert_specie_naming_from_h5_to_ppn(fake_iso)
            Z=u.znum_int
        for i in range(len(data[0])):
            zz=str(int(Z[i]))
            data[1][i]='{:.3E}'.format(data[1][i])+' '
            data[-1].append(zz)


    if sldir.endswith(os.sep):
        filename = str(sldir)+str(filename)
    else:
        filename = str(sldir)+os.sep+str(filename)
    tmp=[] #temp variable
    lines=[]#list of the data lines
    lengthList=[]# list of the longest element (data or column name)
             # in each column

    if os.path.exists(filename):
        print('This method will add table to existing file '+ filename)

    if len(data)!=len(dcols):
        print('The number of data columns does not equal the number of Data attributes')
        print('returning none')
        return None
    for i in range(len(headers)):
        tmp.append(header_char+' '+headers[i]+'\n')
    headers=tmp
    tmp=''

    for i in range(len(data)): #Line length stuff
        length=len(dcols[i])+1
        for j in range(len(data[i])):
            tmp2=data[i][j]
            if isinstance(data[i][j],float):
                tmp2='{:.3E}'.format(data[i][j])+' '
                data[i][j] = tmp2
            if len(str(tmp2))>length:
                length=len(str(tmp2))
        lengthList.append(length)

    tmp=''
    tmp1=''
    for i in range(len(dcols)):
        tmp1=dcols[i]
        if len(dcols[i]) < lengthList[i]:
            j=lengthList[i]-len(dcols[i])
            for k in range(j):
                tmp1+=' '
        tmp+=sep+tmp1
    tmp+='\n'
    dcols=tmp
    tmp=''
    for i in range(len(data[0])):
        for j in range(len(data)):
            if type(data[j][i]) == str:
                #match = re.match(r"([a-z]+)([0-9]+)",data[j][i], re.I)
                                    #items = match.groups()
                tmp1=data[j][i]#items[0].capitalize()+'-'+items[1]
                if len(str(data[j][i])) < lengthList[j]:
                        l=lengthList[j]-len(tmp1)
                        for k in range(l):
                                tmp1+=' '
                extra=''
            #else:
                            #        tmp1=data[j][i]
                            #        if len(data[j][i]) < lengthList[j]:
                            #                l=lengthList[j]-len(data[j][i]))
                            #                for k in xrange(l):
                            #                        tmp1+=' '
            tmp+=sep+tmp1
        lines.append(tmp+'\n')
        tmp=''

    f=open(filename,'a')
    for i in range(len(headers)):
        f.write(headers[i])
    f.write(dcols)
    for i in range(len(lines)):
        f.write(lines[i])

    f.close()
    return None
