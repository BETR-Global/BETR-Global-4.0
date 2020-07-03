################################################################################
'''Module "emissions.py", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module reads emission files and hands out emission-vectors'''
################################################################################
from numpy import *
import scipy.sparse as sp
import inspect
from globalz import *
from helpers import *
import sys
import pdb
    
class Emission:
    """ reads emission file and saves as static dictionary
        The compartments start with 1 """  
       
    def __init__(self,fn):
        self.em = {}
        try:
            f=open(fn, 'r')
        except IOError:
            print("emissions.py: emission file %s not found. Aborting!") % (fn)
            sys.exit(1)
        lines=f.readlines()
        f.close
        for lin in lines:
            if (lin[0] == '#' or lin == ''):
                continue
            #Fangyuan changed [m,r,c,val] to [m,r,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,val1,val2,val3,val4,val5,val6,val7,val8,val9,val10]
            #Fangyuan changed append((int(float(r)),int(float(c)),float(val))) to
            #append((int(float(r)),int(float(c1)),int(float(c2)),int(float(c3)),int(float(c4)),
            #int(float(c5)),int(float(c6)),int(float(c7)),int(float(c8)),int(float(c9)),
            #int(float(c10)),float(val1),float(val2),float(val3),float(val4),float(val5),
            #float(val6),float(val7),float(val8),float(val9),float(val10)))
            [m,r,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,val1,val2,val3,val4,val5,val6,val7,val8,val9,val10]=[x.lstrip() for x in \
                                             [x.rstrip() for x in lin.split()]]
            self.em.setdefault(int(float(m)),[]).\
                             append((int(float(r)),int(float(c1)),int(float(c2)),int(float(c3)),int(float(c4)),
                                     int(float(c5)),int(float(c6)),int(float(c7)),int(float(c8)),int(float(c9)),
                                     int(float(c10)),float(val1),float(val2),float(val3),float(val4),float(val5),
                                     float(val6),float(val7),float(val8),float(val9),float(val10)))
            
    def get_emission(self,timeidx,m,type='csc'):
        """ hands out emission vector for timeindex (starts with 0).
        default type csc returns vector in csc-format
        type=array returns 1D-array."""
        emvec=zeros(m.matdim)
        try:
            ems=self.em[timeidx+1]
            for source in ems:
                #Fangyuan separate emissions to gas phase, fine and coarse particles. 
                cell1 = tocell(source[0],source[1],m); emvec[cell1]=source[11]                
                cell2 = tocell(source[0],source[2],m); emvec[cell2]=source[12]              
                cell3 = tocell(source[0],source[3],m); emvec[cell3]=source[13]
                cell4 = tocell(source[0],source[4],m); emvec[cell4]=source[14]                
                cell5 = tocell(source[0],source[5],m); emvec[cell5]=source[15]              
                cell6 = tocell(source[0],source[6],m); emvec[cell6]=source[16]
                cell7 = tocell(source[0],source[7],m); emvec[cell7]=source[17]                
                cell8 = tocell(source[0],source[8],m); emvec[cell8]=source[18]              
                cell9 = tocell(source[0],source[9],m); emvec[cell9]=source[19]
                cell10 = tocell(source[0],source[10],m); emvec[cell10]=source[20]
        except KeyError:
            pass
        ## compress
        emvec=compressvec(emvec,m.killidx)
        if type=='csc':
             emvec=sp.csr_matrix(emvec).T
        return emvec
    
    # def get_csc(self,timeidx,m):
    #     """ hands out emission vector for timeindex (starts with 0) as csc-matrix"""
    #     emvec=zeros(m.matdim)
    #     try:
    #         ems=self.em[timeidx+1]
    #         for source in ems:
    #             cell = tocell(source[0],source[1],m)
    #             emvec[cell]=source[2]
    #     except KeyError:
    #         pass
    #     ## compress and change to csc-vector
    #     emvec=sp.csr_matrix(compressvec(emvec,m.killidx)).T
    #     # emvec=compressvec(emvec,m.killidx)
    #     ## change units of mol/h to fugaity [Pa/h]
    #     # season=mod(timeidx,m.nots)
    #     # emvec=(m.ZVinvlist[season]*sp.csr_matrix(emvec).T).asformat('csc')
    #     return emvec

