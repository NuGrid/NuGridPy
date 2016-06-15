from NuGridPy import ppn as p
import os
for cycle in range(0,39):
    cycle_str = str(cycle).zfill(2)
    os.system("wget --content-disposition  'http://www.canfar.phys.uvic.ca/vospace/synctrans?TARGET="\
              +"vos%3A%2F%2Fcadc.nrc.ca%21vospace%2Fnugrid%2Fdata%2Fprojects%2Fppn%2Fexamples%2F"\
              +"ppn_Hburn_simple%2Fiso_massf000"+cycle_str+".DAT&DIRECTION=pullFromVoSpace&PROTOCOL"\
              "=ivo%3A%2F%2Fivoa.net%2Fvospace%2Fcore%23httpget'")
a=p.abu_vector('.')
a.abu_chart(range(0,39),plotaxis=[-1,16,-1,15],savefig=True)
