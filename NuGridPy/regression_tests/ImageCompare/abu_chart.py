from builtins import str
from builtins import range
from ... import ppn as p
import os
import os.path

def load_chart_files(path = '.'):
    n = 39
    nsparsity = 10
    for cycle in range(0,n,nsparsity):
        cycle_str = str(cycle).zfill(2)
        os.system("wget -q --content-disposition --directory '" + path + "' "
                  +"'http://www.canfar.phys.uvic.ca/vospace/synctrans?TARGET="\
                  +"vos%3A%2F%2Fcadc.nrc.ca%21vospace%2Fnugrid%2Fdata%2Fprojects%2Fppn%2Fexamples%2F"\
                  +"ppn_Hburn_simple%2Fiso_massf000"+cycle_str+".DAT&DIRECTION=pullFromVoSpace&PROTOCOL"\
                  "=ivo%3A%2F%2Fivoa.net%2Fvospace%2Fcore%23httpget'")
        os.system("wget -q --content-disposition --directory '" + path + "' "
                  +"'http://www.canfar.phys.uvic.ca/vospace/synctrans?TARGET="\
                  +"vos%3A%2F%2Fcadc.nrc.ca%21vospace%2Fnugrid%2Fdata%2Fprojects%2Fppn%2Fexamples%2F"\
                  +"ppn_Hburn_simple%2FMasterAbuChart"+cycle_str+".png&DIRECTION=pullFromVoSpace&PROTOCOL"\
                  "=ivo%3A%2F%2Fivoa.net%2Fvospace%2Fcore%23httpget'")
    a=p.abu_vector(path)
    a.abu_chart(list(range(0,n,nsparsity)),plotaxis=[-1,16,-1,15], savefig=True, path=path)

if __name__ == "__main__":
    load_chart_files()
