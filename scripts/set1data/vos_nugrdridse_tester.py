# mount nugrid CADC VOspace: 
# $ mountvofs --readonly --cache_nodes --cache_dir=/tmp/vosCache --mountpoint=/tmp/nugrid --vospace=vos:nugrid
# authenticate vos with getCert
# !!! I think one should be able to read-mount vos:nugrid without authentication !!!
# execute on commandline: 
# $ python   nugridse_demo.py 

nugrid_path = '/tmp/nugrid_fh'
#nugrid_path = '/apod2/fherwig/cadc'
nugrid_path='/nfs/rpod2/critter/CADC/NuGrid'
use_seeker = False

from NuGridPy import nugridse as ns
ns.set_nugrid_path(nugrid_path)
if use_seeker:
    p2=ns.se(mass=2, Z=0.001)
else:
    print("Opening directory: "+nugrid_path+'/data/set1/set1.3a/ppd_wind/M2.000Z0.0060/H5_out')
    p2 = ns.se(nugrid_path+'/data/set1/set1.3a/ppd_wind/M2.000Z0.0060/H5_out')
p2.plot('mass','Ba-138',fname=24600)
