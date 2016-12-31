#datadir='/tmp/orcinusB/NuGrid'
datadir='/apod2/NuGrid'
datadir='/nfs/rpod2/critter/CADC/NuGrid'
znum = 10
Mnum = 10
from NuGridPy import mesa as ms    
ms.__file__
from NuGridPy import nugridse as mp
import random
mp.__file__

ms.set_nugrid_path(datadir)
mp.set_nugrid_path(datadir)
Mminmax=(0.9,26)
Ms = (Mminmax[1]-Mminmax[0])*random.random(Mnum) + Mminmax[0]
Zminmax=(0.00005,0.021)
Zs = 10**((log10(Zminmax[1])-log10(Zminmax[0]))*random.random(znum) + log10(Zminmax[0]))
for Z in Zs:
    for M in Ms:
        # mesa history
        figure(1)
        mesa = ms.star_log(mass=M,Z=Z)
        xthing=random.choice(mesa.cols.keys())
        ything=random.choice(mesa.cols.keys())
        mesa.plot(xthing,ything)
        title("MESA history")
        savefig("mesa_histM"+'%5.3f'%M+'Z'+'%5.3e'%Z+".png")
        close(1)
        # mesa profile
        figure(2)
        random_cycle = int(random.random()*mesa.get('model_number')[-1])
        profile = ms.mesa_profile(mass=M,Z=Z,num=random_cycle)
        xthing=random.choice(profile.cols.keys())
        ything=random.choice(profile.cols.keys())
        profile.plot(xthing,ything)
        title("MESA profile: "+str(random_cycle))
        savefig("mesa_prof"+'%5.3f'%M+'Z'+'%5.3e'%Z+'Cyc'+str(random_cycle)+".png")
        close(2)
        # post-processing files
        figure(3)
        pp = mp.se(mass=M,Z=Z)
        iso_index = int(random.random()*len(pp.se.isotopes))
        random_cycle = int(round(random_cycle/20.)*20)
        pp.plot('mass',pp.se.isotopes[iso_index],fname=random_cycle)
        title("PPD profile: "+str(random_cycle))
        savefig("ppd_prof"+'%5.3f'%M+'Z'+'%5.3e'%Z+'Cyc'+str(random_cycle)+".png")
        close(3)


        
