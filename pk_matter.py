from nbodykit.io.gadget import Gadget1File
import numpy as np
from nbodykit.lab import UniformCatalog, FFTPower
from nbodykit.source.catalog.array import ArrayCatalog
import sys

if len(sys.argv) != 8:
    print("Program needs 8 parameters!")
    print("exit...")
    sys.exit()

nome_script, name_sim, cosmo, directory, snap, Nfiles, boxsize, Nmesh = sys.argv

N_files = int(Nfiles)
print(N_files)
box_size = float(boxsize)
print(box_size)
N_mesh = int(Nmesh)
print(Nmesh)

num=0

for i in range(N_files):
    dirin=directory+'/snapdir_0{}/snap_0{}.{}'.format(snap, snap, i)
    file = Gadget1File(dirin)  ##, columndefs=[('Position', ('auto', 3), 'all'), ('GadgetVelocity', ('auto', 3), 'all'), ('ID', 'auto', 'all'), ('Mass', 'f8', (0))], ptype=1)
    N = file.size
    num=num+N
    print('Read snap ', i)
    print(dirin)
    if (i==0):
        pos = file.read(file[['Position']],0, N)
        pos['Position']=pos['Position']/1000.0
    else:
        pos1=file.read(file[['Position']],0, N)
        pos1['Position']=pos1['Position']/1000.0
        pos=np.append(pos, pos1, axis=0)


cat = ArrayCatalog(pos)
mesh = cat.to_mesh(Nmesh=N_mesh, BoxSize=box_size, window=None, resampler='tsc', interlaced=True, compensated=True)
print(mesh)

mesh.save('mesh-real-matter-cdm-{}-{}-{}.bigfile'.format(name_sim, cosmo, snap), mode='real', dataset='Field')

r = FFTPower(mesh, mode='1d')
r.save('fftpower-matter-cdm-{}-{}-{}.json'.format(name_sim, cosmo, snap))

r2 = FFTPower.load('fftpower-matter-cdm-{}-{}-{}.json'.format(name_sim, cosmo, snap))

k= r2.power['k']
pk=r2.power['power'].real
shot_noise = r2.attrs['volume'] / r2.attrs['N1']
sn_array = np.repeat(shot_noise, len(pk))
np.savetxt('pk-matter-cdm-{}-{}-{}.txt'.format(name_sim, cosmo, snap), np.column_stack([k, pk, sn_array]), header='# 0:k[h/Mpc], 1:pk[(Mpc/h)^3], 2:shot_noise')

