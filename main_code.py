import time
import os, sys
import numpy as np 
start_time = time.time()


from utils.HELD import Force_Constants
from utils.default_dict import *
from utils.config import *

for defect in defect_percentages:
    for temperature in range(n_temp_variations):
        for lattice_parameter in range(m_lattice_variations):
            lat_par = np.round(initial_lattice_parameter + 0.01 * lattice_parameter,2)
            temp = initial_temperature + 100 * temperature  
            box_length = lat_par * system_size
            print(lat_par,temp,defect)
            positions = np.zeros((natoms, 3, ntsteps))
            forces = np.zeros((natoms, 3, ntsteps))


            for t in range(ntsteps):
                data = np.genfromtxt(simulations_path + f'Simulation_Deffect_{defect}%/'+ f"NiTi_bcc_{lat_par:.2f}_{temp}/"  + dump_file_root + str(t*50), delimiter=' ', skip_header=9)#, skip_header=, skip_header=9
                # print( np.max(data[:, 0]))

                for row in data:
                    positions[int(row[0])-1, :, t] = row[2:5]
                    forces[int(row[0])-1, :, t] = row[5:8]           
            id_lattice = np.zeros((natoms, 3)) # id_lattice is the ideal underlying lattice
            data = np.genfromtxt(simulations_path + f'Simulation_Deffect_{defect}%/'+ f"NiTi_bcc_{lat_par:.2f}_{temp}/"  +pos_file, delimiter=' ', skip_header=9)#, delimiter=' ', skip_header=11)

            for row in data:
                id_lattice[int(row[0])-1, :] = row[2:5]

            id_lattice = positions[:, :, 0]


            out_path =simulations_path + f'Simulation_Deffect_{defect}%/'+ f"NiTi_bcc_{lat_par:.2f}_{temp}/"
            bvk_fcc = Force_Constants(positions=positions, forces=forces,alat=box_length, ideal_lattice=id_lattice,out_path =out_path,lat_par = lat_par )
            ideal_distances, ideal_dist_sca =bvk_fcc.ideal_lat_dist()
            dictn = bvk_fcc.get_fc_steps(normal_fit=True)
            
            with open(out_path+fc_filename, 'w') as f:
                for c in range(bvk_fcc.fc):
                    print(bvk_fcc.index[c] + ', ', end='', file=f)

                print('\n', file=f)

                for c in range(bvk_fcc.fc):
                    print(str(dictn[bvk_fcc.index[c]]) + ', ', end='', file=f)

                print('\n', file=f)

                for t in range(ntsteps):
                    for c in range(bvk_fcc.fc):
                        print(str(bvk_fcc.phiv_steps[c, t]) + ', ', end='', file=f)

                    print('', file=f)
            end_time = time.time()
            mnt = int((end_time - start_time) / 60)
            if mnt >= 60:
                hr = int(mnt / 60)
                mnt = int(mnt % 60)

            else:
                hr = 0

            scn = int((end_time - start_time) % 60.)
            # print(str(hr).zfill(2) + ':' + str(mnt).zfill(2) + ':' + str(scn).zfill(2))
            print(f"dictn: {dictn}")
            