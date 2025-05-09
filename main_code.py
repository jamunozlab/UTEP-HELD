import time
import os, sys
import numpy as np 
start_time = time.time()
# --------- HELD
from utils.HELD import Force_Constants
from utils.default_dict import *
from utils.config import *
# --------- Phonopy
from utils.reorder_atoms_per_id import reorder_and_renumber_atoms
from utils.poscar_gen import write_poscar
from utils.dynmical_matrix_gen_phon import generate_force_constants_and_run_phonopy


# --------- HELD
for defect in defect_percentages:
    for temperature in range(n_temp_variations):
        for lattice_parameter in range(m_lattice_variations):
            lat_par = np.round(initial_lattice_parameter + 0.01 * lattice_parameter,2)
            temp = initial_temperature + 100 * temperature  
            box_length = lat_par * system_size
            print(f'Now we are going to HELD together, simulation for your material begins for defect {defect}% lattice parameter: {lat_par} Angs and temperature {temp}K')
            positions = np.zeros((natoms, 3, ntsteps))
            forces = np.zeros((natoms, 3, ntsteps))


            for t in range(ntsteps):
                data = np.genfromtxt(simulations_path + f'Simulation_Deffect_{defect}%/'+ f"NiTi_bcc_{lat_par:.2f}_{temp}/"  + dump_file_root + str(t*50), delimiter=' ', skip_header=9)

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
            print(f"dictn: {dictn}")

            # --------- Phonopy
            print('Dynamical matrix calculation begins, please wait')
            dir_path =  simulations_path + f'Simulation_Deffect_{defect}%/' +  f'NiTi_bcc_{lat_par:.2f}_{temp}/'
            input_file_path = dir_path + 'atoms_positions.data'
            output_file_path = dir_path + 'input.txt'
            n_atom_type1, n_atom_type2 = reorder_and_renumber_atoms(input_file_path, output_file_path)
            write_poscar(poscar_data, lat_par, path=dir_path, filename='POSCAR')
            command_copy_phonopy_files = f'cp {input_files_phonopy}'+'band.conf' '  ' f'{dir_path}' '&&' f' cp {input_files_phonopy}' + 'mesh.conf ' '  ' f'{dir_path}'
            os.system(command_copy_phonopy_files)   
            generate_force_constants_and_run_phonopy(dir_path,natoms,lat_par*system_size,lat_par,n_atom_type1,n_atom_type2,0,1,root="NiTi_bcc") 