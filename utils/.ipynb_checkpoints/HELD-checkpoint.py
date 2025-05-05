#!/usr/bin/python
# -*- coding: utf-8 -*-

from utils.config import simulations_path
import numpy as np
from math import sqrt, isclose
from scipy.stats import norm
import matplotlib.pyplot as plt

#class part BCC-B2
class Force_Constants:
    #Values of Born-von Kármán model for up to fifth nearest neighbors of a fcc structure
    #for a hdf5 file. args = np array[n, 3, t], np array[n, 3, t], np array[n, 3],  float, np array[n, 3], int[1, 5]
    def __init__(self, positions, forces, alat, lat_par, ideal_lattice,out_path, under_lattice=np.zeros(0), NN=5):
        if NN < 1 or NN > 5 or not isinstance(NN, int):
            raise TypeError('Nearest neighbors argument must be integer in range [1, 5].')

            #/Users/dajuarez4/Downloads/Gamma-Iron_HELD/LAMMPS_MD/Test
        self.out_path = out_path
        self.positions = positions
        self.forces = forces
        self.alat = alat
        self.ide_lat = ideal_lattice
        self.nn = NN
        self.n = positions.shape[0] #Number of atoms in supercell
        self.n_steps = positions.shape[2] #Number of time steps
        self.steps = range(self.n_steps)
        # self.a_val = self.alat / 4
        self.a_val = lat_par  # self.alat / ((self.n)**(1./3.)) #Side length of each cube in fcc
        if under_lattice.size == 0:
            self.und_lat = self.ide_lat

        else:
            self.und_lat = under_lattice


        self.fc = [4, 8, 14, 18, 22][self.nn-1] #Amount of constants according to NN
        self.index = [r'$A^{xx}_0$', r'$B^{xx}_0$', r'$c^{xx}_1$', r'$c_1^{xy}$', r'$A_2^{xx}$', r'$A_2^{yy}$', r'$B_2^{xx}$', r'$B_2^{yy}$',
                      r'$A_3^{xx}$', r'$A_3^{yy}$', r'$A_3^{yz}$', r'$B_3^{xx}$', r'$B_3^{yy}$',
                      r'$B_3^{yz}$', r'$c_4^{xx}$', r'$c_4^{yy}$', r'$c_4^{xy}$', r'$c_4^{yz}$', r'$A_5^{xx}$',
                      r'$A_5^{xy}$', r'$B_5^{xx}$', r'$B_5^{xy}$'] #Index of force constants
        self.sm = np.array([[1., 0., 8., 0., 2., 4., 0., 0., 4., 8., 0., 0., 0., 0., 8., 16., 0., 0., 8., 0., 0., 0.], 
                        [0., 1., 8., 0., 0., 0., 2., 4., 0., 0., 0., 4., 8., 0., 8., 16., 0., 0., 0., 0., 8., 0.]])

        self.disp_mat = np.zeros((3 * self.n, self.fc, self.n_steps)) #Matrices of displacements from equilibrium by time step
        self.phiv_steps = np.zeros((self.fc, self.n_steps)) #Vector of all the unique values of phi



    def ideal_lat_dist(self):
        ideal_distances = np.zeros((self.n, 3, self.n))
        ideal_dist_sca = np.zeros((self.n, self.n))

        for i in range(self.n):
            for j in range(self.n):
                if i != j:
                    ideal_distances[i, :, j] = self.ide_lat[j, :] - self.ide_lat[i, :] #Ideal distance vector from i to j

                    #Apply mic:
                    for d in range(3):
                        if ideal_distances[i, d, j] > (self.alat / 2):
                            ideal_distances[i, d, j] -= self.alat

                        elif ideal_distances[i, d, j] <= (-self.alat / 2):
                            ideal_distances[i, d, j] += self.alat

                    ideal_dist_sca[i, j] = np.linalg.norm(ideal_distances[i, :, j]) #Distance scalar

        return ideal_distances, ideal_dist_sca

    #Identify up to fifth nearest neighbors for every atom
    def nearest_neighbors(self):
        if not hasattr(self, 'ideal_dist_sca'):
            self.ideal_distances, self.ideal_dist_sca = self.ideal_lat_dist()

        neighbors = [] #List of dictionaries holding a list of the IDs for all types of nn

        #Distances to each type of neighbor:
        first = (sqrt(3) * self.a_val) / 2
        second = self.a_val
        third = sqrt(2) * self.a_val
        fourth = (sqrt(11) * self.a_val) / 2
        fifth = sqrt(3) * self.a_val
        self.nn_dist = [first, second, third, fourth, fifth]

        for i in range(self.n):
            i_neigh = [[i]] 
            for prox in range(self.nn):
                i_neigh.append([]) #Empty placeholder lists for neighbors of i

            #Add each list of nn to the placeholder dictionary, append the dictionary to the list:
            for j in range(self.n):
                if i != j:
                    for prox in range(1,self.nn+1):
                        if isclose(self.ideal_dist_sca[i, j], self.nn_dist[prox-1], rel_tol=.0001):
                            i_neigh[prox].append(j)

            neighbors.append(i_neigh)

        return neighbors

    #Arrange part of U corresponding to time t in F=U*phi
    def BvK_matrix(self, t):
        if not hasattr(self, 'neighbors'):
            self.neighbors = self.nearest_neighbors()

        disp_mat = np.zeros((3 * self.n, self.fc)) #Empty matrix U of displacements from equilibrium
        for i in range(self.n):
            for prox in range(self.nn+1):
                for j in self.neighbors[i][prox]:
                    displacements = self.positions[j, :, t] - self.und_lat[j, :] #Distance vector from i to j at t

                    #Apply mic:
                    for d in range(3):
                        if displacements[d] > (self.alat / 2):
                            displacements[d] -= self.alat

                        elif displacements[d] <= (-self.alat / 2):
                            displacements[d] += self.alat

                    #Specify the type of pair according to proximity and type of atoms:
                    if prox == 0:
                        if i < self.n/2:
                            tpe = 0

                        if i >= self.n/2:
                            tpe = 1

                    if prox == 1:
                        tpe = 2

                    if prox == 2:
                        if i < self.n/2:
                            tpe = 3

                        if i >= self.n/2:
                            tpe = 4

                    if prox == 3:
                        if i < self.n/2:
                            tpe = 5

                        if i >= self.n/2:
                            tpe = 6

                    if prox == 4:
                        tpe = 7

                    if prox == 5:
                        if i < self.n/2:
                            tpe = 8

                        if i >= self.n/2:
                            tpe = 9

                    #Arrange matrix U according to similar terms:
                    m1 = np.array([[displacements[0], 0., 0.],
                                   [0., displacements[1], 0.],
                                   [0., 0., displacements[2]]])
                    m2 = np.array([[0., displacements[1], displacements[2]],
                                   [displacements[0], 0., displacements[2]],
                                   [displacements[0], displacements[1], 0.]])
                    for a in range(3):
                        if self.ideal_distances[i, a, j] < 0:
                            m2[a, :] *= -1
                            m2[:, a] *= -1

                    if prox == 0:
                        m1r = np.diag(m1)
                        if tpe == 0:
                            disp_mat[3 * i: 3 * (i+1), 0] += m1r

                        if tpe == 1:
                            disp_mat[3 * i: 3 * (i+1), 1] += m1r

                    if prox == 1 or prox == 5:
                        m1r = np.diag(m1)
                        m2r = np.sum(m2, axis=1)
                        if tpe == 2:
                            disp_mat[3 * i: 3 * (i+1), 2] += m1r
                            disp_mat[3 * i: 3 * (i+1), 3] += m2r

                        if tpe == 8:
                            disp_mat[3 * i: 3 * (i+1), 18] += m1r
                            disp_mat[3 * i: 3 * (i+1), 19] += m2r

                        if tpe == 9:
                            disp_mat[3 * i: 3 * (i+1), 20] += m1r
                            disp_mat[3 * i: 3 * (i+1), 21] += m2r

                    if prox == 2:
                        m1r = np.zeros((3, 2))
                        for a in range(3):
                            if self.ideal_distances[i, a, j] == 0.:
                                m1r[:, 1] += m1[:, a]

                            else:
                                m1r[:, 0] = m1[:, a]

                        if tpe == 3:
                            disp_mat[3 * i: 3 * (i+1), 4:6] += m1r

                        if tpe == 4:
                            disp_mat[3 * i: 3 * (i+1), 6:8] += m1r

                    if prox == 3:
                        m1r = np.zeros((3, 2))
                        for a in range(3):
                            if self.ideal_distances[i, a, j] == 0.:
                                m1r[:, 0] = m1[:, a]

                            else:
                                m1r[:, 1] += m1[:, a]

                            if self.ideal_distances[i, a, j] == 0.:
                                m2[a, :] = np.zeros(3)
                                m2[:, a] = np.zeros(3)

                        m2r = np.sum(m2, axis = 1)

                        if tpe == 5:
                            disp_mat[3 * i: 3 * (i+1), 8:10] += m1r
                            disp_mat[3 * i: 3 * (i+1), 10] += m2r

                        if tpe == 6:
                            disp_mat[3 * i: 3 * (i+1), 11:13] += m1r
                            disp_mat[3 * i: 3 * (i+1), 13] += m2r

                    if prox == 4:
                        m1r = np.zeros((3, 2))
                        m2r = np.zeros((3, 2))
                        for a in range(3):
                            if abs(self.ideal_distances[i, a, j]) == 1.5 * self.a_val:
                                m1r[:, 0] = m1[:, a]
                                m2r[a, 0] = np.sum(m2[a, :])
                                m2r[(a+1) % 3, 0] = m2[(a+1) % 3, a]
                                m2r[(a+2) % 3, 0] = m2[(a+2) % 3, a]
                                m2r[(a+1) % 3, 1] = m2[(a+1) % 3, (a+2) % 3]
                                m2r[(a+2) % 3, 1] = m2[(a+2) % 3, (a+1) % 3]

                            else:
                                m1r[:, 1] += m1[:, a]

                        disp_mat[3 * i: 3 * (i+1), 14:16] += m1r
                        disp_mat[3 * i: 3 * (i+1), 16:18] += m2r

        return disp_mat #Matrix of sums of displacements from equilibrium

    #Calculate force constant values for list of steps, normal_fit defines wether the mean of the distribution is given 
    #instead of all time steps
    def get_fc_steps(self, steps_list=[], normal_fit=False):
        if not steps_list:
            steps_list = self.steps

        Fv = np.zeros((3 * self.n) + 2) #Vector of all the force vetors at every time step appended
        Disp = np.zeros(((3 * self.n) + 2, self.fc)) #Empty matrix of displacements from equilibrium
        for t in steps_list:
            if self.phiv_steps[:, t].all() == 0.:
                if self.disp_mat[:, :, t].all() == 0.:
                    self.disp_mat[:, :, t] = self.BvK_matrix(t)

                Disp[:3 * self.n, :] = self.disp_mat[:, :, t] #Fill equations related to force

                #Append last row so that the sum of all phi_{alfa, alfa} equals 0:
                Disp[3 * self.n:, :] = self.sm[:, :self.fc]

                #Arrange all force vectors at t in a single vector:
                for i in range(self.n):
                    Fv[3 * i: 3 * (i + 1)] = -self.forces[i, :, t]

                Disp_inv = np.linalg.pinv(Disp)
                self.phiv_steps[:, t] = np.matmul(Disp_inv, Fv)

                self.phiv_steps[:2, t] = -np.matmul(self.sm[:, 2:], self.phiv_steps[2:, t])

            print(t)

        if not normal_fit:
            dicts_list = {}
            for t in steps_list:
                phiv_dict = {}
                for tpe in range(self.fc):
                    phiv_dict[self.index[tpe]] = self.phiv_steps[tpe, t]

                dicts_list[str(t)] = phiv_dict

            return dicts_list #Dictionary of dictionaries with unique force constant values at each time step

        else:
            phiv_fit = np.zeros((self.fc, 2)) #Fit to be made for the vector of phi
            for tpe in range(self.fc):
                phiv_fit[tpe, 0], phiv_fit[tpe, 1] = norm.fit(np.take(self.phiv_steps[tpe, :], steps_list)) #Fit the Gaussian
                plt.hist(np.take(self.phiv_steps[tpe, :], steps_list), bins=40, density=True)
                xmin, xmax = plt.xlim()
                x = np.linspace(xmin, xmax, 40)
                pdf = norm.pdf(x, phiv_fit[tpe, 0], phiv_fit[tpe, 1])
                plt.plot(x, pdf, 'r',color='red')
                plt.title('Distribution of ' + str(self.index[tpe]),fontsize=34)
                plt.xlabel(r'Force constant value $(eV/ \AA^2)$',fontsize=24)
                plt.ylabel('Density',fontsize=30)
                plt.tick_params(axis='both', which='major', labelsize=30)
                plt.savefig(self.out_path+'Normal fit A2' + self.index[tpe] + '.png',bbox_inches='tight', pad_inches=0.1, transparent=True)
                plt.show()
                plt.close()

            fit_dict = {}
            for tpe in range(self.fc):
                fit_dict[self.index[tpe]] = phiv_fit[tpe, 0]

            return fit_dict #Dictionary with unique force constant values
