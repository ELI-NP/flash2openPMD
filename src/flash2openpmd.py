#!/usr/bin/env python
"""
This file is part of Flash2OpenPMD software, which
converts the FLASH files to openPMD
Copyright (C) 2022, Flash2OpenPMD contributors
Author: Borsos Andrei Paul, Jian Fuh Ong

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import yt
import os
import numpy as np
import sys
from scipy import constants as sc
import openpmd_api  as io
import glob
from scipy import interpolate

class Convert(object):
    def __init__(self,run_directory,filename):
        """
        Parameters
        ----------
        run_directory : string
            The path to the directory where the output of FLASH are.

        filename : string
            The output filename of FLASH. (e.g. lasslab_hdf5_plt_cnt_0043)
        """

        if run_directory is None:
            raise ValueError('The run_directory parameter can not be None!')

        self.path_to_sim = os.path.join(run_directory,filename)

    def read4Flash(self,level=4):
        """
        Read the FLASH output

        Parameters
        ----------
        level : int
            The level of refinement using yt-project. Currently fixed to level=4

        Returns
        -------
        A 2d array containing the required density, shape and maximum boundary in each direction.

        """

        ds = yt.load(self.path_to_sim)

        if level==0:
            scale = 1
        else:
            scale = np.array([2**level,2**level,1])

        ds.force_periodicity()
        all_data = ds.covering_grid(level=level,
                                    left_edge=ds.domain_left_edge,
                                    dims=ds.domain_dimensions * scale)

        density = all_data["gas", "density"] # density in [g/cm3]

        x_min, y_min, z_min = ds.domain_left_edge
        x_max, y_max, z_max = ds.domain_right_edge

        nx, ny, nz = density.shape

        return density, nx, ny, nz, x_max.v, y_max.v, z_max.v, x_min.v, y_min.v, z_min.v

    def get_data(self):

        dens, nx, ny, nz, x_max, y_max, z_max, x_min, y_min, zmin = self.read4Flash()

        """
        Selecting the data that will be wrote to be used in openopmd

        Parameters
        ----------
        density : array 
            Includes the values of density in different coordinates in the XoY plane as a 2D array
            
        nx, ny, nz : int
            The length of the x,y 1D arrays and the numbers of lines and columns of the density matrix/z array
            
        
        x_max, y_max, z_max :float
            The maximum boundary in x, y, and z directions in the unit of centimeter (unit of FLASH code)
            
        refinex, refiney : int
            The level of refinement we use for the interpolation in order to obtain a more clear simulation
            (How many times more points we will have in the final figure on the x and y axis)
            
        xnew, ynew, znew :array
            The new coordonates on the x and y axis and the znew value at each intersecton of x column and y line after interpolation
            
        Returns
        -------
        A 2d array containing the required density after the interpolation

        """
        refinex = input("Level of x-axis refinement (int):\n")
        refiney = input("level of y-axis refinement (int):\n")
        refinex = int(refinex)
        refiney = int(refiney)
        density = np.zeros([ny-1,nx-1])

        for i in range(0,nx-1,1):
            for j in range(0,ny-1,1):
                density[j,i] = dens[i,j,0]

        x = np.linspace(x_min, x_max, nx-1)
        y = np.linspace(y_min, y_max, ny-1)

        z = density

        f = interpolate.interp2d(x, y, z, kind='linear')

        xnew = np.linspace(x_min, x_max, (nx-1)*refinex)
        ynew = np.linspace(y_min, y_max, (ny-1)*refiney)

        znew = f(xnew, ynew)

        return np.array(znew, dtype='float32')

    def write2openpmd(self,n_e_input):

        """
        Writing the data to be rady to use in openopmd

        Parameters
        ----------
        n_e_input : array
            The input data we want to write for openpmd(in our case density)

        series_out :record
            The file that contains the data used for openpmd

        k :int
            The first iteration of the system/the starting point

        n_e_out : array
            The output array used for openpmd that will contain the interpolated density values

        n_e_mrc :array
            The array used to record the data in the output array

        dataset :array
            An array with the same dimensions and values a as n_e_input

        Returns
        -------
        Creates an hdf5 file that contains the input needed for openpmd

        """
        n_e_input = n_e_input.T
        
        n_e_input = n_e_input/np.max(n_e_input) # Normalized to 1

        nume_fisier = input("Output file name :\n")
        series_out = io.Series("output/%s.h5" %(nume_fisier),io.Access.create)

        k = series_out.iterations[0]

        series_out.author = "Your Name <Your@email>"
        # record - again,important to specify as scalar
        n_e_out = k.meshes["e_density"]
        n_e_mrc = n_e_out[io.Mesh_Record_Component.SCALAR]

        dataset = io.Dataset(
            n_e_input.dtype,
            n_e_input.shape)
        n_e_mrc.reset_dataset(dataset)

        n_e_out.set_attribute('dataOrder','C')
        n_e_out.set_attribute('axisLabels',['z','y','x'])
        n_e_mrc.store_chunk(n_e_input)
        # After registering a data chunk such as x_data and y_data,
        # it MUST NOT be modified or deleted until the flush() step is performed!
        series_out.flush()

        del series_out
