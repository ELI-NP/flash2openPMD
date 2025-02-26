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
            
        self.run_directory = run_directory
        self.path_to_sim = os.path.join(run_directory,filename)       

    def read4Flash(self,x_min,y_min,z_min,fields,level):
        """
        Read the FLASH output

        Parameters
        ----------
        x_min : float
            The minimum boundary in x direction in the unit of centimeter (unit of FLASH code)
        
        y_min : float
            The minimum boundary in y direction in the unit of centimeter (unit of FLASH code)

        z_min : float
            The minimum boundary in y direction in the unit of centimeter (unit of FLASH code)
            z_min=0 for 2D simulation
        
        field : string
        	The field data (unit of FLASH code), i.e. fields="El_number_density"
        
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
            if z_min == 0:
                scale = np.array([2**level,2**level,1])
            else:
                scale = np.array([2**level,2**level,2**level])

        ds.force_periodicity()
        all_data = ds.covering_grid(level=level, 
                                    left_edge=[x_min, y_min, z_min], 
                                    dims=ds.domain_dimensions * scale)

        density = all_data[('gas',f"{fields}")]

        x_max, y_max, z_max = ds.domain_right_edge

        nx, ny, nz = density.shape

        return density, nx, ny, nz, x_max.v, y_max.v, z_max.v

    def get_data(self,x_min,y_min,z_min,fields,level):

        dens, nx, ny, nz, x_max, y_max, z_max = self.read4Flash(x_min,y_min,z_min,fields,level)

        if z_min == 0.0:
            density = np.array(dens[:,:,0], order='C', dtype='float32')
        else:
            density = np.array(dens.T, order='C', dtype='float32')
		
        return density

    def write2openpmd(self,density_input,species,output_name):

        series_out = io.Series(
        self.run_directory + f"{species}_{output_name}_0.h5",
        io.Access.create)

        k = series_out.iterations[0]

        series_out.author = "Your Name <Your@email>"
        # record - again,important to specify as scalar
        n_e_out = k.meshes[f"{species}_density"]
        n_e_mrc = n_e_out[io.Mesh_Record_Component.SCALAR]

        dataset = io.Dataset(
            density_input.dtype,
            density_input.shape)
        n_e_mrc.reset_dataset(dataset)

        n_e_out.set_attribute('dataOrder','C')
        n_e_out.set_attribute('axisLabels',['z','y','x'])
        n_e_mrc.store_chunk(density_input)
        # After registering a data chunk such as x_data and y_data,
        # it MUST NOT be modified or deleted until the flush() step is performed!
        series_out.flush()

        del series_out
