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

    def read4Flash(self,fields,level):
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

        ND = ds.dimensionality

        print(f"Data geometry: {ds.geometry}")
        print(f"Data dimensionality: {ND}D")

        if ND == 2:
            scale = np.array([2**level,2**level,1])
        elif ND == 3:
            scale = np.array([2**level,2**level,2**level])

        ds.force_periodicity()
        all_data = ds.smoothed_covering_grid(level=level, 
                                    left_edge=ds.domain_left_edge, 
                                    dims=ds.domain_dimensions * scale,
                                    )
        
        if ds.geometry == 'cylindrical':
            density = self.cylindricalRotateAxial(all_data,fields=f"{fields}")
        else:
            density = all_data[('gas',f"{fields}")]

        return density, ND, ds.geometry
    
    def cylindricalRotateAxial(self,array,fields):
        """
        Rotate the 2D cylindrical data to 3D Cartesian data

        Parameters
        ----------
        array : 2D array
            The 2D cylindrical data

        Returns
        -------
        A 3D array containing the rotated Cartesian data.

        """

        rho = array[(f'{fields}')]
        r = array[("index","r")].v[:,0].reshape(-1)
        z = array[("index","z")].v[0,:].reshape(-1)
        R, Z = np.meshgrid(r, z)

        Nx = len(2*r)
        Ny = len(2*r)
        X = np.linspace(-r.max(), r.max(), Nx)
        Y = np.linspace(-r.max(), r.max(), Ny)
        Z_out = z.copy()

        xx, yy = np.meshgrid(X, Y, indexing='xy')   # shape (Ny, Nx)
        r_grid_xy = np.sqrt(xx**2 + yy**2)          # radial distance at each (x,y)

        Nz = len(Z_out)
        vol = np.zeros((Nz, Ny, Nx), dtype=rho.dtype)  # (z, y, x)

        if not np.all(np.diff(r) > 0):
            raise ValueError("r must be strictly increasing")
        
        print("Converting cylindrical to 3D Cartesian coordinates...")

        rflat = r_grid_xy.ravel()
        for iz in range(Nz):
            slice_r = rho[:, iz, 0]             
            vals_flat = np.interp(rflat, r, slice_r, left=slice_r[0], right=slice_r[-1])
            vol[iz, :, :] = vals_flat.reshape(Ny, Nx)

        return vol

    def get_data(self,fields,level,dtype=None):

        dens, ND, geometry = self.read4Flash(fields,level)

        if ND == 2:
            if geometry == 'cylindrical':
                density = np.array(dens[2:-2,2:-2,2:-2].T, order='C', dtype=dtype)
            else:
                density = np.array(dens[2:-2,2:-2,0], order='C', dtype=dtype)
        elif ND == 3:
            density = np.array(dens[2:-2,2:-2,2:-2].T, order='C', dtype=dtype)
		
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
