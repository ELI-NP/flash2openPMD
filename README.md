Conversion from FLASH to openPMD
================================

FLASH => openPMD (ELI-NP, Romania)

## QuickStart
To begin, the following dependencies are required:

#### Dependencies:

- numpy, scipy, yt
- openPMD-api

#### Running the conversion:

See `Convert.ipynb` for the flow of conversion

The converted file has the extension `.h5` with custom name

Simulation output from FLASH code.
<img src="Data/Figure1.png" alt="text" width="800"/>


The openPMD output after 6 level refinement in both X and Y directions.
<img src="src/output/visit0001.png" alt="text" width="400"/>


PIConGPU simulation by using the profile converted from FLASH. The figure shows the electron density profile duing the interaction with the rising edge of the laser main pulse.
<img src="src/output/dens_H_1384.png" alt="text" width="400"/> <img src="src/output/Ex_1384.png" alt="text" width="400"/>

#### Cylindrical geometry:

The output from cylindrical geometry simulation is rotated by 360 degree around its axis of symmetry to obtain 3D Cartesian geometry array.

<img src="Data/water_hdf5_plt_cnt_0050_Slice_theta_density.png" alt="text" width="300"/>

<img src="Data/3DRotated.png" alt="text" width="300"/>
