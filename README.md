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

The converted file will be located in the output folder with the extension `.h5`

Simulation output from FLASH code.
<img src="Data/Figure1.png" alt="text" width="800"/>


The openPMD output after 6 times refinement in both X and Y directions.
<img src="src/output/visit0001.png" alt="text" width="400"/>


PIConGPU simulation by using the profile converted from FLASH. The figure shows the electron density profile duing the interaction with the rising edge of the laser main pulse.
<img src="src/output/dens_H_1384.png" alt="text" width="400"/> <img src="src/output/Ex_1384.png" alt="text" width="400"/>
