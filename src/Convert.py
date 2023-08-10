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

import sys
import os
import pwd
user = pwd.getpwuid(os.getuid())[0]
sys.path.append(rf"{user}/Flash2OpenPMD/src")
import yt
import numpy as np
import glob
from flash2openpmd import Convert


run_directory = input("Enter path of the simulation data:")  # path of the simulation data

#filename = "lasslab_hdf5_plt_cnt_0043"
        
step = input("Iteration: ")

filename = glob.glob1(run_directory,"*_hdf5_plt_cnt_"+str(int(step)).zfill(4))

print(filename[0])

ts = Convert(run_directory,filename[0])

density = ts.get_data() #minimum valueas for x and y

ts.write2openpmd(density)