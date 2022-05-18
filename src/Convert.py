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

from flash2openpmd import Convert


run_directory = rf"~/Flash2OpenPMD/Data"  #path of the flash2openpmd script you have on your PC

#filename = "lasslab_hdf5_plt_cnt_0043"
filename = input("Enter the name of the file you want to convert: ")
ts = Convert(run_directory,filename)

density = ts.get_data(x_min=0,y_min=0) #minimum valueas for x and y

ts.write2openpmd(density)