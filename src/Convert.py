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