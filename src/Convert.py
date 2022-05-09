import sys
import os
sys.path.append(rf"C:\Users\WIZZARD100\F2OPMD\flash2openpmd\src")
import yt
import numpy as np

from flash2openpmd import Convert


run_directory = r'C:\Users\WIZZARD100\F2OPMD\flash2openpmd\Data'  #path of the flash2openpmd script you have on your PC

#filename = "lasslab_hdf5_plt_cnt_0043"
filename=input("Enter the name of the file you want to convert: ")
ts = Convert(run_directory,filename)

density = ts.get_data(x_min=-3e-4,y_min=0) #minimum valueas for x and y

ts.write2openpmd(density)