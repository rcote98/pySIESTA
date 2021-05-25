"""
Classes execute SIESTA simulations
"""

# third party imports
import numpy as np          # matrix support
import pandas as pd         # .out file loading

# standard library imports
from shutil import move,rmtree,copy # remove output folder
from pathlib import Path            # general folder management
import os, sys, csv                 # remove files, get pwd
import pickle                       # store parameter vectors
import time                         # check simulation run time
import re                           # regular expressions

# package imports
from pySIESTA.structures import SolidGeometry

SIESTA_CORES = 20
SIESTA_EXEC = os.getenv("SIESTA_EXEC", default = "None") 

if not os.path.exists(SIESTA_EXEC):
    print("WARNING: SIESTA executable provided does not exist.")
    #exit

def log(log_file, message):

    with open(log_file, "a+") as f:
        # Move read cursor to the start of file.
        f.seek(0)
        # If file is not empty then append '\n'
        data = f.read(100)
        if len(data) > 0 :
            f.write("\n")
        # Append text at the end of file
        f.write(message)


class FDFSettings():

    """
        docs
    """

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def reset(self):

        self.settings = {}
        pass

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def read_fdf(self, fdf_file):


        # empty previous settings
        self.reset()

        with open(fdf_file, "r") as f:

            lines = f.readlines()

            # do a general cleanup
            lines = [l.strip() for l in lines]          # remove trailing spaces
            lines = [l.split("#",1)[0] for l in lines]  # removes comments
            lines = [l for l in lines if l]             # remove empty lines
            lines = [l.lower() for l in lines]          # make everything lowercase
            lines = [l.split() for l in lines]          # split lines into lists

        i = 0
        block = None
        blockname = None

        while i < len(lines): 
            
            line = lines[i]
            
            # check if its a block setting
            if line[0] == r"%block":

                block = []
                blockname = line[1]
                
                i+=1
                line = lines[i]

                while line[0] != r"%endblock":
                    block.append(line)
                    i += 1
                    line = lines[i]

                self.settings[blockname] = block

            # if its not, then just read it
            else:
                if len(line) == 1:
                    # setting without value = true
                    self.settings[line[0]] = [".true."]
                else:
                    self.settings[line[0]] = line[1:]

            i += 1
        
        pass

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def write_fdf(self, fname):

        if len(self.settings) == None:
            print("WARNING: No settings file loaded. Aborting.")
            return 0

        f = open(fname, "w")

        for k in self.settings:

            if isinstance(self.settings[k][0], list):

                f.write(r"%block " + k + "\n")

                # turn array into string
                string = ""
                for row in self.settings[k]:
                    for n in row:
                        string += str(n) + " "
                    string += "\n"

                f.write(string)
                f.write(r"%endblock " + k + "\n")
            
            else:
                
                line = k + " " + " ".join([str(el) for el in self.settings[k]]) + "\n" 
                f.write(line)

            f.write("\n")
                
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                         simple jobs                           #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def single_point_run(folder, geom, fdf):

    original_dir = os.getcwd()
    working_dir = os.path.abspath(folder)

    try:
        os.makedirs(working_dir)
    except:
        pass

    # create input geometry file
    geom_file = os.path.join(working_dir, fdf.settings["systemlabel"] + ".STRUCT_IN")
    geom.write_STRUCT(geom_file)

    # create input file
    fdf_file  = os.path.join(working_dir, "pySIESTAinput.fdf")
    fdf.settings['md.numcgsteps'] = 0
    fdf.write_fdf(fdf_file)

    # move to directory 
    os.chdir(working_dir)

    # run the command
    command = f"mpirun -n {SIESTA_CORES} {SIESTA_EXEC} < pySIESTAinput.fdf > SIESTA_output.txt"
    os.system(command)

    # go back to 505
    os.chdir(original_dir)

    return True

def geometry_run(folder, geom, fdf, max_steps=100):

    original_dir = os.getcwd()
    working_dir = os.path.abspath(folder)

    try:
        os.makedirs(working_dir)
    except:
        pass

    # create input geometry file
    geom_file = os.path.join(working_dir, fdf.settings["systemlabel"][0] + ".STRUCT_IN")
    geom.write_STRUCT(geom_file)

    # create input file
    fdf_file  = os.path.join(working_dir, "pySIESTAinput.fdf")
    fdf.settings['md.numcgsteps'] =  [max_steps]
    fdf.write_fdf(fdf_file)

    # move to directory 
    os.chdir(working_dir)

    # run the command
    command = f"\nmpirun -n {SIESTA_CORES} {SIESTA_EXEC} < pySIESTAinput.fdf > SIESTA_output.txt\n"
    os.system(command)

    # go back to 505
    os.chdir(original_dir)

    return True