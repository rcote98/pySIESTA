
import numpy as np
import pandas as pd

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def generate_lines_that_contain(string, fp):
    for line in fp:
        if string in line:
            yield line.rstrip()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def parse_energy(fname):

    energy = {}

    with open(fname, "r") as f:

        en_parse = False
        counter = 12
        for line in generate_lines_that_contain("siesta:", f):
            
            if "Final energy (eV):" in line:
                en_parse = True

            if en_parse and counter > 0 and not ("Final energy (eV):" in line):
                
                line = remove_prefix(line, "siesta:")
                line = line.split("=")
                line = [w.replace(" ", "") for w in line]
                line[1] = float(line[1])
                energy[line[0]] = float(line[1])

                counter -= 1

    return energy

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #