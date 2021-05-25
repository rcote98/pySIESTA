
import re

# auxiliary funcs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text 

def generate_lines_that_match(string, fp):
    for line in fp:
        if re.search(string, line):
            yield line

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def parse_energy(fname):

    energy = {}

    with open(fname, "r") as f:

        en_parse = False
        counter = 12
        for line in generate_lines_that_match("siesta:", f):
            
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


def check_last_cg_step(fname):
    
    lcg = None
    with open(fname, "r") as f:
        string = "Begin CG opt. move ="
        for l in generate_lines_that_match(string, f):
            lcg = l

    if lcg is None:
        return None
    else:
        return int(lcg.split()[5])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #