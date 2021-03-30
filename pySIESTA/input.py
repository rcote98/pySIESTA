
import os

SIESTA_EXEC = os.getenv("SIESTA_EXEC", default = None) 
SIESTA_CORES = 4

class FDF():

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

    def launch(self, output_file=None):

        """
        Execute a SCALE-UP simulation with the current FDF settings.

        Parameters:
        ----------
        - output_file: human output filename. Defaults to [system_name].out.

        """

        #if not os.path.exists(self.scup_exec):
        #    print("WARNING: SCUP executable provided does not exist.")
        #    exit

        if output_file == None:
            # set default value after we know settings are loaded
            output_file = self.settings["System_name"] + ".out"

        # create temporary input
        self.read_fdf("_pySIESTAinput.fdf")

        command = f"mpirun -n {SIESTA_CORES} {SIESTA_EXEC} < _pySIESTAinput.fdf > {output_file}"

        # execute simulation
        print(command)
        #os.system(command)

        # remove temporary input
        os.remove("_pySIESTAinput.fdf")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #