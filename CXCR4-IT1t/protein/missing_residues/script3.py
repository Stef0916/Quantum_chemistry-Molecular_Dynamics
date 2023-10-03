from modeller import *
from modeller.automodel import *    # Load the AutoModel class

log.verbose()
env = Environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

class MyModel(AutoModel):
    def select_atoms(self):
        return Selection(self.residue_range('1:A', '44:A'),
                         self.residue_range('237:A', '249:A'), 
			 self.residue_range('355:A', '356:A'),
			 self.residue_range('400:A', '408:A'),
			 self.residue_range('487:A', '508:A'),)

a = MyModel(env, alnfile = 'alignment.ali',
            knowns = '3oe6', sequence = '3oe6_fill')
a.starting_model= 1
a.ending_model  = 1

a.make()