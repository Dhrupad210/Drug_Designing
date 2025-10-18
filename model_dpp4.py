# Import MODELLER
from modeller import *
from modeller.automodel import *

# Initialize the environment
env = environ()
env.io.atom_files_directory = ['.']  # Current directory with PDB

# Create an automodel class
a = automodel(env,
              alnfile  = 'align.pir',   # Your PIR alignment
              knowns   = 'ij2e_renum',      # Template code from PIR
              sequence = 'dpp4',        # Query code from PIR
              assess_methods=(assess.DOPE, assess.GA341))

# Generate 5 models (you can increase if needed)
a.starting_model = 1
a.ending_model   = 3

# Run the modeling
a.make()
