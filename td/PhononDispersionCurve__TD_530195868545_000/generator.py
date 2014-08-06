import json
from ase.data import chemical_symbols
from ase.data import reference_states
from random import randint

#get fcc stuff
chemical_symbols = [ sym for pk,sym in enumerate(chemical_symbols) if reference_states[pk] and reference_states[pk].get('symmetry') == 'fcc' ]

with open("test_generator.json","w") as f:
    for lattice in ['fcc']:
        for element in chemical_symbols:
            kimnum = "{:012d}".format(randint(0,10**12-1))
            f.write( json.dumps( {"symbol": element, "lattice": lattice, "kimnum": kimnum} ) +"\n")

