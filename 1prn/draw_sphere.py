from pymol.cgo import *
from pymol import cmd

spherelist = [
   COLOR,    1.100,    1.000,    1.000,
   SPHERE,   38.0,   41.0,   41.0,    0.80,
    ]

cmd.load_cgo(spherelist, 'point',   1)
