import pytraj as pt
traj = pt.load("../1prn.nc", top="../dry_1prn.prmtop")

for x in range(0, 599):

    pdb="snapshot"
    pdb+=str(x)
    pdb+=".pdb"
    pt.write_traj(pdb, traj, frame_indices=[x], overwrite=False)
