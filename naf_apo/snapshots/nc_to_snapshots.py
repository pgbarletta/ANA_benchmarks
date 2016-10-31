import pytraj as pt
traj = pt.load("../../naf_apo.nc", top="../../dry_naf_apo_154.prmtop")

for x in range(0, 599):

    pdb="snapshot"
    pdb+=str(x)
    pdb+=".pdb"
    pt.write_traj(pdb, traj, frame_indices=[x], overwrite=False)

