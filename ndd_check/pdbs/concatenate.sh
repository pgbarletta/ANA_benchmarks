#!/bin/bash

ls *prn.pdb | sort -n | xargs -i cat {} > 1prn_traj.pdb
sed -i '/END$/d' 1prn_traj.pdb 
echo "END" >> 1prn_traj.pdb
exit 0
