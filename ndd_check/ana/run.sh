#! /bin/bash
#/home/german/labo/16/ANA/Debug/ANA ../1prn.pdb -d ../1prn.nc -c 1_1prn.cfg
#/home/german/labo/16/ANA/Debug/ANA ../1prn.pdb -d ../1prn.nc -c 1_1prn.cfg -o salida_1_1prn > volumen_hiprec 
/home/german/labo/16/ANA/Debug/ANA ../pdbs/1prn_avg.pdb -c 1_1prn.cfg --NDD_input=in_ndd --NDD_output=out_ndd
exit 0
