#! /bin/bash

soft=epock

#########################################################

echo "Esperando a " $soft

nproc=1
while [ $nproc -eq 1 ]
do
    nproc=$(ps aux | grep $soft | wc | awk '{print $1}')
done

echo "Arrancó " $soft

start_at=$(date +%s,%N)
_s1=$(echo $start_at | cut -d',' -f1)   # sec
_s2=$(echo $start_at | cut -d',' -f2)   # nano sec

nproc=2
while [ $nproc -eq 2 ]
do
    nproc=$(ps aux | grep epock | wc | awk '{print $1}')
done

end_at=$(date +%s,%N)
_e1=$(echo $end_at | cut -d',' -f1)
_e2=$(echo $end_at | cut -d',' -f2)
time_cost=$(bc <<< "scale=3; $_e1 - $_s1 + ($_e2 -$_s2)/1000000000")

echo $time_cost

exit 0
