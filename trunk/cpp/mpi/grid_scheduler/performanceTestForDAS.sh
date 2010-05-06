#!/bin/sh
rm *.txt
rm *.out

r=15
t=100
p=1600
b=8
i=0
f=1
echo "#Resources: "$r ", Tasks: "$t ", Population: "$p ", Batch: "$b ", Initilization: "$i ", Fitness: "$f >> timeMeasure.out

prun -v -1 -sge-script /usr/local/sitedep/reserve.sge/sge-script "/home2/skoulouz/scheduler_10-12/simulation -r $r -t $t -p $p -b $b -i $i -f $f " 1 >> timeMeasure.out
# tail timeMeasure.out
prun -v -1 -sge-script /usr/local/sitedep/reserve.sge/sge-script "/home2/skoulouz/scheduler_10-12/simulation -r $r -t $t -p $p -b $b -i $i -f $f " 2 >> timeMeasure.out
# tail timeMeasure.out
prun -v -1 -sge-script /usr/local/sitedep/reserve.sge/sge-script "/home2/skoulouz/scheduler_10-12/simulation -r $r -t $t -p $p -b $b -i $i -f $f " 4 >> timeMeasure.out
# tail timeMeasure.out
prun -v -1 -sge-script /usr/local/sitedep/reserve.sge/sge-script "/home2/skoulouz/scheduler_10-12/simulation -r $r -t $t -p $p -b $b -i $i -f $f " 8 >> timeMeasure.out
# tail timeMeasure.out
prun -v -1 -sge-script /usr/local/sitedep/reserve.sge/sge-script "/home2/skoulouz/scheduler_10-12/simulation -r $r -t $t -p $p -b $b -i $i -f $f " 16 >> timeMeasure.out
# tail timeMeasure.out

echo "Done!!"
