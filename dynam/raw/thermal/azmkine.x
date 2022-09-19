#!/bin/csh

set photodyn = '/home/carlotta/azometano/dynam/thermal/400traj_boltz/photo.dyn'

/home/mau/azometano/dynam/azmkine << EOF > azmkine.out
&DAT tmin=-1.0, tmax=600.0, filedyn="$photodyn" &END
EOF

/home/mau/azometano/dynam/azmkine << EOF >> azmkine.out
&DAT tmin=5000.0, tmax=20000.0, filedyn="$photodyn" &END
EOF

/home/mau/azometano/dynam/azmkine << EOF >> azmkine.out
&DAT tmin=-1.0, tmax=20000.0, filedyn="$photodyn" &END
EOF

set ntraj = `tail -9 azmkine.out | grep 'Total' | awk '{print $5}'`
set nbad = `tail -9 azmkine.out | grep 'bad traj' | awk '{print $5}'`
set ndiss = `tail -9 azmkine.out | grep 'dissoc' | awk '{print $5}'`
set ekaz = `tail -9 azmkine.out | grep 'N2' | awk '{print $6}'`
set ekfir = `tail -9 azmkine.out | grep 'first CH3' | awk '{print $7}'`
set eksec = `tail -9 azmkine.out | grep 'second CH3' | awk '{print $7}'`
set ekazlast = `grep -A5 '5000.00    20000.00' azmkine.out | \
   grep 'N2' | awk '{print $6}'`
set ekfirlast = `grep -A5 '5000.00    20000.00' azmkine.out | \
   grep 'first CH3' | awk '{print $7}'`
set ekseclast = `grep -A5 '5000.00    20000.00' azmkine.out | \
   grep 'second CH3' | awk '{print $7}'`

set ntraj = `echo "$ntraj - $nbad" | bc -l`
set nodis = `echo "$ntraj - $ndiss" | bc -l`
set ekaz = \
   `echo "scale=6; ($ekaz *$ndiss +$ekazlast *$nodis)/$ntraj" | bc -l`
set ekfir = \
   `echo "scale=6; ($ekfir *$ndiss +$ekfirlast *$nodis)/$ntraj" | bc -l`
set eksec = \
   `echo "scale=6; ($eksec *$ndiss +$ekseclast *$nodis)/$ntraj" | bc -l`

echo ' ' >> azmkine.out
echo 'Extrapolated values' >> azmkine.out
echo 'Number of undissociated trajectories:' $nodis >> azmkine.out
echo 'Average transl. energy of N2          ' $ekaz 'eV' >> azmkine.out
echo 'Average transl. energy of first CH3   ' $ekfir 'eV' >> azmkine.out
echo 'Average transl. energy of second CH3  ' $eksec 'eV' >> azmkine.out
