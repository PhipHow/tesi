#!/bin/csh

set samplegeo = 'wigner.mc_geo01'
set mopacgeo = 'trans_wigner.out'

rm nm1011.out
set ntraj=`head -1 $samplegeo | awk '{print $1}'`
grep 'Geometry number' $samplegeo >geom.out
set itraj=0
iter:
set itraj=`echo "$itraj + 1" |bc -l`
set nwig=`head -$itraj geom.out | tail -1 |awk '{print $4}'`
set r10=`grep -B50 " $nwig XXX" $mopacgeo | grep -A15 'dimensionless' | grep ' 10 '`  
set r11=`grep -B50 " $nwig XXX" $mopacgeo | grep -A15 'dimensionless' | grep ' 11 '` 
echo $itraj $nwig $r10 $r11 >>nm1011.out 
if ($itraj < $ntraj) goto iter
