scenario="BezAM"
subscenario="NotSingle"
t=17000
exptype="SuperTest"
expnum=0
for (( spacing=311; spacing>=8; spacing-=20 ))
    do
    ((expnum++))
    python3 ../AircraftSpacing.py -s1 $scenario -s2 $subscenario -d $spacing -t $t -en $expnum -et $exptype
done
