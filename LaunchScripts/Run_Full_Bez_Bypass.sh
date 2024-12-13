scenario="BezAM"
subscenario="Single"
t=5000
exptype="SuperTest"
expnum=0
for spacing in {1..20..1}
    do
    ((expnum++))
    python3 ../AircraftSpacing.py -s1 $scenario -s2 $subscenario -d $spacing -t $t -en $expnum -et $exptype
done
