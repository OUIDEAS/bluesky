scenario="Hold"
subscenario="Single"
t=5000
exptype="SingleHold"
expnum=1
for (( spacing=160; spacing>=8; spacing-=8 ))
    do
    ((expnum++))
    python3 ../AircraftSpacing.py -s1 $scenario -s2 $subscenario -d $spacing -t $t -en $expnum -et $exptype
done
