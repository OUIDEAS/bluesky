scenario="Hold"
exptype="SingleHold"
for (( spacing=160; spacing>=8; spacing-=8 ))
    do
    python3 ../SQL/json_to_sql.py -s1 $scenario -d $spacing -et $exptype
done
