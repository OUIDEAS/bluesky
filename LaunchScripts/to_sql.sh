scenario="BezAM"
exptype="SuperTest"
for (( spacing=311; spacing>=8; spacing-=20 ))
    do
    python3 ../SQL/json_to_sql.py -s1 $scenario -d $spacing -et $exptype
done
