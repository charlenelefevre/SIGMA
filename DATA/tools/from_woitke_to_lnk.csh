awk '$1!~"data"{print $0}' Mg01Fe09S_Begemann1994.dat | tr '&' "\n" | tr -d " " | tr ',' "\n" | sed '/^$/d' | tr -d "/" | split -l 70
echo 70 2.16 > Mg01Fe09S_Begemann1994.dat
paste xaa xab > temp
paste temp xac >> Mg01Fe09S_Begemann1994.dat
