rm -f x*
rm -f temp.dat
split -l 5 Mg00Fe10S_Begemann1994.dat 
foreach i (x*)
awk '{print $3}' $i | tr -d '\r' | tr '\n' ', '| sed 's%^M% %g' | awk '{print "     & "$0}' > $i".dat"
sed 's%+00%+0%g' $i".dat" >> temp.dat
sed 's%-00%-0%g' temp.dat >> $i".dat"
end
