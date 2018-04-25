

for files in $(cat raw_5p_rc/raw_5p_rc_list_2d.txt)
do
	f1=$(echo "$files" | awk -F '\t' '{print $1}')
	f2=$(echo "$files" | awk -F '\t' '{print $2}')
	ct=$(echo "$f1" | awk -F '.' '{print $1}')
	time python $script_folder'pknorm_0326.py' -r 'raw_5p_rc/'$f1 -t 'raw_5p_rc/'$f2 -m 1 -i 2 -f 0.05 -n 100000 -l 10000 -a 100000 -b 0 -s $script_folder -p nb
	time python $script_folder'pknorm_check_dif.py' -r $f1 -t $ct'.pknorm.txt' -f 0.05 -s $script_folder -p nb
done


