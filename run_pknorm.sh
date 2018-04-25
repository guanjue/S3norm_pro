script_folder='/storage/home/gzx103/group/software/PKnorm/pknorm_scripts/'

time python $script_folder'pknorm_0326.py' -r reference_signal_track.txt -t target_signal_track.txt -m 1 -i 2 -f 0.05 -n 5000 -l 1000 -a 100 -b 0 -s $script_folder -p nb


time python $script_folder'pknorm_0326.py' -r ERY_fl.ctcfrep.100162.bamtobed5endintersect.signal -t ERY_fl.ctcfrep.100163.bamtobed5endintersect.signal -m 1 -i 2 -f 0.05 -n 100000 -l 10000 -a 100000 -b 0 -s $script_folder -p z

time python $script_folder'pknorm_0326.py' -r LSK_BM_987.atacrep.987.bamtobed5endintersect.signal -t LSK_BM_1196.atacrep.1196.bamtobed5endintersect.signal -m 1 -i 2 -f 0.05 -n 100000 -l 10000 -a 100000 -b 0 -s $script_folder -p z


time python $script_folder'pknorm_check_dif.py' -r LSK_BM_987.atacrep.987.bamtobed5endintersect.signal -t LSK_BM_1196.atacrep.1196.bamtobed5endintersect.signal -f 0.05 -s $script_folder -p z




