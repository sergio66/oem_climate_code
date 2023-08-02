## Howard said it is fine to use 36 threads .... though if I launch 8 jpbs should that be 36/8 = 4.5 threads???
## so lets use 8 threads

matlab -nodesktop -nodisplay -singleCompThread -r "LASTN = maxNumCompThreads(8); JOBM=1; local_do_quick_8jobs; quit "  > out1 &
matlab -nodesktop -nodisplay -singleCompThread -r "LASTN = maxNumCompThreads(8); JOBM=2; local_do_quick_8jobs; quit "  > out2 &
matlab -nodesktop -nodisplay -singleCompThread -r "LASTN = maxNumCompThreads(8); JOBM=3; local_do_quick_8jobs; quit "  > out3 &
matlab -nodesktop -nodisplay -singleCompThread -r "LASTN = maxNumCompThreads(8); JOBM=4; local_do_quick_8jobs; quit "  > out4 &
matlab -nodesktop -nodisplay -singleCompThread -r "LASTN = maxNumCompThreads(8); JOBM=5; local_do_quick_8jobs; quit "  > out5 &
matlab -nodesktop -nodisplay -singleCompThread -r "LASTN = maxNumCompThreads(8); JOBM=6; local_do_quick_8jobs; quit "  > out6 &
matlab -nodesktop -nodisplay -singleCompThread -r "LASTN = maxNumCompThreads(8); JOBM=7; local_do_quick_8jobs; quit "  > out7 &
matlab -nodesktop -nodisplay -singleCompThread -r "LASTN = maxNumCompThreads(8); JOBM=8; local_do_quick_8jobs; quit "  > out8 &
