ls -lt Output/Quantile05/test*.mat | wc -l
echo "removing latbins 01 .. 11 == first 11x72 = 792 mat files"
## 11 x 72 === 792
## 11 x 72 === 792
## 11 x 72 === 792
## remove test 1,test2, ... test9
/bin/rm Output/Quantile05/test[123456789].mat;  
## remove test 10,.19,test20..29, ... test90..99
/bin/rm Output/Quantile05/test[123456789][0123456789].mat
## remove test 100..199,200..299,300..399,400..499,500..599,600..699
/bin/rm Output/Quantile05/test[123456][0123456789][0123456789].mat
## remove test 700..789
/bin/rm Output/Quantile05/test7[012345678][0123456789].mat
## remove test 790..792
/bin/rm Output/Quantile05/test79[012].mat

ls -lt Output/Quantile05/test*.mat | wc -l
echo "removing latbins 54 .. 64 == last 11x72 = 792 mat files"
## 11 x 72 === 792
## 11 x 72 === 792
## 11 x 72 === 792
## remove test 3817 ..3819
/bin/rm Output/Quantile05/test381[789].mat
## remove test 3820:3899
/bin/rm Output/Quantile05/test38[23456789][0123456789].mat
## remove test 3900:3999
/bin/rm Output/Quantile05/test39[0123456789][0123456789].mat
## remove test 4000:4599
/bin/rm Output/Quantile05/test4[012345][0123456789][0123456789].mat
## remove test 4600 .. 4608
/bin/rm Output/Quantile05/test460[012345678].mat

ls -lt Output/Quantile05/test*.mat | wc -l

/bin/rm slurm*.out
## loop forwards 1 -- 72
sbatch  -p high_mem  --exclude=cnode018,cnode019 --array=01-11 sergio_matlab_jobB.sbatch 0
sbatch  -p high_mem  --exclude=cnode018,cnode019 --array=54-64 sergio_matlab_jobB.sbatch 0
## loop backwards 72 -- 1
sbatch  -p high_mem  --exclude=cnode018,cnode019 --array=01-11 sergio_matlab_jobB.sbatch 1
sbatch  -p high_mem  --exclude=cnode018,cnode019 --array=54-64 sergio_matlab_jobB.sbatch 1
