function olr = compute_olr(hoem,poem);

addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD

poemx = do_clear_cloud_calc(hoem,poem);
olr = poemx.sarta_rclearcalc;


