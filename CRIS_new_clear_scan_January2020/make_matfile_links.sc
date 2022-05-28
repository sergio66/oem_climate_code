# 12:lrwxrwxrwx 1 sergio pi_strow   49 Jan 16  2020 kbad_chans0.mat -> ../AIRS_new_clear_scan_August2019/kbad_chans0.mat
# 13:lrwxrwxrwx 1 sergio pi_strow   53 Jan 16  2020 good_chans_2016.mat -> ../AIRS_new_clear_scan_August2019/good_chans_2016.mat
# 14:lrwxrwxrwx 1 sergio pi_strow   59 Jan 16  2020 indices_of_l1b_in_l1c.mat -> ../AIRS_new_clear_scan_August2019/indices_of_l1b_in_l1c.mat
# 15:lrwxrwxrwx 1 sergio pi_strow   45 Jan 16  2020 stratSW.mat -> ../AIRS_new_clear_scan_August2019/stratSW.mat
# 16:lrwxrwxrwx 1 sergio pi_strow   57 Jan 16  2020 sarta_chans_for_l1c.mat -> ../AIRS_new_clear_scan_August2019/sarta_chans_for_l1c.mat
# 17:lrwxrwxrwx 1 sergio pi_strow   45 Jan 16  2020 btn_avg.mat -> ../AIRS_new_clear_scan_August2019/btn_avg.mat

rm kbad_chans0.mat good_chans_2016.mat indices_of_l1b_in_l1c.mat stratSW.mat sarta_chans_for_l1c.mat btn_avg.mat
ln -s ../../AIRS_new_clear_scan_August2019_AMT2020PAPER/kbad_chans0.mat .
ln -s ../../AIRS_new_clear_scan_August2019_AMT2020PAPER/good_chans_2016.mat .
ln -s ../../AIRS_new_clear_scan_August2019_AMT2020PAPER/indices_of_l1b_in_l1c.mat .
ln -s ../../AIRS_new_clear_scan_August2019_AMT2020PAPER/stratSW.mat .
ln -s ../../AIRS_new_clear_scan_August2019_AMT2020PAPER/sarta_chans_for_l1c.mat .
ln -s ../../AIRS_new_clear_scan_August2019_AMT2020PAPER/btn_avg.mat .
