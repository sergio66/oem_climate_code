%-rw-rw-r-- 1 sergio pi_strow 16639 May 11 22:25 smo_F12_MM.dat
%-rw-rw-r-- 1 sergio pi_strow 14683 May 11 22:25 sum_F12_MM.dat
%-rw-rw-r-- 1 sergio pi_strow 16642 May 11 22:24 nwr_F12_MM.dat
%-rw-rw-r-- 1 sergio pi_strow 16631 May 11 22:23 mlo_F12_MM.dat
%-rw-rw-r-- 1 sergio pi_strow 16626 May 11 22:23 brw_F12_MM.dat
%-rw-rw-r-- 1 sergio pi_strow 16640 Mar  8 16:47 spo_F12_MM.dat

addpath /home/sergio/MATLABCODE/FIND_TRENDS

junk = load('CFC12/brw_F12_MM.dat');
  yy = junk(:,1) + (junk(:,2)-1)/12;
  boo = find(yy >= 2002+9/12 & yy <= 2019+8/12 & isfinite(junk(:,2))); 
  plot(yy(boo),junk(boo,3))
  doy = yy-2002; doy = doy*365.25; B = Math_tsfit_lin_robust(doy(boo),junk(boo,3),4); cfc12(1) = B(2);

junk = load('CFC12/mlo_F12_MM.dat');
  yy = junk(:,1) + (junk(:,2)-1)/12;
  boo = find(yy >= 2002+9/12 & yy <= 2019+8/12 & isfinite(junk(:,2))); 
  plot(yy(boo),junk(boo,3))
  doy = yy-2002; doy = doy*365.25; B = Math_tsfit_lin_robust(doy(boo),junk(boo,3),4); cfc12(2) = B(2);

junk = load('CFC12/nwr_F12_MM.dat');
  yy = junk(:,1) + (junk(:,2)-1)/12;
  boo = find(yy >= 2002+9/12 & yy <= 2019+8/12 & isfinite(junk(:,2))); 
  plot(yy(boo),junk(boo,3))
  doy = yy-2002; doy = doy*365.25; B = Math_tsfit_lin_robust(doy(boo),junk(boo,3),4); cfc12(3) = B(2);

junk = load('CFC12/smo_F12_MM.dat');
  yy = junk(:,1) + (junk(:,2)-1)/12;
  boo = find(yy >= 2002+9/12 & yy <= 2019+8/12 & isfinite(junk(:,2))); 
  plot(yy(boo),junk(boo,3))
  doy = yy-2002; doy = doy*365.25; B = Math_tsfit_lin_robust(doy(boo),junk(boo,3),4); cfc12(4) = B(2);

junk = load('CFC12/spo_F12_MM.dat');
  yy = junk(:,1) + (junk(:,2)-1)/12;
  boo = find(yy >= 2002+9/12 & yy <= 2019+8/12 & isfinite(junk(:,2))); 
  plot(yy(boo),junk(boo,3))
  doy = yy-2002; doy = doy*365.25; B = Math_tsfit_lin_robust(doy(boo),junk(boo,3),4); cfc12(5) = B(2);

junk = load('CFC12/sum_F12_MM.dat');
  yy = junk(:,1) + (junk(:,2)-1)/12;
  boo = find(yy >= 2002+9/12 & yy <= 2019+8/12 & isfinite(junk(:,2))); 
  plot(yy(boo),junk(boo,3))
  doy = yy-2002; doy = doy*365.25; B = Math_tsfit_lin_robust(doy(boo),junk(boo,3),4); cfc12(6) = B(2);

