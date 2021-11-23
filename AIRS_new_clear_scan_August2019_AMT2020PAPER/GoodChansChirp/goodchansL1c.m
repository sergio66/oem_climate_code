addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CLOUD

%% simple code to translate my index list of L1b goodchannels to L1c indices

fairs2378 = instr_chans;

g2378 = dogoodchan;

load ../indices_of_l1b_in_l1c.mat
load ../f2645

clf; plot(fairs2378,'b'); hold on; plot(f2645(l1b_ind_in_l1c),'r'); hold off; title('DEFINITELY WRONG')
clf; plot(fairs2378,'b'); hold on; plot(f2645(l1c_ind_for_l1b),'r'); hold off; title('CORRECT map 2645 to 2378')

for ii = 1 : length(g2378)
  junk = abs(fairs2378(g2378(ii))-f2645);
  g2645(ii) = find(junk == min(junk),1);
end;

clf; plot(fairs2378(g2378),'bo'); hold; plot(f2645(g2645),'r'); hold off

plot(fairs2378(g2378),fairs2378(g2378)-f2645(g2645),'r'); hold off

bad = find(fairs2378(g2378)-f2645(g2645) > 0);
good = find(fairs2378(g2378)-f2645(g2645) < 0);
g2645 = g2645(good);
clf; plot(fairs2378(g2378),'bo'); hold; plot(f2645(g2645),'r'); hold off

comment = 'list of L1b good channels changed to L1c; see ~/MATLABCODE/oem_pkg_run/AIRS_new_clear_scan_August2019/GoodChansChirp/goodchansL1c.m';
save good2645.mat g2645 comment
