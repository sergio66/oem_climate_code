hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
vchan2834 = hdfread(hdffile,'freq');
f = vchan2834;
load sarta_chans_for_l1c.mat
fairs = f(ichan);

xyz = load('f1305.mat');
fcris = xyz.f1305;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 1305
  woo = abs(fcris(ii)-fairs);
  woo = find(woo == min(woo));
  closest_airs2645_to_cris1305_chan(ii) = woo;
end

figure(1)
subplot(211); plot(fcris,fairs(closest_airs2645_to_cris1305_chan),'b.',fcris,fcris,'k')
  title('Closest airs2cris'); ylabel('fairs')
subplot(212); plot(fcris,fcris'-fairs(closest_airs2645_to_cris1305_chan),'b.')
  xlabel('fcris'); ylabel('delta')
delta_closest_airs2645_to_cris1305_chan = fcris'-fairs(closest_airs2645_to_cris1305_chan);

save closest_airs2645_to_cris1305_chan.mat *closest_airs2645_to_cris1305_chan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 2645
  woo = abs(fairs(ii)-fcris);
  woo = find(woo == min(woo));
  closest_cris1305_to_airs2645_chan(ii) = woo;
end

figure(2)
subplot(211); plot(fairs,fcris(closest_cris1305_to_airs2645_chan),'b.',fairs,fairs,'k')
  title('Closest cris2airs'); ylabel('fcris')
subplot(212); plot(fairs,fairs'-fcris(closest_cris1305_to_airs2645_chan),'b.')
  xlabel('fcris'); ylabel('delta')
delta_closest_cris1305_to_airs2645_chan = fairs'-fcris(closest_cris1305_to_airs2645_chan);

whos *closest_airs2645_to_cris1305_chan *closest_cris1305_to_airs2645_chan
comment = 'see closest_airs2645_to_cris1305.m';
save closest_airs2645_to_cris1305_chan.mat *closest_airs2645_to_cris1305_chan *closest_cris1305_to_airs2645_chan comment
