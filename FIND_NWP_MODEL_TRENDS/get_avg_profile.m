[hh2,hha2,pp20,ppa2] = rtpread('/home/sergio/KCARTA/WORK/feb2002_raw_op_airs.rad.constemiss.rtp');
junk2378 = hh2.vchan;
[hh2,hha2,pp20,ppa2] = rtpread('/home/sergio/KCARTA/WORK/pin_feb2002_sea_airsnadir_op.so2.latlon.const_emiss.rtp');
hh2.vchan = junk2378;

hh2 = h; 
if h.nchan ~= 2645
  error('ohoh need 2645 chans')
end
h.pfields = 1;

[hh2,pp2]   = subset_rtp(hh2,pp20,[],[],[1 2 4]); %% TRP, MLS, SAS
[hh2,pp2_1] = subset_rtp(hh2,pp20,[],[],1); %% TRP, MLS, SAS
[hh2,pp2_2] = subset_rtp(hh2,pp20,[],[],2); %% TRP, MLS, SAS
[hh2,pp2_4] = subset_rtp(hh2,pp20,[],[],4); %% TRP, MLS, SAS

sarta = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';
fop = 'junk.op.rtp';
frp = 'junk.rp.rtp';
sartaer = ['!' sarta ' fin=' fop ' fout=' frp ' >& ugh'];

pp0 = pp2_1;
[hh2,pp0] = cat_rtp(hh2,pp0,hh2,pp2_2);
[hh2,pp0] = cat_rtp(hh2,pp0,hh2,pp2_2);
[hh2,pp0] = cat_rtp(hh2,pp0,hh2,pp2_4);
[hh2,pp0] = cat_rtp(hh2,pp0,hh2,pp2_4);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% ocean
junk = find(p.landfrac == 0 & abs(p.rlat) < 30 & p.nemis == 19);
  nemis = floor(mean(p.nemis(junk)));
  efreq = mean(p.efreq(1:nemis,junk),2);
  emis = mean(p.emis(1:nemis,junk),2);

pp0.nemis = nemis * ones(size(pp0.stemp));
pp0.emis = emis * ones(size(pp0.stemp));
pp0.efreq = efreq * ones(size(pp0.stemp));
pp0.rho = (1-pp0.emis)/pi;

rtpwrite(fop,hh2,hha2,pp0,ppa2);
eval(sartaer);

[~,~,pprad0_ocean,~] = rtpread(frp);

disp('STEMP and MMW for TRP, NML,SML, NP,SP')
fprintf(1,' %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',[pprad0_ocean.stemp; mmwater_rtp(hh2,pprad0_ocean)])

%%%%%%%%%%%%%%%%%%%%%%%%%
%% land
pp0 = rmfield(pp0,'emis');
pp0 = rmfield(pp0,'efreq');
pp0 = rmfield(pp0,'rho');

pp0.emis = zeros(100,5);
pp0.efreq = zeros(100,5);
pp0.rho = zeros(100,5);

%% land
junk = find(p.landfrac == 1 & abs(p.rlat) < 30 & p.nemis > 30);
  nemis = floor(mean(p.nemis(junk)));
  efreq = mean(p.efreq(1:nemis,junk),2);
  emis = mean(p.emis(1:nemis,junk),2);
  [Y,I] = sort(efreq); efreq = efreq(I); emis = emis(I);
  pp0.nemis(1) = nemis;
  pp0.emis(1:nemis,1) = emis;
  pp0.efreq(1:nemis,1) = efreq;
  pp0.rho = (1-pp0.emis)/pi;
junk = find(p.landfrac == 1 & p.rlat > -60 & p.rlat < -30 & p.nemis > 30);
  nemis = floor(mean(p.nemis(junk)));
  efreq = mean(p.efreq(1:nemis,junk),2);
  emis = mean(p.emis(1:nemis,junk),2);
  [Y,I] = sort(efreq); efreq = efreq(I); emis = emis(I);
  pp0.nemis(2) = nemis;
  pp0.emis(1:nemis,2) = emis;
  pp0.efreq(1:nemis,2) = efreq;
  pp0.rho = (1-pp0.emis)/pi;
junk = find(p.landfrac == 1 & p.rlat > +30 & p.rlat < +60 & p.nemis > 30);
  nemis = floor(mean(p.nemis(junk)));
  efreq = mean(p.efreq(1:nemis,junk),2);
  emis = mean(p.emis(1:nemis,junk),2);
  [Y,I] = sort(efreq); efreq = efreq(I); emis = emis(I);
  pp0.nemis(3) = nemis;
  pp0.emis(1:nemis,3) = emis;
  pp0.efreq(1:nemis,3) = efreq;
  pp0.rho = (1-pp0.emis)/pi;
junk = find(p.landfrac == 1 & p.rlat < -60 & p.nemis > 30);
  nemis = floor(mean(p.nemis(junk)));
  efreq = mean(p.efreq(1:nemis,junk),2);
  emis = mean(p.emis(1:nemis,junk),2);
  [Y,I] = sort(efreq); efreq = efreq(I); emis = emis(I);
  pp0.nemis(4) = nemis;
  pp0.emis(1:nemis,4) = emis;
  pp0.efreq(1:nemis,4) = efreq;
  pp0.rho = (1-pp0.emis)/pi;
junk = find(p.landfrac == 1 & p.rlat > +60 & p.nemis > 30);
  nemis = floor(mean(p.nemis(junk)));
  efreq = mean(p.efreq(1:nemis,junk),2);
  emis = mean(p.emis(1:nemis,junk),2);
  [Y,I] = sort(efreq); efreq = efreq(I); emis = emis(I);
  pp0.nemis(5) = nemis;
  pp0.emis(1:nemis,5) = emis;
  pp0.efreq(1:nemis,5) = efreq;
  pp0.rho = (1-pp0.emis)/pi;

rtpwrite(fop,hh2,hha2,pp0,ppa2);
eval(sartaer);

[~,~,pprad0_land,~] = rtpread(frp);

