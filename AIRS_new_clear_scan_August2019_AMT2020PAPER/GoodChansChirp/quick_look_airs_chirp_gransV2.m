%
% quick look at chirp granules
% copied from /home/motteler/shome/chirp_test/quick_look.m
%

% addpath /home/motteler/cris/ccast/source

p1 = '/asl/hpcnfs1/chirp/airs_L1c_src/2019/061';
% g1 = 'SNDR.SS1330.CHIRP.20190302T0011.m06.g002.L1_AIR.std.v01_07.U.2004132227.nc';
  g1 = 'SNDR.SS1330.CHIRP.20190302T0017.m06.g003.L1_AIR.std.v01_07.U.2004132229.nc';
  g1 = 'SNDR.SS1330.CHIRP.20190302T0211.m06.g022.L1_AIR.std.v01_07.U.2004132258.nc';

p2 = '/asl/hpcnfs1/chirp/cris_npp_src/2019/061';
g2 = 'SNDR.SS1330.CHIRP.20190302T0453.m06.g050.L1_CNP.std.v01_07.U.2004132235.nc';
g2 = 'SNDR.SS1330.CHIRP.20190302T0253.m06.g030.L1_CNP.std.v01_07.U.2004132231.nc';

[d1, a1] = read_netcdf_h5(fullfile(p1, g1));  % AIRS
[d2, a2] = read_netcdf_h5(fullfile(p2, g2));  % CrIS

bt1 = real(rad2bt(d1.wnum, d1.rad));
bt2 = real(rad2bt(d2.wnum, d2.rad));

figure(1); clf
subplot(2,1,1)
plot(d1.wnum, bt1(:, 2001:2020))
title('AIRS to CHIRP sample BT spectra')
ylabel('BT (K)')
grid on

subplot(2,1,2)
plot(d2.wnum, bt2(:, 2001:2020))
title('CrIS NPP to CHIRP sample BT spectra')
xlabel('wavenumber (cm-1)')
ylabel('BT (K)')
grid on

figure(2); clf
subplot(2,1,1)
[x1, y1] = pen_lift(d1.wnum, d1.nedn);
semilogy(x1, y1)
ylim([0.001, 1.0])
title('AIRS L1C to CHIRP NEdN')
ylabel('mw sr-1 m-2')
grid on

subplot(2,1,2)
[x2, y2] = pen_lift(d2.wnum, d2.nedn);
semilogy(x2, y2)
ylim([0.001, 1.0])
title('CrIS NPP to CHIRP NEdN')
xlabel('wavenumber (cm-1)')
ylabel('mw sr-1 m-2')
grid on

chirp.wnum = d1.wnum;
chirp.nednA = nanmean(d1.nedn,2);
chirp.nednC = nanmean(d2.nedn,2);

comment = 'see quick_look_airs_chirp_grans.m';
save chirpNeDN.mat chirp comment

addpath /home/sergio/MATLABCODE/PCA_NLTE
addpath /home/sergio/MATLABCODE/ClimateToolbox/EOF
%[eof_maps_ssw,pc_ssw,expvar_ssw] = eof_sergio(double(bt1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load "same" AIRS and CrIS granules for an NEdN comparison
apath = '/asl/hpcnfs1/airs/L1C/2019/061';
agran = 'AIRS.2019.03.02.005.L1C.AIRS_Rad.v6.1.2.0.G19061124436.hdf';
agran = 'AIRS.2019.03.02.022.L1C.AIRS_Rad.v6.1.2.0.G19061125315.hdf';
d3 = read_airs_h4(fullfile(apath, agran));

cpath = '/home/motteler/shome/daac_test/SNPPCrISL1B.2/2019/061';
cgran = 'SNDR.SNPP.CRIS.20190302T0054.m06.g010.L1B.std.v02_05.G.190302083725.nc';
cgran = 'SNDR.SNPP.CRIS.20190302T0300.m06.g031.L1B.std.v02_05.G.190302103326.nc';
cgran = 'SNDR.SNPP.CRIS.20190302T0254.m06.g030.L1B.std.v02_05.G.190302102953.nc';
d4 = read_cris_h5(fullfile(cpath, cgran));

plot(d4.lon(:),d4.lat(:),'b.',d3.Longitude(:),d3.Latitude(:),'r.'); 
  hl = legend('CRIS','AIRS');

bt3 = real(rad2bt(d3.nominal_freq,reshape(d3.radiances,2645,90*135)));

x4 = [d4.wnum_lw; d4.wnum_mw; d4.wnum_sw];
r4 = [d4.rad_lw; d4.rad_mw; d4.rad_sw];
bt4 = real(rad2bt(x4,reshape(r4,2223,9*30*45)));

[eof_maps1,pc1,expvar1,eigv1] = eof_sergiocode2(double(bt1'));  %% AIRS2CHIRP
[eof_maps2,pc2,expvar2,eigv2] = eof_sergiocode2(double(bt2'));  %% CRIS2CHIRP
[eof_maps3,pc3,expvar3,eigv3] = eof_sergiocode2(double(bt3'));  %% AIRS
[eof_maps4,pc4,expvar4,eigv4] = eof_sergiocode2(double(bt4'));  %% CRIS

figure(4); 
  semilogy(1:length(expvar3),expvar3,'b',1:length(expvar4),expvar4,'g',...
           1:length(expvar1),expvar1,'r',1:length(expvar2),expvar2,'m',...
           'linewidth',2)
  hl = legend('AIRS','CRIS','AIRS2CHIRP','CRIS2CHIRP');
  xlim([0 40])
grid

[eof_maps1,pc1,expvar1,eigv1] = eof_sergiocode2(double(bt1')+10*randn(size(bt1')));  %% AIRS2CHIRP
[eof_maps2,pc2,expvar2,eigv2] = eof_sergiocode2(double(bt2')+20*randn(size(bt2')));  %% CRIS2CHIRP
[eof_maps3,pc3,expvar3,eigv3] = eof_sergiocode2(double(bt3')+30*randn(size(bt3')));  %% AIRS
[eof_maps4,pc4,expvar4,eigv4] = eof_sergiocode2(double(bt4')+40*randn(size(bt4')));  %% CRIS

figure(5); 
  semilogy(1:length(expvar3),expvar3,'b',1:length(expvar4),expvar4,'g',...
           1:length(expvar1),expvar1,'r',1:length(expvar2),expvar2,'m',...
           'linewidth',2)
  hl = legend('AIRS + 30Noise','CRIS + 40noise','AIRS2CHIRP + 10noise','CRIS2CHIRP +20noise');
  xlim([0 40])
grid

figure(6); 
  boo = find(d1.wnum >= 1231,1); scatter_coast(d1.lon,d1.lat,10,d1.rad(boo,:)); title('AIRS 2 CHIRP')
figure(7); 
  boo = find(d2.wnum >= 1231,1); scatter_coast(d2.lon,d2.lat,10,d2.rad(boo,:)); title('CRIS 2 CHIRP')
figure(8)
  d3.lon = d3.Longitude;
  d3.lat = d3.Latitude;
  d3.wnum = d3.nominal_freq;
  boo = find(d3.wnum >= 1231,1);  wonk = d3.radiances(boo,:); wonk = wonk(:);
  scatter_coast(d3.lon(:),d3.lat(:),10,wonk); title('AIRS')
figure(9)
  boo = find(x4 >= 1231,1);  wonk = reshape(r4,2223,9*30*45); wonk = wonk(boo,:);
  scatter_coast(d4.lon(:),d4.lat(:),10,wonk); title('CRIS')

figure(6); colormap jet
figure(7); colormap jet
figure(8); colormap jet
figure(9); colormap jet

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AIRS L1c NEdN
nchan_airs = 2645;
nobs_airs = 90 * 135;
nedn_airs = d3.NeN;

% take the mean of valid NEdN values over the full granule
nOK = zeros(nchan_airs, 1);
sOK = zeros(nchan_airs, 1);
for j = 1 : nobs_airs
  iOK = nedn_airs(:, j) < 2;  % flag per-obs valid NEdN values
  nOK = nOK + iOK;
  sOK = sOK + iOK .* nedn_airs(:, j);
end
jOK = nOK > 0;         % flag valid AIRS NEdN values
ntmp1 = sOK ./ nOK;    % mean of all AIRS NEdN values

figure(3); clf
subplot(2,1,1)
x1 = d3.nominal_freq;
y1 = ntmp1;
[x1, y1] = pen_lift(x1, y1);
semilogy(x1, y1)
ylim([0.001, 1.0])
title('AIRS L1c NEdN')
xlabel('wavenumber (cm-1)')
ylabel('mw sr-1 m-2')
grid on

subplot(2,1,2)
x1 = [d4.wnum_lw; d4.wnum_mw; d4.wnum_sw];
y1 = [d4.nedn_lw; d4.nedn_mw; d4.nedn_sw];
[x1, y1] = pen_lift(x1, y1);
semilogy(x1, y1)
ylim([0.001, 1.0])
title('CrIS NEdN')
ylabel('mw sr-1 m-2')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
