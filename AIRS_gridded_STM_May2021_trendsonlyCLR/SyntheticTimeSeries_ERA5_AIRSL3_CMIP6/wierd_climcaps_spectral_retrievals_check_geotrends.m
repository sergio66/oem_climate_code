addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS

[h,ha,p,pa] = rtpread('../summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_PERTv1.rtp');
z11 = p.stemp;
z12 = z11 * 2;
zz = z11 * 3;
z21 = z11 * 4;
z22 = z11 * 5;
aslmap_2x1x2tiledlayout(z11,z12,zz,z21,z22,1);

[h,ha,p,pa] = rtpread('../summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_PERTv1.rtp');
z11 = p.ptemp;
z11 = reshape(p.ptemp,101,72,64);
z11 = squeeze(nanmean(z11,2));
z12 = z11 * 2;
zz = z11 * 3;
z21 = z11 * 4;
z22 = z11 * 5;
y = p.plevs(:,2300);
x = p.rlat; x = reshape(x,72,64); x = nanmean(x,1);
profile_plots_2x1x2tiledlayout_wide(x,y,z11,z12,zz,z21,z22,2);
profile_plots_2x1x2tiledlayout_tall(x,y,z11,z12,zz,z21,z22,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 64
  iaBad(ii,1:72) = 0;
  fin = ['SimulateTimeSeries/CLIMCAPSL3/reconstruct_climcapsL3_spectra_geo_rlat' num2str(ii,'%02d') '_2002_09_2022_08.mat'];
  a = load(fin);
  plays                = a.zonalplays;
  rlat(ii)             = a.rlatx;
  stemprate(ii,:)      = a.thesave.xst_trend;
  ptempnwprate(:,:,ii) = a.thesave.t2d_xtrendnwp;
  wvnwprate(:,:,ii)    = a.thesave.wv2d_xtrendnwp;
  ptemprate(:,:,ii)    = a.thesave.t2d_xtrend(1:100,:);
  wvrate(:,:,ii)       = a.thesave.wv2d_xtrend(1:100,:);
  if isfield(a.thesave,'bad_lonbins')
    iaBad(ii,a.thesave.bad_lonbins) = 1;
  end
end

figure(1); clf; pcolor(rlat,plays,squeeze(nanmean(wvrate,2)));    shading interp; set(gca,'ydir','reverse'); colorbar; caxis([-0.015 0.015]); title('WV rate'); ylim([100 1000]); set(gca,'yscale','linear')
figure(2); clf; pcolor(rlat,plays,squeeze(nanmean(ptemprate,2))); shading interp; set(gca,'ydir','reverse'); colorbar; caxis([-0.15 0.15]);   title('TZ rate'); ylim([10 1000]);  set(gca,'yscale','log')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iLonBin = 1;
for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'+');
  else
    fprintf(1,'.');
  end

  iaBadCLIMCAPSAIRS(ii,1:72) = 0;
  fin = ['SimulateTimeSeries/CLIMCAPSL3/simulate64binsCLIM_' num2str(ii) '_2002_09_2022_08.op.rtp'];
  [h,ha,p,pa] = rtpread(fin);
  wah = (iLonBin-1) + (1 : 72 : 17280); %% this is time series for iLonBin of 72 latbins
  allTz_CLIMCAPS(ii,:) = p.ptemp(:,1);  %% so just check lonbin 1, timestep 1 of 240
  allWV_CLIMCAPS(ii,:) = p.gas_1(:,1);
  allO3_CLIMCAPS(ii,:) = p.gas_3(:,1);
  allST_CLIMCAPS(ii)   = p.stemp(1);
  if isfield(p,'bad_lonbins')
    iaBadCLIMCAPSAIRS(ii,p.bad_lonbins) = 1;
  end

  iaBadAIRS(ii,1:72) = 0;
  fin = ['SimulateTimeSeries/AIRSL3/simulate64binsAIRSL3_' num2str(ii) '_2002_09_2022_08.op.rtp'];
  [h,ha,p,pa] = rtpread(fin);
  wah = (iLonBin-1) + (1 : 72 : 17280); %% this is time series for iLonBin of 72 latbins
  allTz_AIRS(ii,:) = p.ptemp(:,1);      %% so just check lonbin 1, timestep 1 of 240
  allWV_AIRS(ii,:) = p.gas_1(:,1);
  allO3_AIRS(ii,:) = p.gas_3(:,1);
  allST_AIRS(ii)   = p.stemp(1);
  if isfield(p,'bad_lonbins')
    iaBadAIRS(ii,p.bad_lonbins) = 1;
  end

  fin = ['SimulateTimeSeries/ERA5/simulate64binsERA5_' num2str(ii) '_2002_09_2022_08.op.rtp'];
  [h,ha,p,pa] = rtpread(fin);
  wah = (iLonBin-1) + (1 : 72 : 17280); %% this is time series for iLonBin of 72 latbins
  allTz_ERA5(ii,:) = p.ptemp(:,1);      %% so just check lonbin 1, timestep 1 of 240
  allWV_ERA5(ii,:) = p.gas_1(:,1);
  allO3_ERA5(ii,:) = p.gas_3(:,1);
  allST_ERA5(ii)   = p.stemp(1);
end
fprintf(1,'\n');

figure(3); clf; pcolor(rlat,plays,allTz_CLIMCAPS(:,1:100)'-allTz_AIRS(:,1:100)');        shading interp; colorbar; caxis([-5 +5]);    colormap(usa2); set(gca,'ydir','reverse'); title('T(z) : CLIMCAPS-AIRS');
figure(4); clf; pcolor(rlat,plays,1 - (allWV_CLIMCAPS(:,1:100)'./allWV_AIRS(:,1:100)')); shading interp; colorbar; caxis([-5 +5]/10); colormap(usa2); set(gca,'ydir','reverse'); title('WV(z) : 1-(CLIMCAPS/AIRS)');
