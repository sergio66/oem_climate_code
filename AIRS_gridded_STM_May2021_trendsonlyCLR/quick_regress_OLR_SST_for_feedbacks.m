%% Direct observation of Earth’s spectral long-wave feedback parameter
%% Florian E. Roemer   1,2 , Stefan A. Buehler   1, Manfred Brath   1, Lukas Kluft3 & Viju O. John   4
%% Nature Geoscience 2023, https://doi.org/10.1038/s41561-023-01175-6

roemer_2023_band = [100 570 770 990 1080 1250 2000 2760];
  roemer_farIR_wv  = [100 570];
  roemer_co2       = [570 770];
  roemer_win1      = [770 990];
  roemer_o3        = [990 1080];
  roemer_win2      = [1080 1250];
  roemer_midIR_wv  = [1250 2000]; 
rrtm_bands       = [10 250 500 630 700 820 980 1080 1180 1390 1480 1800 2080 2250 2600 3000];
  rrtm_farIR_wv  = [  1 2];
  rrtm_co2       = [  3 4 5];
  rrtm_win       = [  6 9 ];
  rrtm_o3        = [  7 8];
  rrtm_midIR_wv  = [10 11 12]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indSST    = results(:,6)';
boo0 = -(umbc_spectral_olr.allperts_ecRad.clr - umbc_spectral_olr.olr0_ecRad.clr);
scatter_coast(p.rlon,p.rlat,100,boo0./indSST); colormap jet
junk0U = polyfit(indSST,boo0,1);
[nn,nx,ny,nmean,nstd] = myhist2d(indSST,boo0,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
xlabel('dSST'); ylabel('d(OLR)'); title(['UMBC all and GHG \newline d(OLR) = ' num2str(junk0U(1)) ' d(SST) + ' num2str(junk0U(2))])

boo1 = -(umbc_spectral_olr.allperts_no_tracegas_ecRad.clr - umbc_spectral_olr.olr0_ecRad.clr);
scatter_coast(p.rlon,p.rlat,100,boo1./indSST); colormap jet
junk1U = polyfit(indSST,boo1,1);
[nn,nx,ny,nmean,nstd] = myhist2d(indSST,boo1,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
xlabel('dSST'); ylabel('d(OLR)'); title(['UMBC all but GHG \newline d(OLR) = ' num2str(junk1U(1)) ' d(SST) + ' num2str(junk1U(2))])

%%%%%%%%%%%%%%%%%%%%%%%%%

umbc_spectral_olr.olr0_ecRad.bands.fir_wv = sum(umbc_spectral_olr.olr0_ecRad.bands.clr(rrtm_farIR_wv,:),1);
umbc_spectral_olr.olr0_ecRad.bands.co2    = sum(umbc_spectral_olr.olr0_ecRad.bands.clr(rrtm_co2,:),1);
umbc_spectral_olr.olr0_ecRad.bands.win    = sum(umbc_spectral_olr.olr0_ecRad.bands.clr(rrtm_win,:),1);
umbc_spectral_olr.olr0_ecRad.bands.o3     = sum(umbc_spectral_olr.olr0_ecRad.bands.clr(rrtm_o3,:),1);
umbc_spectral_olr.olr0_ecRad.bands.mid_wv = sum(umbc_spectral_olr.olr0_ecRad.bands.clr(rrtm_midIR_wv,:),1);

umbc_spectral_olr.allperts_ecRad.bands.fir_wv = sum(umbc_spectral_olr.allperts_ecRad.bands.clr(rrtm_farIR_wv,:),1);
umbc_spectral_olr.allperts_ecRad.bands.co2    = sum(umbc_spectral_olr.allperts_ecRad.bands.clr(rrtm_co2,:),1);
umbc_spectral_olr.allperts_ecRad.bands.win    = sum(umbc_spectral_olr.allperts_ecRad.bands.clr(rrtm_win,:),1);
umbc_spectral_olr.allperts_ecRad.bands.o3     = sum(umbc_spectral_olr.allperts_ecRad.bands.clr(rrtm_o3,:),1);
umbc_spectral_olr.allperts_ecRad.bands.mid_wv = sum(umbc_spectral_olr.allperts_ecRad.bands.clr(rrtm_midIR_wv,:),1);

booX = -(umbc_spectral_olr.allperts_ecRad.bands.fir_wv - umbc_spectral_olr.olr0_ecRad.bands.fir_wv);
scatter_coast(p.rlon,p.rlat,100,booX./indSST); colormap jet
junk_fir_wvU = polyfit(indSST,booX,1);
[nn,nx,ny,nmean,nstd] = myhist2d(indSST,booX,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
xlabel('dSST'); ylabel('d(OLR)'); title(['UMBC all and GHG \newline d(OLR) = ' num2str(junk_fir_wvU(1)) ' d(SST) + ' num2str(junk_fir_wvU(2))])

booX = -(umbc_spectral_olr.allperts_ecRad.bands.co2 - umbc_spectral_olr.olr0_ecRad.bands.co2);
scatter_coast(p.rlon,p.rlat,100,booX./indSST); colormap jet
junk_co2U = polyfit(indSST,booX,1);
[nn,nx,ny,nmean,nstd] = myhist2d(indSST,booX,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
xlabel('dSST'); ylabel('d(OLR)'); title(['UMBC all and GHG \newline d(OLR) = ' num2str(junk_co2U(1)) ' d(SST) + ' num2str(junk_co2U(2))])

booX = -(umbc_spectral_olr.allperts_ecRad.bands.win - umbc_spectral_olr.olr0_ecRad.bands.win);
scatter_coast(p.rlon,p.rlat,100,booX./indSST); colormap jet
junk_winU = polyfit(indSST,booX,1);
[nn,nx,ny,nmean,nstd] = myhist2d(indSST,booX,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
xlabel('dSST'); ylabel('d(OLR)'); title(['UMBC all and GHG \newline d(OLR) = ' num2str(junk_winU(1)) ' d(SST) + ' num2str(junk_winU(2))])

booX = -(umbc_spectral_olr.allperts_ecRad.bands.o3 - umbc_spectral_olr.olr0_ecRad.bands.o3);
scatter_coast(p.rlon,p.rlat,100,booX./indSST); colormap jet
junk_o3U = polyfit(indSST,booX,1);
[nn,nx,ny,nmean,nstd] = myhist2d(indSST,booX,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
xlabel('dSST'); ylabel('d(OLR)'); title(['UMBC all and GHG \newline d(OLR) = ' num2str(junk_o3U(1)) ' d(SST) + ' num2str(junk_o3U(2))])

booX = -(umbc_spectral_olr.allperts_ecRad.bands.mid_wv - umbc_spectral_olr.olr0_ecRad.bands.mid_wv);
scatter_coast(p.rlon,p.rlat,100,booX./indSST); colormap jet
junk_mid_wvU = polyfit(indSST,booX,1);
[nn,nx,ny,nmean,nstd] = myhist2d(indSST,booX,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
xlabel('dSST'); ylabel('d(OLR)'); title(['UMBC all and GHG \newline d(OLR) = ' num2str(junk_mid_wvU(1)) ' d(SST) + ' num2str(junk_mid_wvU(2))])

disp(' W/m2/K |  firWV  CO2   WIN   O3   midWV   |  SUM')
junk = [junk_fir_wvU(1) junk_co2U(1) junk_winU(1) junk_o3U(1) junk_mid_wvU(1)];
junk = [junk sum(junk)];
disp('----------------------------------------')
fprintf(1,'UMBC      %5.2f %5.2f %5.2f %5.2f %5.2f | %5.2f \n',junk)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indSST    = era5.trend_stemp;
boo0 = -(era5_spectral_olr.allperts_ecRad.clr - era5_spectral_olr.olr0_ecRad.clr);
scatter_coast(p.rlon,p.rlat,100,boo0./indSST); colormap jet
junk0E = polyfit(indSST,boo0,1);
[nn,nx,ny,nmean,nstd] = myhist2d(indSST,boo0,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
xlabel('dSST'); ylabel('d(OLR)'); title(['ERA5 all and GHG \newline d(OLR) = ' num2str(junk0E(1)) ' d(SST) + ' num2str(junk0E(2))])

boo1 = -(era5_spectral_olr.allperts_no_tracegas_ecRad.clr - era5_spectral_olr.olr0_ecRad.clr);
scatter_coast(p.rlon,p.rlat,100,boo1./indSST); colormap jet
junk1E = polyfit(indSST,boo1,1);
[nn,nx,ny,nmean,nstd] = myhist2d(indSST,boo1,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
xlabel('dSST'); ylabel('d(OLR)'); title(['ERA5 all but GHG \newline d(OLR) = ' num2str(junk1E(1)) ' d(SST) + ' num2str(junk1E(2))])

%%%%%%%%%%%%%%%%%%%%%%%%%

era5_spectral_olr.olr0_ecRad.bands.fir_wv = sum(era5_spectral_olr.olr0_ecRad.bands.clr(rrtm_farIR_wv,:),1);
era5_spectral_olr.olr0_ecRad.bands.co2    = sum(era5_spectral_olr.olr0_ecRad.bands.clr(rrtm_co2,:),1);
era5_spectral_olr.olr0_ecRad.bands.win    = sum(era5_spectral_olr.olr0_ecRad.bands.clr(rrtm_win,:),1);
era5_spectral_olr.olr0_ecRad.bands.o3     = sum(era5_spectral_olr.olr0_ecRad.bands.clr(rrtm_o3,:),1);
era5_spectral_olr.olr0_ecRad.bands.mid_wv = sum(era5_spectral_olr.olr0_ecRad.bands.clr(rrtm_midIR_wv,:),1);

era5_spectral_olr.allperts_ecRad.bands.fir_wv = sum(era5_spectral_olr.allperts_ecRad.bands.clr(rrtm_farIR_wv,:),1);
era5_spectral_olr.allperts_ecRad.bands.co2    = sum(era5_spectral_olr.allperts_ecRad.bands.clr(rrtm_co2,:),1);
era5_spectral_olr.allperts_ecRad.bands.win    = sum(era5_spectral_olr.allperts_ecRad.bands.clr(rrtm_win,:),1);
era5_spectral_olr.allperts_ecRad.bands.o3     = sum(era5_spectral_olr.allperts_ecRad.bands.clr(rrtm_o3,:),1);
era5_spectral_olr.allperts_ecRad.bands.mid_wv = sum(era5_spectral_olr.allperts_ecRad.bands.clr(rrtm_midIR_wv,:),1);

booX = -(era5_spectral_olr.allperts_ecRad.bands.fir_wv - era5_spectral_olr.olr0_ecRad.bands.fir_wv);
scatter_coast(p.rlon,p.rlat,100,booX./indSST); colormap jet
junk_fir_wvE = polyfit(indSST,booX,1);
[nn,nx,ny,nmean,nstd] = myhist2d(indSST,booX,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
xlabel('dSST'); ylabel('d(OLR)'); title(['ERA5 all and GHG \newline d(OLR) = ' num2str(junk_fir_wvE(1)) ' d(SST) + ' num2str(junk_fir_wvE(2))])

booX = -(era5_spectral_olr.allperts_ecRad.bands.co2 - era5_spectral_olr.olr0_ecRad.bands.co2);
scatter_coast(p.rlon,p.rlat,100,booX./indSST); colormap jet
junk_co2E = polyfit(indSST,booX,1);
[nn,nx,ny,nmean,nstd] = myhist2d(indSST,booX,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
xlabel('dSST'); ylabel('d(OLR)'); title(['ERA5 all and GHG \newline d(OLR) = ' num2str(junk_co2E(1)) ' d(SST) + ' num2str(junk_co2E(2))])

booX = -(era5_spectral_olr.allperts_ecRad.bands.win - era5_spectral_olr.olr0_ecRad.bands.win);
scatter_coast(p.rlon,p.rlat,100,booX./indSST); colormap jet
junk_winE = polyfit(indSST,booX,1);
[nn,nx,ny,nmean,nstd] = myhist2d(indSST,booX,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
xlabel('dSST'); ylabel('d(OLR)'); title(['ERA5 all and GHG \newline d(OLR) = ' num2str(junk_winE(1)) ' d(SST) + ' num2str(junk_winE(2))])

booX = -(era5_spectral_olr.allperts_ecRad.bands.o3 - era5_spectral_olr.olr0_ecRad.bands.o3);
scatter_coast(p.rlon,p.rlat,100,booX./indSST); colormap jet
junk_o3E = polyfit(indSST,booX,1);
[nn,nx,ny,nmean,nstd] = myhist2d(indSST,booX,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
xlabel('dSST'); ylabel('d(OLR)'); title(['ERA5 all and GHG \newline d(OLR) = ' num2str(junk_o3E(1)) ' d(SST) + ' num2str(junk_o3E(2))])

booX = -(era5_spectral_olr.allperts_ecRad.bands.mid_wv - era5_spectral_olr.olr0_ecRad.bands.mid_wv);
scatter_coast(p.rlon,p.rlat,100,booX./indSST); colormap jet
junk_mid_wvE = polyfit(indSST,booX,1);
[nn,nx,ny,nmean,nstd] = myhist2d(indSST,booX,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
xlabel('dSST'); ylabel('d(OLR)'); title(['ERA5 all and GHG \newline d(OLR) = ' num2str(junk_mid_wvE(1)) ' d(SST) + ' num2str(junk_mid_wvE(2))])

%disp(' W/m2/K |  firWV  CO2   WIN   O3   midWV   |  SUM')
junk = [junk_fir_wvE(1) junk_co2E(1) junk_winE(1) junk_o3E(1) junk_mid_wvE(1)];
junk = [junk sum(junk)];
disp('----------------------------------------')
fprintf(1,'ERA5      %5.2f %5.2f %5.2f %5.2f %5.2f | %5.2f \n',junk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'regression of OLR vs SST : UMBC allpert (+GHG) and allpert (noGHG) = %8.6f %8.6f \n',junk0U(1),junk1U(1));
fprintf(1,'regression of OLR vs SST : ERA5 allpert (+GHG) and allpert (noGHG) = %8.6f %8.6f \n',junk0E(1),junk1E(1));
