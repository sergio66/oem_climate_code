olr = load('OLR_ecRad/ERA5/all_era5_olr.mat');

figure(51); clf; aslmap(51,rlat65,rlon73,smoothn(reshape(nanmean(olr.xmmw,1),72,64)',1),  [-90 +90],[-180 +180]); colormap(jet); title('ERA5 mmw');
figure(52); clf; aslmap(52,rlat65,rlon73,smoothn(reshape(nanmean(olr.xstemp,1),72,64)',1),[-90 +90],[-180 +180]); colormap(jet); title('ERA5 stemp');
figure(53); clf; aslmap(53,rlat65,rlon73,smoothn(reshape(nanmean(olr.xolr,1),72,64)',1),  [-90 +90],[-180 +180]); colormap(jet); title('ERA5 olr clr');

figure(54); clf; aslmap(54,rlat65,rlon73,smoothn(reshape(olr.trend_mmw,72,64)',1),  [-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('ERA5 trend mmw');
figure(55); clf; aslmap(55,rlat65,rlon73,smoothn(reshape(olr.trend_stemp,72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('ERA5 trend stemp');

umbc_mmw0 = mmwater_rtp(h,p);
pUMBC = p;
pUMBC.stemp = pUMBC.ptemp + umbc.stemprate;
pUMBC.ptemp(1:100,:) = pUMBC.ptemp(1:100,:) + umbc.ptemprate;
pUMBC.gas_1(1:100,:) = pUMBC.gas_1(1:100,:) .*( 1 + umbc.waterrate);
umbc_mmwrate = mmwater_rtp(h,pUMBC) - umbc_mmw0;
figure(56); clf; aslmap(56,rlat65,rlon73,smoothn(reshape(umbc_mmwrate,72,64)',1),  [-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('UMBC trend mmw');
figure(57); clf; aslmap(57,rlat65,rlon73,smoothn(reshape(umbc.stemprate,72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('UMBC trend stemp');

junk = load('MERRA2_atm_data_2002_09_to_2022_08_trends_desc.mat','trend_mmw'); merra2_mmwtrend = junk.trend_mmw;
figure(58); clf; aslmap(58,rlat65,rlon73,smoothn(reshape(merra2_mmwtrend,72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('MERRA2 trend mmw');
junk = load('MERRA2_atm_data_2002_09_to_2022_08_trends_desc.mat','trend_stemp'); merra2_stemptrend = junk.trend_stemp;
figure(59); clf; aslmap(59,rlat65,rlon73,smoothn(reshape(merra2_stemptrend,72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('MERRA2 trend stemp');

figure(60); clf; aslmap(60,rlat65,rlon73,smoothn(reshape(olr.trend_olr,72,64)',1),   [-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('ERA5 trend olr clr');

figure(61); clf; plot(trend_rlat64,nanmean(reshape(olr.trend_mmw,72,64),1),'b',trend_rlat64,nanmean(reshape(merra2_mmwtrend,72,64),1),'g',trend_rlat64,nanmean(reshape(umbc_mmwrate,72,64),1),'r','linewidth',2);
  plotaxis2; hl = legend('ERA5','MERRA2','UMBC','location','best'); title('mmw trends')
figure(62); clf; plot(trend_rlat64,nanmean(reshape(olr.trend_stemp,72,64),1),'b',trend_rlat64,nanmean(reshape(merra2_stemptrend,72,64),1),'g',trend_rlat64,nanmean(reshape(umbc.stemprate,72,64),1),'r','linewidth',2);
  plotaxis2; hl = legend('ERA5','MERRA2','UMBC','location','best'); title('stemp trends')
figure(63); clf; plot(nanmean(reshape(olr.trend_stemp,72,64),1),nanmean(reshape(olr.trend_mmw,72,64),1),'bs',nanmean(reshape(merra2_stemptrend,72,64),1),nanmean(reshape(merra2_mmwtrend,72,64),1),'gd',...
                      nanmean(reshape(umbc.stemprate,72,64),1),nanmean(reshape(umbc_mmwrate,72,64),1),'rx','markersize',5,'linewidth',2)
  plotaxis2; hl = legend('ERA5','MERRA2','UMBC','location','best'); xlabel('stemp trends'); ylabel('mmw trends')
