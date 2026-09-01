addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD
addpath /asl/matlib/maps

for ii = 1 : 64
  fx = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(ii,'%02d') '/sarta_spectral_trends_const_tracegas_latbin' num2str(ii,'%02d') '_2002_09_2022_08.mat'];
  miaos = load(fx);
  ind = (1:72) + (ii-1)*72;
  cleartrend(:,ind) = miaos.thesave.xtrend;
  clearunc(:,ind) = miaos.thesave.xtrendErr;
end

load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q01.mat
moo01 = b_desc;
unc01 = b_err_desc;

load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q02.mat
moo02 = b_desc;
unc02 = b_err_desc;

load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat
moo03 = b_desc;
unc03 = b_err_desc;

load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q04.mat
moo04 = b_desc;
unc04 = b_err_desc;

load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q05.mat
moo05 = b_desc;
unc05 = b_err_desc;

load h2645structure.mat

zoo01 = permute(moo01,[3 1 2]); zoo01 = reshape(zoo01,2645,72*64);
zoo02 = permute(moo02,[3 1 2]); zoo02 = reshape(zoo02,2645,72*64);
zoo03 = permute(moo03,[3 1 2]); zoo03 = reshape(zoo03,2645,72*64);
zoo04 = permute(moo04,[3 1 2]); zoo04 = reshape(zoo04,2645,72*64);
zoo05 = permute(moo05,[3 1 2]); zoo05 = reshape(zoo05,2645,72*64);

uoo01 = permute(unc01,[3 1 2]); uoo01 = reshape(uoo01,2645,72*64);
uoo02 = permute(unc02,[3 1 2]); uoo02 = reshape(uoo02,2645,72*64);
uoo03 = permute(unc03,[3 1 2]); uoo03 = reshape(uoo03,2645,72*64);
uoo04 = permute(unc04,[3 1 2]); uoo04 = reshape(uoo04,2645,72*64);
uoo05 = permute(unc05,[3 1 2]); uoo05 = reshape(uoo05,2645,72*64);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(1:2645,nanmean(zoo01,2))
plot(h.vchan,nanmean(zoo01,2),h.vchan,nanmean(zoo02,2))
plot(h.vchan,nanmean(zoo01,2),h.vchan,nanmean(zoo02,2),h.vchan,nanmean(zoo03,2),h.vchan,nanmean(zoo04,2),h.vchan,nanmean(zoo05,2))
grid;
xlim([640 1640])
hl = legend('Q1','Q2','Q3','Q4','Q5','location','best','fontsize',8);

woof = nanmean(zoo01,2);
subplot(211); plot(h.vchan,woof,'b',h.vchan,nanmean(cleartrend,2),'r'); ylabel('Q01'); grid; xlim([640 1640]); title('Trend (b) cld (r) ERA5')
subplot(212); plot(h.vchan,woof-nanmean(zoo01,2),h.vchan,woof-nanmean(zoo02,2),h.vchan,woof-nanmean(zoo03,2),h.vchan,woof-nanmean(zoo04,2),h.vchan,woof-nanmean(zoo05,2)); ylabel('Q01-Qx')
grid;
xlim([640 1640]); ax = axis; ylim([-1e-3 ax(4)])
hl = legend('Q1','Q2','Q3','Q4','Q5','location','best','fontsize',8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
plot(1:2645,nanmean(uoo01,2))
plot(h.vchan,nanmean(uoo01,2),h.vchan,nanmean(uoo02,2))
plot(h.vchan,nanmean(uoo01,2),h.vchan,nanmean(uoo02,2),h.vchan,nanmean(uoo03,2),h.vchan,nanmean(uoo04,2),h.vchan,nanmean(uoo05,2))
grid;
xlim([640 1640])
hl = legend('Q1','Q2','Q3','Q4','Q5','location','best','fontsize',8);

woof = nanmean(uoo01,2);
subplot(211); plot(h.vchan,woof,'b',h.vchan,nanmean(clearunc,2),'r'); ylabel('Q01'); grid; xlim([640 1640]); title('Trend Unc (b) cld (r) ERA5')
subplot(212); plot(h.vchan,woof-nanmean(uoo01,2),h.vchan,woof-nanmean(uoo02,2),h.vchan,woof-nanmean(uoo03,2),h.vchan,woof-nanmean(uoo04,2),h.vchan,woof-nanmean(uoo05,2)); ylabel('Q01-Qx')
grid;
xlim([640 1640]); ax = axis; ylim([-1e-3 ax(4)])
hl = legend('Q1','Q2','Q3','Q4','Q5','location','best','fontsize',8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i1385 = find(h.vchan >= 1385,1);
i1419 = find(h.vchan >= 1419,1);
figure(3); pcolor((reshape(zoo01(i1419,:),72,64))');      shading interp; colormap(usa2);  caxis([-1 +1]*0.1); colorbar
figure(4); pcolor((reshape(cleartrend(i1419,:),72,64))'); shading interp; colormap(usa2);  caxis([-1 +1]*0.1); colorbar

 load latB64.mat
  rlat65 = latB2; rlon73 = -180 : 5 : +180;
  rlon = -180 : 5 : +180;  rlat = latB2; 
  rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
  rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
aslmap(3,rlat65,rlon73,smoothn((reshape(zoo01(i1419,:),72,64)') ,1), [-90 +90],[-180 +180]);      title('Cld Rates dBT1419/dt');     caxis([-1 +1]*0.1/2); colormap(usa2)
aslmap(4,rlat65,rlon73,smoothn((reshape(cleartrend(i1419,:),72,64)') ,1), [-90 +90],[-180 +180]); title('Clr Rates dBT1419/dt');     caxis([-1 +1]*0.1/2); colormap(usa2)

aslmap(3,rlat65,rlon73,smoothn((reshape(zoo01(i1385,:),72,64)') ,1), [-90 +90],[-180 +180]);      title('Cld Rates dBT1385/dt');     caxis([-1 +1]*0.1/2); colormap(usa2)
aslmap(4,rlat65,rlon73,smoothn((reshape(cleartrend(i1385,:),72,64)') ,1), [-90 +90],[-180 +180]); title('Clr Rates dBT1385/dt');     caxis([-1 +1]*0.1/2); colormap(usa2)
