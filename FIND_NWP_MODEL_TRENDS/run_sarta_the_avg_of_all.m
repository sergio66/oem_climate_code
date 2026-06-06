addpath0

addpath /home/sergio/git/rtpmake/CLUST_RTPMAKE/COMMON_SETTINGS/
set_path_to_execs

if ~exist('dirout0')
  disp('warning .. you seem to be running yhis as stand alone ... setting dirout = MEAN_PROFILES and iNumYears = 23  ... ')
  disp('   change if needed, also change fout_avg_all if needed')  
  disp('ret to continue'); pause
  
  dirout0 = 'MEAN_PROFILES';
  iNumYears = 23;
  fout_avg_pall = [dirout0 '/avg_ERA5_atm_N_cld_data_2002_09_to_2025_08_desc.mat'];
end

if ~exist('avg_pall')
  loader = ['load ' fout_avg_pall];
  eval(loader);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load /home/sergio/git/oem_climate_jacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries//asc_desc_solzen_time_525_64x72.mat

[h0,ha,p0,pa] = rtpread('/home/sergio/git/matlabcode/REGR_PROFILES_SARTA/REGR49_PROFILES_for_kCARTA_breakouts_for_SARTA/ecmwf_co2_400ppm_1100mb.op.rtp');

avg_pall.gas_2  = p0.gas_2(:,83)*ones(1,4608);
avg_pall.gas_4  = p0.gas_4(:,83)*ones(1,4608);
avg_pall.gas_5  = p0.gas_5(:,83)*ones(1,4608);
avg_pall.gas_6  = p0.gas_6(:,83)*ones(1,4608);
avg_pall.gas_9  = p0.gas_9(:,83)*ones(1,4608);
% avg_pall.gas_11 = p0.gas_11(:,83)*ones(1,4608);
avg_pall.gas_12 = p0.gas_12(:,83)*ones(1,4608);

avg_pall.plat = avg_pall.rlat;
avg_pall.plon = avg_pall.rlon;

junk = reshape(thedata.solzen_desc,4608,526); avg_pall.solzen = nanmean(junk,2)';
junk = reshape(thedata.satzen_desc,4608,526); avg_pall.satzen = nanmean(junk,2)';
avg_pall.zobs = 705000 * ones(size(avg_pall.stemp));

load /home/sergio/git/matlabcode/h2645structure.mat

h0.pfields = 1;
h0.nchan = h.nchan;
h0.ichan = h.ichan;
h0.vchan = h.vchan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fop = ['MEAN_PROFILES/summary_' num2str(iNumYears) 'yrs_era5_monthly_cld.op.rtp'];
frp = ['MEAN_PROFILES/summary_' num2str(iNumYears) 'yrs_era5_monthly_cld.rp.rtp'];

%%%%%
%{
%% do_the_avg_of_all.m has already fixed the clouds
%% avg_pall = fix_clouds_as_needed(avg_pall); %% see /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/

bad = find(avg_pall.cprtop > avg_pall.spres | avg_pall.cprbot > avg_pall.spres | avg_pall.cprtop > avg_pall.cprbot | avg_pall.cprtop < 0 | avg_pall.cprbot < 0 | avg_pall.cfrac > 1 | avg_pall.cfrac12 > 1);
  avg_pall.cfrac(bad) = 0;
  avg_pall.cfrac12(bad) = 0;  
  avg_pall.cngwat(bad) = 0;  
bad = find(avg_pall.cprtop2 > avg_pall.spres | avg_pall.cprbot2 > avg_pall.spres | avg_pall.cprtop2 > avg_pall.cprbot2 | avg_pall.cprtop2 < 0 | avg_pall.cprbot2 < 0 | avg_pall.cfrac2 > 1 | avg_pall.cfrac12 > 1);
  avg_pall.cfrac2(bad) = 0;
  avg_pall.cfrac12(bad) = 0;    
  avg_pall.cngwat2(bad) = 0;  
%}
%%%%%

%{
do_the_avg_of_all.m has already fixed the ctypes
avg_pall.ctype(avg_pall.ctype >  100) = 201;
avg_pall.ctype(avg_pall.ctype <= 101) = 101;
avg_pall.ctype2(avg_pall.ctype2 >  100) = 201;
avg_pall.ctype2(avg_pall.ctype2 <= 101) = 101;
%}

avg_pall = fix_clouds_as_needed(avg_pall); %% see /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
rtpwrite(fop,h0,ha,avg_pall,pa);

sartaer = ['!time ' sartaCld ' fin=' fop ' fout=' frp];
eval(sartaer);

[hx,ha,px,pa] = rtpread(frp);
tcld = rad2bt(hx.vchan,px.rcalc);

colorstr = 'jet';
do_XX_YY_from_X_Y

figure(1); clf; scatter_coast(px.rlon,px.rlat,40,px.stemp); title('ERA5 stemp K');  colormap(colorstr);
figure(1); clf; aslmap(1,rlat65,rlon73,smoothn((reshape(px.stemp,72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 SKT');      colormap(colorstr)
figure(2); clf; aslmap(2,rlat65,rlon73,smoothn((reshape(tcld(1520,:),72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 CLD 1231');  colormap(colorstr)
figure(3); clf; aslmap(3,rlat65,rlon73,smoothn((reshape(px.stemp-tcld(1520,:),72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 SKT - CLD 1231');  colormap(colorstr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fop = ['MEAN_PROFILES/summary_' num2str(iNumYears) 'yrs_era5_monthly_clr.op.rtp'];
frp = ['MEAN_PROFILES/summary_' num2str(iNumYears) 'yrs_era5_monthly_clr.rp.rtp'];

avg_pall.cngwat  = 0 * avg_pall.cngwat;
avg_pall.cngwat2 = 0 * avg_pall.cngwat2;
avg_pall.cfrac   = 0 * avg_pall.cfrac;
avg_pall.cfrac2  = 0 * avg_pall.cfrac2;
avg_pall.cfrac12 = 0 * avg_pall.cfrac12;

rtpwrite(fop,h0,ha,avg_pall,pa);

sartaer = ['!time ' sartaCld ' fin=' fop ' fout=' frp];
eval(sartaer);

[hx,ha,px,pa] = rtpread(frp);
tclr = rad2bt(hx.vchan,px.rcalc);

figure(4); clf; aslmap(4,rlat65,rlon73,smoothn((reshape(tclr(1520,:),72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 CLR 1231');  colormap(colorstr)
figure(5); clf; aslmap(5,rlat65,rlon73,smoothn((reshape(px.stemp - tclr(1520,:),72,64)') ,1), [-90 +90],[-180 +180]); title('SKT - ERA5 CLR 1231');  colormap(colorstr)

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6); clf;
subplot(211);
  plot(meanvaluebin(rlat65),nanmean(reshape(px.stemp,72,64),1),'k',meanvaluebin(rlat65),nanmean(reshape(tclr(1520,:),72,64),1),'b',meanvaluebin(rlat65),nanmean(reshape(tcld(1520,:),72,64),1),'r','linewidth',2)
  legend('SKT','BT 1231 clr','BT 1231 cld','location','best'); ylabel('K'); xlabel('latitude'); xlim([-1 +1]*90)
subplot(212);
  plot(meanvaluebin(rlat65),nanmean(reshape(px.stemp-tclr(1520,:),72,64),1),'b',meanvaluebin(rlat65),nanmean(reshape(px.stemp-tcld(1520,:),72,64),1),'r')
  legend('SKT - BT 1231 clr','SKT - BT 1231 cld','location','best'); ylabel('K'); xlabel('latitude'); xlim([-1 +1]*90)
