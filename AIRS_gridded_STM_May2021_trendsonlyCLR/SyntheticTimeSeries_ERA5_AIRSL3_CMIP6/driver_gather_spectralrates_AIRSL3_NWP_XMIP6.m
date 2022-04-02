addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

airsL3_geo_rates         = getdata_AIRSL3vsCLIMCAPSL3(1);
airsCLIMCAPSL3_geo_rates = getdata_AIRSL3vsCLIMCAPSL3(-1);
merra2_geo_rates         = getdata_NWP(2);
era5_geo_rates           = getdata_NWP(5);
cmip6_geo_rates          = getdata_XMIP6(-1);
amip6_geo_rates          = getdata_XMIP6(+1);

airsL3_geo_rates         = airsL3_geo_rates.thestats64x72.stemprate;         
  airsL3_geo_rates = airsL3_geo_rates(:)';
airsCLIMCAPSL3_geo_rates = airsCLIMCAPSL3_geo_rates.thestats64x72.stemprate;
  airsCLIMCAPSL3_geo_rates = airsCLIMCAPSL3_geo_rates(:)';
merra2_geo_rates         = merra2_geo_rates.trend_stemp;
era5_geo_rates           = era5_geo_rates.trend_stemp;
amip6_geo_rates          = amip6_geo_rates.trend_stemp;
cmip6_geo_rates          = cmip6_geo_rates.trend_stemp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

airsL3         = zeros(2645,64*72);
airsCLIMCAPSL3 = zeros(2645,64*72);

merra2 = zeros(2645,64*72);
era5   = zeros(2645,64*72);

amip6 = zeros(2645,64*72);
cmip6 = zeros(2645,64*72);

foundfile = zeros(6,64);
for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'+');
  else
    fprintf(1,'.');
  end
  ind = (ii-1)*72 + (1:72);

  fname = ['SimulateTimeSeries/MERRA2/reconstruct_merra2_spectra_geo_rlat' num2str(ii,'%02d') '.mat'];
  if exist(fname)
    foundfile(1,ii) = 1;
    a = load(fname);
    if ~exist('fchanx')
      fchanx = a.fchanx;
    end
    a = a.thesave.xtrendSpectral;
    merra2(:,ind) = a;
  end
  
  fname = ['SimulateTimeSeries/ERA5/reconstruct_era5_spectra_geo_rlat' num2str(ii,'%02d') '.mat'];
  if exist(fname)
    foundfile(2,ii) = 1;
    a = load(fname);
    a = a.thesave.xtrendSpectral;
    era5(:,ind) = a;
  end
  
  fname = ['SimulateTimeSeries/CMIP6/reconstruct_cmip6_spectra_geo_rlat' num2str(ii,'%02d') '.mat'];
  if exist(fname)
    foundfile(3,ii) = 1;
    a = load(fname);
    a = a.thesave.xtrendSpectral;
    cmip6(:,ind) = a;
  end
  
  fname = ['SimulateTimeSeries/AMIP6/reconstruct_amip6_spectra_geo_rlat' num2str(ii,'%02d') '.mat'];
  if exist(fname)
    foundfile(4,ii) = 1;
    a = load(fname);
    a = a.thesave.xtrendSpectral;
    amip6(:,ind) = a;
  end

  fname = ['SimulateTimeSeries/AIRSL3/reconstruct_airsL3_spectra_geo_rlat' num2str(ii,'%02d') '.mat'];
  if exist(fname)
    foundfile(5,ii) = 1;
    a = load(fname);
    a = a.thesave.xtrendSpectral;
    airsL3(:,ind) = a;
  end
  
  fname = ['SimulateTimeSeries/CLIMCAPSL3/reconstruct_climcapsL3_spectra_geo_rlat' num2str(ii,'%02d') '.mat'];
  if exist(fname)
    foundfile(6,ii) = 1;
    a = load(fname);
    a = a.thesave.xtrendSpectral;
    airsclimcapsL3(:,ind) = a;
  end  
end
fprintf(1,' done \n');
figure(1); imagesc(foundfile); colorbar; title('foundfile'); colormap jet
sum(foundfile,2)

x = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithERA5_uncX3.mat','rates');

%disp('merra2 --> nan'); merra2 = merra2 * 0;
%disp('merra2 --> some nan'); ix = 4608-260+1:4608; merra2(:,ix) = merra2(:,ix) * 0;

figure(2); clf;   plot(fchanx,nanmean(era5'),'b',fchanx,nanmean(merra2'),'r',fchanx,nanmean(x.rates'),'k');           plotaxis2;  hl = legend('era5','merra2','AIRS obs','location','best'); ylim([-0.1 +0.1])
figure(3); clf;   plot(fchanx,nanmean(airsL3'),'b',fchanx,nanmean(airsclimcapsL3'),'r',fchanx,nanmean(x.rates'),'k'); plotaxis2;  hl = legend('airsL3','airsclimcapsL3','AIRS obs','location','best'); ylim([-0.1 +0.1])
figure(4); clf;   plot(fchanx,nanmean(amip6'),'b',fchanx,nanmean(cmip6'),'r',fchanx,nanmean(x.rates'),'k');           plotaxis2;  hl = legend('amip6','cmip6','AIRS obs','location','best'); ylim([-0.1 +0.1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load('llsmap5');
if length(llsmap5) == 64
  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
  llsmap5 = llsmap5(2:64,:);
end

plotoptions.cx = [-1 +1]*0.15; plotoptions.maintitle = 'dBT1231/dt'; plotoptions.plotcolors = llsmap5;
plotoptions.str1 = 'ERA5';    plotoptions.str2 = 'MERRA2';    figure(5); clf; aslmap_2tiledlayout(era5(1520,:),merra2(1520,:),5,plotoptions);
plotoptions.str1 = 'AIRS L3'; plotoptions.str2 = 'CLIMCAPS2'; figure(6); clf; aslmap_2tiledlayout(airsL3(1520,:),airsclimcapsL3(1520,:),6,plotoptions);
plotoptions.str1 = 'AMIP6';   plotoptions.str2 = 'CMIP6';     figure(7); clf; aslmap_2tiledlayout(amip6(1520,:),cmip6(1520,:),7,plotoptions);

plotoptions.cx = [-1 +1]*0.15; plotoptions.maintitle = '(top) dBT1231/dt (bot) dST/dt'; plotoptions.plotcolors = llsmap5;
plotoptions.str1 = 'ERA5';    plotoptions.str2 = 'MERRA2';    plotoptions.str3 = 'ERA5';    plotoptions.str4 = 'MERRA2';    
  figure(5); clf; aslmap_4tiledlayout(era5(1520,:),merra2(1520,:),era5_geo_rates,merra2_geo_rates,5,plotoptions);
plotoptions.str1 = 'AIRS L3'; plotoptions.str2 = 'CLIMCAPS'; plotoptions.str3 = 'AIRS L3'; plotoptions.str4 = 'CLIMCAPS'; 
  figure(6); clf; aslmap_4tiledlayout(airsL3(1520,:),airsclimcapsL3(1520,:),airsL3_geo_rates,airsCLIMCAPSL3_geo_rates,6,plotoptions);
plotoptions.str1 = 'AMIP6';   plotoptions.str2 = 'CMIP6';     plotoptions.str3 = 'AMIP6';   plotoptions.str4 = 'CMIP6';     
  figure(7); clf; aslmap_4tiledlayout(amip6(1520,:),cmip6(1520,:),amip6_geo_rates,cmip6_geo_rates,7,plotoptions);

