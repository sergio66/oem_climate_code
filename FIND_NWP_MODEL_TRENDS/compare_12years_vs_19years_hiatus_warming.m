if ~exist('era5_12')
  merra2_12 = getdata_NWP(2,1,0,12); figure(1); title('MERRA2 : 2002/09-2014/08 dST/dt K/yr'); disp('ret to continue'); pause; 
  merra2_19 = getdata_NWP(2,1,0,19); figure(1); title('MERRA2 : 2002/09-2021/08 dST/dt K/yr'); disp('ret to continue'); pause;

  era5_12 = getdata_NWP(5,1,0,12); figure(1); title('ERA5 : 2002/09-2014/08 dST/dt K/yr'); disp('ret to continue'); pause;
  era5_19 = getdata_NWP(5,1,0,19); figure(1); title('ERA5 : 2002/09-2021/08 dST/dt K/yr'); disp('ret to continue'); pause;

  airsV7_12 = getdata_AIRSL3vsCLIMCAPSL3(1,1,0,12); figure(1); title('AIRSv7 : 2002/09-2014/08 dST/dt K/yr'); disp('ret to continue'); pause;
  airsV7_19 = getdata_AIRSL3vsCLIMCAPSL3(1,1,0,19); figure(1); title('AIRSv7 : 2002/09-2021/08 dST/dt K/yr'); disp('ret to continue'); pause;

  climcaps_12 = getdata_AIRSL3vsCLIMCAPSL3(-1,1,0,12); figure(1); title('CLIMCAPS : 2002/09-2014/08 dST/dt K/yr'); disp('ret to continue'); pause;
  climcaps_19 = getdata_AIRSL3vsCLIMCAPSL3(-1,1,0,19); figure(1); title('CLIMCAPS : 2002/09-2021/08 dST/dt K/yr'); disp('ret to continue'); pause;

  amip6_12 = getdata_XMIP6(+1);    figure(1); title('AMIP6 : 2002/09-2014/08 dST/dt K/yr'); disp('ret to continue'); pause;
  cmip6_12 = getdata_XMIP6(-1);    figure(1); title('CMIP6 : 2002/09-2014/08 dST/dt K/yr'); disp('ret to continue'); pause;
  
  era5_warming.trend_gas_1 = (era5_19.trend_gas_1 * 19 - era5_12.trend_gas_1 * 12)/(19-12);
  era5_warming.trend_gas_3 = (era5_19.trend_gas_3 * 19 - era5_12.trend_gas_3 * 12)/(19-12);
  era5_warming.trend_RH = (era5_19.trend_RH * 19 - era5_12.trend_RH * 12)/(19-12);
  era5_warming.trend_ptemp = (era5_19.trend_ptemp * 19 - era5_12.trend_ptemp * 12)/(19-12);
  era5_warming.trend_stemp = (era5_19.trend_stemp * 19 - era5_12.trend_stemp * 12)/(19-12);

  merra2_warming.trend_gas_1 = (merra2_19.trend_gas_1 * 19 - merra2_12.trend_gas_1 * 12)/(19-12);
  merra2_warming.trend_gas_3 = (merra2_19.trend_gas_3 * 19 - merra2_12.trend_gas_3 * 12)/(19-12);
  merra2_warming.trend_RH = (merra2_19.trend_RH * 19 - merra2_12.trend_RH * 12)/(19-12);
  merra2_warming.trend_ptemp = (merra2_19.trend_ptemp * 19 - merra2_12.trend_ptemp * 12)/(19-12);
  merra2_warming.trend_stemp = (merra2_19.trend_stemp * 19 - merra2_12.trend_stemp * 12)/(19-12);
end

addpath /home/sergio/MATLABCODE/PLOTTER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iCntr = 3000;
iNorD = +1;
iAorOorL = 0;

  maskLF = zeros(1,4608);
  if iAorOorL == 0
    maskLF = ones(1,4608);
  elseif iAorOorL == -1
    maskLF(landfrac == 1) = 1;
  elseif iAorOorL == +1
    maskLF(landfrac == 0) = 1;
  end
  maskLFmatr = reshape(maskLF,72,64)';

if ~exist('rlat')
  load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
  rlat65 = latB2; rlon73 = -180 : 5 : +180;
  rlon = -180 : 5 : +180;  rlat = latB2;
  rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
  rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
end

if ~exist('pavgLAY')
  boo = load('/home/sergio/MATLABCODE/airslevels.dat');
  pjunkN = boo(1:100)-boo(2:101);
  pjunkD = log(boo(1:100)./boo(2:101));
  pavgLAY = pjunkN./pjunkD;
  pavgLAY = flipud(pavgLAY);
  pavgLAY = pavgLAY*ones(1,iCntr);
end

if ~exist('maskLF')
  maskLF = zeros(1,4608);
  if iAorOorL == 0
    maskLF = ones(1,4608);
  elseif iAorOorL == -1
    maskLF(landfrac == 1) = 1;
  elseif iAorOorL == +1
    maskLF(landfrac == 0) = 1;
  end
  maskLFmatr = reshape(maskLF,72,64)';
end

if iNorD > 0
  strNorD = 'NIGHT';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strChoice  = 'ERA5';

iFig = 0;
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,maskLFmatr.*smoothn(reshape(era5_warming.trend_stemp,72,64)',1),[-90 +90],[-180 +180]);   colormap(llsmap5); caxis([-0.15 +0.15]); title([strNorD ' stemp d/dt ' strChoice ' K/yr']); 

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era5_warming.trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([10 1000]); title([strNorD ' dT/dt ' strChoice ' K/yr']);

%%%%%%%%%%%%%%%%%%%%%%%%%

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era5_warming.trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([100 1000]); title([strNorD ' dRH/dt ' strChoice ' percent/yr']);

%%%%%%%%%%%%%%%%%%%%%%%%%

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era5_warming.trend_gas_1,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.01); ylim([100 1000]); title([strNorD ' d(fracWV)/dt ' strChoice ' 1/yr']);

figure(1); title('ERA5 : 2014/09-2021/08 dST/dt K/yr'); disp('ret to continue'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strChoice  = 'MERRA2';

iFig = 0;
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,maskLFmatr.*smoothn(reshape(merra2_warming.trend_stemp,72,64)',1),[-90 +90],[-180 +180]);   colormap(llsmap5); caxis([-0.15 +0.15]); title([strNorD ' stemp d/dt ' strChoice ' K/yr']); 

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(merra2_warming.trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([10 1000]); title([strNorD ' dT/dt ' strChoice ' K/yr']);

%%%%%%%%%%%%%%%%%%%%%%%%%

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(merra2_warming.trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); ylim([100 1000]); title([strNorD ' dRH/dt ' strChoice ' percent/yr']);

%%%%%%%%%%%%%%%%%%%%%%%%%

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(merra2_warming.trend_gas_1,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,iCntr),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.01); ylim([100 1000]); title([strNorD ' d(fracWV)/dt ' strChoice ' 1/yr']);

figure(1); title('MERRA2 : 2014/09-2021/08 dST/dt K/yr'); disp('ret to continue'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


