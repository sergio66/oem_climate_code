function [co2x,n2ox,ch4x] = get_co2_n2o_ch4_for_strow_override_backwards(driver,settings,iVersJac); %% sets co2x,n2ox,ch4x
  
%{
clear junk
junk.iibin = ((1:64)-1)*72 + 36; %% choose the middle of each latitude set so these are at longitude = 0,latitude = -85 : 3: +85
llll = 0;
for years = -(1:20);
  junk.iNumYears = years;
  llll = llll + 1;
  [co2xx,n2oxx,ch4xx] = get_co2_n2o_ch4_for_strow_override_backwards(junk,[],[]);
  co2x(llll,:) = co2xx;
  n2ox(llll,:) = n2oxx;
  ch4x(llll,:) = ch4xx;
end
figure(1); clf; pcolor(-(1:20) + 2022,1:64,co2x'); colorbar; shading interp; ylabel('Latitude'); colormap(jet); xlabel('Time'); title('CO2 trend between 20XY --> 2022')
figure(2); clf; pcolor(-(1:20) + 2022,1:64,n2ox'); colorbar; shading interp; ylabel('Latitude'); colormap(jet); xlabel('Time'); title('N2O trend between 20XY --> 2022')
figure(3); clf; pcolor(-(1:20) + 2022,1:64,ch4x'); colorbar; shading interp; ylabel('Latitude'); colormap(jet); xlabel('Time'); title('CH4 trend between 20XY --> 2022')
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE_Git/PLOTTER
%lats = -90 : 5 : +90;
%lats = equal_area_spherical_bands(20);
%rlat = 0.5*(lats(1:end-1) + lats(2:end));

ilatjunkxy = floor((driver.iibin-1)/72)+1;
ilonjunkxy = driver.iibin - 72 * (ilatjunkxy - 1);

load /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
  
rlatx = rlat(ilatjunkxy);

co2x = 2.2;
n2ox = 1.0;
ch4x = 5.0;

co2x = 2.2;
n2ox = 0.8;
ch4x = 4.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('get_co2_n2o_ch4_for_strow_override_backwards.m : AIRS 2002/09 - 2022/08')
  co2x = 2.262445;
  n2ox = 0.925496;
  ch4x = 6.613673;

  iVersESRL = 0;
  iVersESRL = 4;
  if iVersESRL == 0
    co2x = co2x;
    n2ox = n2ox;
    ch4x = ch4x;
  elseif iVersESRL == 4
    esrl_trend = load('/home/sergio/MATLABCODE_Git/ESRL_TRACE_GAS/esrl_co2_ch4_trends_vs_lat_2002_Nyears_2022_backwards.mat');
    iNX = 2022 - abs(driver.iNumYears);
    if iNX <= esrl_trend.yyE(end)      
      [~,iNX] = intersect(esrl_trend.yyS,iNX);
    else
      iNX = length(esrl_trend.yyS);
    end
    n2ox = interp1(esrl_trend.rlat0,esrl_trend.n2otrend(iNX,:),rlatx);
    co2x = interp1(esrl_trend.rlat0,esrl_trend.co2trend(iNX,:),rlatx);
    ch4x = interp1(esrl_trend.rlat0,esrl_trend.ch4trend(iNX,:),rlatx);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'iiBin %4i lonbin %2i latbin %2i -->  latitude %8.6f : ESRL trends co2x = %8.6f ppm/yr n2ox = %8.6f ppb/yr ch4x = %8.6f ppb/yr \n',driver.iibin,ilonjunkxy,ilatjunkxy,rlatx,co2x,n2ox,ch4x)

