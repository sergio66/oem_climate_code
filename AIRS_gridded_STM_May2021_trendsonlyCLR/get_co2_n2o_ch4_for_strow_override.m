function [co2x,n2ox,ch4x] = get_co2_n2o_ch4_for_strow_override(driver,iVersJac); %% sets co2x,n2ox,ch4x

%strlatbin = num2str(floor((driver.iibin-1)/72)+1,'%02d');
rlatxy = floor((driver.iibin-1)/72)+1;

load /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

rlatx = rlat(rlatxy);

co2x = 2.2;
n2ox = 1.0;
ch4x = 5.0;

co2x = 2.2;
n2ox = 0.8;
ch4x = 4.5;

if iVersJac == 2019
  %% CRIS 2012/05 - 2019/04
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
    esrl_trend = load('/home/sergio/MATLABCODE/ESRL_TRACE_GAS/esrl_co2_ch4_trends_vs_lat_2002_2014_2021.mat');
    n2ox = n2ox;
    co2x = interp1(esrl_trend.rlat,esrl_trend.co2trend_cris07,rlatx);
    ch4x = interp1(esrl_trend.rlat,esrl_trend.ch4trend_cris07,rlatx);
  end

elseif iVersJac == 2021
  %% AIRS 2002/09 - 2021/08
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
    esrl_trend = load('/home/sergio/MATLABCODE/ESRL_TRACE_GAS/esrl_co2_ch4_trends_vs_lat_2002_2014_2021.mat');
    n2ox = n2ox;
    co2x = interp1(esrl_trend.rlat,esrl_trend.co2trend_19,rlatx);
    ch4x = interp1(esrl_trend.rlat,esrl_trend.ch4trend_19,rlatx);
  end

elseif iVersJac == 2014
  %% CMIP6 2002/09 - 2014/08
  iVersESRL = 3;
  iVersESRL = 4;
  if iVersESRL == 0
    %% this saved as no_tracegas_spectral_rate_12years0.mat
    co2x = 2.2;
    n2ox = 1.0;
    ch4x = 5.0;

    co2x = co2x;
    n2ox = n2ox;
    ch4x = ch4x;
  
  else
    co2x = 2.065319;
    n2ox = 0.863650;
    ch4x = 5.316252;

    if iVersESRL == 1
      %% this saved as no_tracegas_spectral_rate_12years1.mat
      co2x = 2.065319;
      n2ox = 0.863650;
      ch4x = 5.316252;

      co2x = co2x;
      n2ox = n2ox;
      ch4x = ch4x;
    
    elseif iVersESRL == 2
      %% this saved as no_tracegas_spectral_rate_12years2.mat
      co2x = co2x * 0.98;
      n2ox = n2ox * 0.98;
      ch4x = ch4x * 0.98;

      co2x = co2x;
      n2ox = n2ox;
      ch4x = ch4x;
    
    elseif iVersESRL == 3
      %% this saved as no_tracegas_spectral_rate_12years3.mat
      co2x = co2x * 0.95;
      n2ox = n2ox * 0.95;
      ch4x = ch4x * 0.95;

      co2x = co2x;
      n2ox = n2ox;
      ch4x = ch4x;

    elseif iVersESRL == 4
      %% this saved as no_tracegas_spectral_rate_12years4.mat
      co2x = co2x * 0.95;
      n2ox = n2ox * 0.95;
      ch4x = ch4x * 0.95;

      esrl_trend = load('/home/sergio/MATLABCODE/ESRL_TRACE_GAS/esrl_co2_ch4_trends_vs_lat_2002_2014_2021.mat');
      n2ox = n2ox;
      co2x = interp1(esrl_trend.rlat,esrl_trend.co2trend_12,rlatx);
      ch4x = interp1(esrl_trend.rlat,esrl_trend.ch4trend_12,rlatx);
    end
  end

end

fprintf(1,'iiBin %4i latbin %2i latitude %8.6f : ESRL trends co2x = %8.6f ppm/yr n2ox = %8.6f ppb/yr ch4x = %8.6f ppb/yr \n',driver.iibin,rlatxy,rlatx,co2x,n2ox,ch4x)