addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP

load /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iNumYears = 19;  %% 2002-2021
iNumYears = 12;  %% 2002-2014
iNumYears = 07;  %% 2015-2021
iNumYears = -20;  %% 2002-2022, allsky
iNumYears = +20;  %% 2002-2022, clrsky

iNumYears = input('Enter iNumYears : (12 or 19 or 20 or 07) : ');

% see ../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLD/put_together_kcarta_jacs.m
% normer.normCO2 = 2.2/370;   %% ppm/yr out of 370
% normer.normO3  = 0.01;      %% frac/yr
% normer.normN2O = 1.0/300;   %% ppb/yr outof 300
% normer.normCH4 = 5/1860;    %% ppb/yr out of 1800
% normer.normCFC = 1/1300;    %% ppt/yr out of 1300
% normer.normST  = 0.1;       %% K/yr

%% see driver_quickget_mauna_loaa_rate.m for rates
%% see get_co2_n2o_ch4_for_strow_override.m
if iNumYears == 07
  load iType_8_convert_sergio_clearskygrid_obsonly_Q16.mat

  co2x = 2.262445;
  n2ox = 0.925496;
  ch4x = 6.613673;

  co2x = co2x * ones(1,64);
  n2ox = n2ox * ones(1,64);
  ch4x = ch4x * ones(1,64);

  iVersESRL = 0;
  iVersESRL = 4;
  if iVersESRL == 0
    co2x = co2x;
    n2ox = n2ox;
    ch4x = ch4x;

    co2x = co2x * ones(1,64);
    n2ox = n2ox * ones(1,64);
    ch4x = ch4x * ones(1,64);

  elseif iVersESRL == 4
    esrl_trend = load('/home/sergio/MATLABCODE/ESRL_TRACE_GAS/esrl_co2_ch4_trends_vs_lat_2002_2014_2021.mat');
    n2ox = n2ox;
    co2x = interp1(esrl_trend.rlat,esrl_trend.co2trend_oco2_07,rlat);
    ch4x = interp1(esrl_trend.rlat,esrl_trend.ch4trend_oco2_07,rlat);
  end

  [mean(co2x) mean(n2ox) mean(ch4x)]

elseif iNumYears == 12
  load iType_5_convert_sergio_clearskygrid_obsonly_Q16.mat

  iVersESRL = 3
  iVersESRL = 4
  if iVersESRL == 0
    %% this saved as no_tracegas_spectral_rate_12years0.mat
    co2x = 2.2;
    n2ox = 1.0;
    ch4x = 5.0;

    co2x = co2x * ones(1,64);
    n2ox = n2ox * ones(1,64);
    ch4x = ch4x * ones(1,64);
  
  else
    co2x = 2.065319;
    n2ox = 0.863650;
    ch4x = 5.316252;

    if iVersESRL == 1
      %% this saved as no_tracegas_spectral_rate_12years1.mat
      co2x = 2.065319;
      n2ox = 0.863650;
      ch4x = 5.316252;

      co2x = co2x * ones(1,64);
      n2ox = n2ox * ones(1,64);
      ch4x = ch4x * ones(1,64);
    
    elseif iVersESRL == 2
      %% this saved as no_tracegas_spectral_rate_12years2.mat
      co2x = co2x * 0.98;
      n2ox = n2ox * 0.98;
      ch4x = ch4x * 0.98;

      co2x = co2x * ones(1,64);
      n2ox = n2ox * ones(1,64);
      ch4x = ch4x * ones(1,64);
    
    elseif iVersESRL == 3
      %% this saved as no_tracegas_spectral_rate_12years3.mat
      co2x = co2x * 0.95;
      n2ox = n2ox * 0.95;
      ch4x = ch4x * 0.95;

      co2x = co2x * ones(1,64);
      n2ox = n2ox * ones(1,64);
      ch4x = ch4x * ones(1,64);

    elseif iVersESRL == 4
      %% this saved as no_tracegas_spectral_rate_12years4.mat
      co2x = co2x * 0.95;
      n2ox = n2ox * 0.95;
      ch4x = ch4x * 0.95;

      esrl_trend = load('/home/sergio/MATLABCODE/ESRL_TRACE_GAS/esrl_co2_ch4_trends_vs_lat_2002_2014_2021.mat');
      n2ox = n2ox * ones(1,64);
      co2x = interp1(esrl_trend.rlat,esrl_trend.co2trend_12,rlat);
      ch4x = interp1(esrl_trend.rlat,esrl_trend.ch4trend_12,rlat);
    end
  end
  [mean(co2x) mean(n2ox) mean(ch4x)]

elseif iNumYears == 19
  co2x = 2.262445;
  n2ox = 0.925496;
  ch4x = 6.613673;

  iVersESRL = 0
  iVersESRL = 4
  if iVersESRL == 0
    co2x = co2x * ones(1,64);
    n2ox = n2ox * ones(1,64);
    ch4x = ch4x * ones(1,64);
  elseif iVersESRL == 4
    esrl_trend = load('/home/sergio/MATLABCODE/ESRL_TRACE_GAS/esrl_co2_ch4_trends_vs_lat_2002_2014_2021.mat');
    n2ox = n2ox * ones(1,64);
    co2x = interp1(esrl_trend.rlat,esrl_trend.co2trend_19,rlat);
    ch4x = interp1(esrl_trend.rlat,esrl_trend.ch4trend_19,rlat);
  end

  load iType_4_convert_sergio_clearskygrid_obsonly_Q16.mat

elseif iNumYears == -20
  co2x = 2.262445;
  n2ox = 0.925496;
  ch4x = 6.613673;

  iVersESRL = 0
  iVersESRL = 4
  if iVersESRL == 0
    co2x = co2x * ones(1,64);
    n2ox = n2ox * ones(1,64);
    ch4x = ch4x * ones(1,64);
  elseif iVersESRL == 4
    esrl_trend = load('/home/sergio/MATLABCODE/ESRL_TRACE_GAS/esrl_co2_ch4_trends_vs_lat_2002_2014_2021.mat');
    n2ox = n2ox * ones(1,64);
    co2x = interp1(esrl_trend.rlat,esrl_trend.co2trend_19,rlat);
    ch4x = interp1(esrl_trend.rlat,esrl_trend.ch4trend_19,rlat);
  end

  load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q05.mat

elseif iNumYears == +20
  co2x = 2.262445;
  n2ox = 0.925496;
  ch4x = 6.613673;

  iVersESRL = 0
  iVersESRL = 4
  if iVersESRL == 0
    co2x = co2x * ones(1,64);
    n2ox = n2ox * ones(1,64);
    ch4x = ch4x * ones(1,64);
  elseif iVersESRL == 4
    esrl_trend = load('/home/sergio/MATLABCODE/ESRL_TRACE_GAS/esrl_co2_ch4_trends_vs_lat_2002_2014_2021.mat');
    n2ox = n2ox * ones(1,64);
    co2x = interp1(esrl_trend.rlat,esrl_trend.co2trend_19,rlat);
    ch4x = interp1(esrl_trend.rlat,esrl_trend.ch4trend_19,rlat);
  end

  %load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q05.mat
  load iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat
else
  error('unknown iNumYears')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%keyboard_nowindow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jj = 1 : 64
  fprintf(1,'>>>>>>>>>>>>>>>>>>>>>> latbin %2i of 64 \n',jj);
  if iNumYears == 07
    jacname = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin//kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_07yr_' num2str(jj) '.mat'];
  elseif iNumYears == 12
    jacname = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin//kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_12yr_' num2str(jj) '.mat'];
  elseif iNumYears == 19
    jacname = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin//kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_Dec2021_' num2str(jj) '.mat'];
  elseif iNumYears == -20
    jacname = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin//kcarta_cld_subjac_nostruct_LatBin_kCARTA_ERA5_20yr_CLD_Q09_' num2str(jj,'%02d') '.mat'];
  elseif iNumYears == +20
    jacname = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin//kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_20yr_' num2str(jj,'%02d') '.mat'];
  else
    error('unknown iNumYears')
  end
  for ii = 1 : 72
    iibin = (jj-1)*72 + ii;

    fprintf(1,'\n\n    lat/lon = %2i/%2i \n',jj,ii)
    if iNumYears == 12 | iNumYears == 19 | abs(iNumYears) == 20
      topts.iVersQRenorm = 1;
      topts.iLatX = jj;
      if iNumYears == 20
        [m_ts_jac0,nlays,qrenorm,freq2645]  = get_jac_fast(jacname,iibin,ii,jj,2002+iNumYears,+5,topts);
      else
        [m_ts_jac0,nlays,qrenorm,freq2645]  = get_jac_fast(jacname,iibin,ii,jj,2002+iNumYears,+9,topts);
      end
    elseif iNumYears == 07
      [m_ts_jac0,nlays,qrenorm,freq2645]  = get_jac_fast(jacname,iibin,ii,jj,2015);
    else
      error('unknown iNumYears')
    end
    co2 = m_ts_jac0(:,1) * co2x(jj)/2.2;
    n20 = m_ts_jac0(:,2) * n2ox(jj)/1.0;
    ch4 = m_ts_jac0(:,3) * ch4x(jj)/5.0;
    cfc11 = m_ts_jac0(:,4);
    cfc12 = m_ts_jac0(:,5);
    newrad(:,iibin) = squeeze(b_desc(ii,jj,:)) - co2 - ch4 - n20 - cfc11 - cfc12;
    nedt_b(:,iibin) = squeeze(b_err_desc(ii,jj,:));
    nedt_true(:,iibin) = squeeze(airs_noiseTtrue(ii,jj,:));
    %plot(1:2645,newrad(:,iibin),1:2645,squeeze(b_desc(ii,jj,:)),1:2645,co2)
  end
end

oldrad = reshape(permute(b_desc,[3 1 2]),2645,72*64);

comment = 'see driver_remove_CO2_CH4_N20.m';
if iNumYears == 07
  saver = ['save no_tracegas_spectral_rate_07years_ESRLv' num2str(iVersESRL) '.mat newrad nedt_b nedt_true oldrad comment freq2645 co2x n2ox ch4x'];
elseif iNumYears == 12
  saver = ['save no_tracegas_spectral_rate_12years_ESRLv' num2str(iVersESRL) '.mat newrad nedt_b nedt_true oldrad comment freq2645 co2x n2ox ch4x'];
elseif iNumYears == 19
  saver = ['save no_tracegas_spectral_rate_19years_ESRLv' num2str(iVersESRL) '.mat newrad nedt_b nedt_true oldrad comment freq2645 co2x n2ox ch4x'];
elseif iNumYears == -20
  saver = ['save no_tracegas_spectral_rate_cld20years_ESRLv' num2str(iVersESRL) '.mat newrad nedt_b nedt_true oldrad comment freq2645 co2x n2ox ch4x'];
elseif iNumYears == +20
  saver = ['save no_tracegas_spectral_rate_20years_ESRLv' num2str(iVersESRL) '.mat newrad nedt_b nedt_true oldrad comment freq2645 co2x n2ox ch4x'];
else
  error('unknown iNumYears')
end
eval(saver)

figure(1);
pcolor(reshape(oldrad(50,:)-newrad(50,:),72,64)'); colormap(usa2); colorbar; shading interp; caxis([-1 +1]*0.15)

figure(2)
plot(freq2645,nanmean(oldrad,2),freq2645,nanmean(newrad,2))
plotaxis2; xlim([645 1640])

