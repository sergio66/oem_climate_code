%[Tw, Tw1km, Tdew, WBGT, RH, RH1km, colwater, TwSurf, RHSurf, TdewSurf] = layeramt2RH_wet_bulb_dew_point_temperature(h,p);
%[TwX,Tw1kmX,TdewX,WBGTX,RHX,RH1kmX,colwaterX,TwSurfX,RHSurfX,TdewSurfX] = layeramt2RH_wet_bulb_dew_point_temperature(h,pert);

% Eqn 9, see Bk 48 of notes
% Earth Syst. Dynam., 14, 1363–1374, 2023 https://doi.org/10.5194/esd-14-1363-2023
% Understanding variations in downwelling longwave radiation using Brutsaert’s equation
% Yinglin Tian1,2, Deyu Zhong1, Sarosh Alam Ghausi2,3, Guangqian Wang1, and Axel Kleidon2

ceres_ilr = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/ceres_trends_20year_T_ilr.mat');

if ~exist('nlays_straight_from_results')
  load ../AIRS_gridded_STM_May2021_trendsonlyCLR/nlays_straight_from_results_50lays.mat
end

if ~exist('era5_day_file')
  era5_day_file   = '../FIND_NWP_MODEL_TRENDS/ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_asc.mat';
  era5_night_file = '../FIND_NWP_MODEL_TRENDS/ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc.mat';
end
junk = load(era5_day_file,'trend_gas_1','trend_ptemp','trend_stemp','trend_nwp_frac');
junk = load(era5_night_file,'trend_gas_1','trend_ptemp','trend_stemp','trend_nwp_frac');

if ~exist('resultsWV_5')
  resultsWV_5 = junk.trend_gas_1;
  resultsT_5  = junk.trend_ptemp;
  results_5   = junk.trend_stemp;
end
if ~isfield(p,'plays')
  p.plays = plevs2plays(p.plevs);
end
if ~exist('RHSurf0') & exist('RHSurf')
  RHSurf0 = RHSurf;
end

meth = 'spline';
meth = 'makima';
meth = 'linear';

for iii = 1 : length(p.stemp)
  dWVsurf0_ERA5(iii) = resultsWV_5(p.nlevs(iii)-1,iii);  
  T2m0_ERA5(iii) = p.ptemp(p.nlevs(iii)-1,iii);

  junkind = 1:p.nlevs(iii)-1;
    junk1 = p.ptemp(:,iii);;
    junk = interp1(log(p.plays(junkind,iii)),junk1(junkind),log(p.spres(iii)),meth,'extrap');
    T2mX_ERA5(iii) = junk;

    junk1 = resultsWV_5(:,iii);
    junk = interp1(log(p.plays(junkind,iii)),log(junk1(junkind)),log(p.spres(iii)),meth,'extrap');
    dWVsurfX_ERA5(iii) = exp(real(junk));

    junk1 = resultsT_5(:,iii);
    junk = interp1(log(p.plays(junkind,iii)),junk1(junkind),log(p.spres(iii)),meth,'extrap');
    dT2mX_ERA5(iii) = junk;
end

%dWVsurf_ERA5 = dWVsurf0_ERA5;
dWVsurf_ERA5 = dWVsurfX_ERA5;

%dTx2m_ERA5 = results(:,6)';
dTx2m_ERA5 = dT2mX_ERA5;

%Tx2m_ERA5 = p.stemp;
Tx2m_ERA5 = T2mX_ERA5;

eair0_ERA5 = RHSurf0/100 .* esat2(Tx2m_ERA5)/100;                  %% convert RH from percent to fraction, and esat from Pa to mb
eair_ERA5 = eair0_ERA5;

% https://andrewsforest.oregonstate.edu/sites/default/files/lter/data/studies/ms01/dewpt_vpd_calculations.pdf
% From Tetens Formula, the relation between temperature and the partial pressure of water vapor:
%                e(millibars) = 6.1078 exp( (17.269*T) / (237.3+T) )
%           where,
% e is saturated vapor pressure in millibars
% T is temperature in degrees C
%           and the equation for relative humidity:
%                Rh=(ea/es)*100
%           where,
%                 ea is the actual vapor pressure or vapor pressure at dewpoint temperature
%                 es is the saturation vapor pressure or vapor pressure at air temperature
% 
% eair2 = 6.1079 *exp(17.269.*TdewSurf0./(237.3+TdewSurf0));  %% Eqn 7 of paper
% eair2 = 6.1079 *exp(17.269.*(TdewSurf0-273)./(237.3+(TdewSurf0-273)));  %% Eqn 7 of paper
% eair = eair2;

blah.m = 1/7;
blah.A = 0.75;
blah.kw = 0.44;
blah.k2 = blah.kw + 0.065;
blah.k1 = 4*0.0226/293 + blah.k2; 
blah.Rd = 287;
blah.ans = blah.m *  blah.A * (0.622/blah.k2/blah.Rd)^blah.m * beta(blah.k1/blah.k2,blah.m);

dTx2m_ERA5_orig = dTx2m_ERA5;
dWVsurf_ERA5_orig = dWVsurf_ERA5;

onePtwofour = 1.24;
stefan_boltzmann = 5.67e-8; %% W/m2/K4
deltaILR_ERA5 = stefan_boltzmann*(Tx2m_ERA5 .^4) *onePtwofour/7.*(eair_ERA5./Tx2m_ERA5).^(1/7).*(dWVsurf_ERA5 - dTx2m_ERA5./Tx2m_ERA5);

if isfield(fERA5_night,'frac_e2m')
  deltaILR_ERA5_from_actual_monthly_trends = stefan_boltzmann*(Tx2m_ERA5 .^4) *onePtwofour/7.*(eair_ERA5./Tx2m_ERA5).^(1/7).*(fERA5_night.frac_e2m - fERA5_night.t2m./Tx2m_ERA5);  
  %desc_2m_ERA5_240months = load('../FIND_NWP_MODEL_TRENDS/desc_2m.mat');
  plot(rlat,nanmean(reshape(deltaILR_ERA5_from_actual_monthly_trends,72,64),1),'b',ceres_ilr.ceres_trend.trend_lat,ceres_ilr.ceres_trend.trend_lw_clr,'r',rlat,nanmean(reshape(fERA5_night.ilr_Rld,72,64),1),'c')

  load ../FIND_NWP_MODEL_TRENDS/desc_2m.mat
  deltaILR_ERA5_from_actual_monthly_trends = stefan_boltzmann*(nanmean(desc_2m.t2m,1) .^4) *onePtwofour/7.*(nanmean(desc_2m.e2a,1)./nanmean(desc_2m.t2m,1)).^(1/7).*(fERA5_night.frac_e2m - fERA5_night.t2m./nanmean(desc_2m.t2m,1));  

  deltaILR_ERA5_from_actual_monthly_trends = nanmean(desc_2m.Rld,1) .* (fERA5_night.frac_e2m + 27* fERA5_night.t2m./nanmean(desc_2m.t2m,1))/7;  
  plot(rlat,nanmean(reshape(deltaILR_ERA5_from_actual_monthly_trends,72,64),1),'b',ceres_ilr.ceres_trend.trend_lat,ceres_ilr.ceres_trend.trend_lw_clr,'r')
end

iTryOther = -1;
if iTryOther == +1
  
  dTx2m_ERA5 = dTx2m_ERA5 * 0.5;
  deltaILRpertparam_ERA5 = stefan_boltzmann*(Tx2m_ERA5 .^4) *onePtwofour/7.*(eair_ERA5./Tx2m_ERA5).^(1/7).*(dWVsurf_ERA5 - dTx2m_ERA5./Tx2m_ERA5);
  
  dWVsurf_ERA5 = dWVsurf_ERA5 * 2;
  deltaILRpertparam_ERA5 = stefan_boltzmann*(Tx2m_ERA5 .^4) *onePtwofour/7.*(eair_ERA5./Tx2m_ERA5).^(1/7).*(dWVsurf_ERA5 - dTx2m_ERA5./Tx2m_ERA5);
  
  dWVsurf_ERA5 = dWVsurf_ERA5/2;
  Oneptwofour = onePtwofour * 2;
  deltaILRpertparam_ERA5 = stefan_boltzmann*(Tx2m_ERA5 .^4) *onePtwofour/7.*(eair_ERA5./Tx2m_ERA5).^(1/7).*(dWVsurf_ERA5 - dTx2m_ERA5./Tx2m_ERA5);
  
  aslmap(iFig,rlat65,rlon73,smoothn(reshape(deltaILR_ERA5,72,64)',1), [-90 +90],[-180 +180]); title('ERA5 \deltaILR (W/m2) from Brutsaert'); caxis([-1 +1]*0.25); colormap(llsmap5)
  
  clf; plot(rlat,nanmean(reshape(deltaILR_ERA5,72,64),1),'b',rlat,nanmean(reshape(deltaILRpertparam_ERA5,72,64),1),'c',ceres_ilr.ceres_trend.trend_lat,ceres_ilr.ceres_trend.trend_lw_clr,'r','linewidth',2);
    title('\deltaILR (W/m2) from Brutsaert'); plotaxis2;
    hl = legend('ERA5 Brutsaert','ERA5 Brutsaert with WVfrac x2','CERES','location','best');

else
  aslmap(iFig,rlat65,rlon73,smoothn(reshape(deltaILR_ERA5,72,64)',1), [-90 +90],[-180 +180]); title('ERA5 \deltaILR (W/m2) from Brutsaert'); caxis([-1 +1]*0.25); colormap(llsmap5)
  
  clf; plot(rlat,nanmean(reshape(deltaILR_ERA5,72,64),1),'b',ceres_ilr.ceres_trend.trend_lat,ceres_ilr.ceres_trend.trend_lw_clr,'r','linewidth',2);
    title('\deltaILR (W/m2) from Brutsaert'); plotaxis2;
    hl = legend('ERA5 Brutsaert','CERES','location','best');

  if isfield(fERA5_night,'frac_e2m')
    clf; plot(rlat,nanmean(reshape(deltaILR_ERA5,72,64),1),'b',rlat,nanmean(reshape(deltaILR_ERA5_from_actual_monthly_trends,72,64),1),'co-',...
              ceres_ilr.ceres_trend.trend_lat,ceres_ilr.ceres_trend.trend_lw_clr,'r','linewidth',2);
      title('\deltaILR (W/m2) from Brutsaert'); plotaxis2;
      hl = legend('ERA5 Brutsaert','ERA5 Brutsaert from actual monthly trends','CERES','location','best');
  end


end
  
iY = -1;
if iY == +1 & isfield(fERA5_night,'t2m')
   plot(rlat,nanmean(reshape(dTx2m_ERA5_orig,72,64),1),'m',rlat,nanmean(reshape(fERA5_night.t2m,72,64),1),'r',rlat,nanmean(reshape(dTx2m_orig,72,64),1),'k',...
        rlat,10*nanmean(reshape(dWVsurf_ERA5_orig,72,64),1),'cx-',rlat,10*nanmean(reshape(fERA5_night.frac_e2m,72,64),1),'bx-',rlat,10*nanmean(reshape(dWVsurf_orig,72,64),1),'gx-','linewidth',2);
  ylim([0 0.075])
  plotaxis2;
  hl = legend('ERA5 interped T2m rate','ERA5 actual monthly T2m rate','CHIRP\_A interped T2m rate',...
              'ERA5 interped fracWV 2m rate  x10','ERA5 actual monthly fracWV 2m rate x10','CHIRP\_A interped fracWV 2m rate x10','location','best','fontsize',10);
end

 
  