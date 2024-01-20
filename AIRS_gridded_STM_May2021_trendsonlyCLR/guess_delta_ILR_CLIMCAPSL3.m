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

if ~exist('climcapsL3_night_file')
  climcapsL3_night_file     = '/asl/s1/sergio/CLIMCAPS_L3/climcapsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat'
  climcapsL3_day_file       = '/asl/s1/sergio/CLIMCAPS_L3/climcapsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_asc.mat';
  climcapsL3_night_file = '/asl/s1/sergio/CLIMCAPS_L3/climcapsclimcaps_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat'
  climcapsL3_day_file   = '/asl/s1/sergio/CLIMCAPS_L3/climcapsclimcaps_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_asc.mat';
end
%junk = load(climcapsL3_day_file,'waterrate','ptemprate','stemprate','trend_nwp_frac');
%junk = load(climcapsL3_night_file,'waterrate','ptemprate','stemprate','trend_nwp_frac');

%% see eg ../FIND_NWP_MODEL_TRENDS/prep_colWV_T_WV_trends_Day_vs_Night.m
if ~exist('resultsWV_C3')
  resultsWV_C3 = fCLIMCAPSL3_night.gas_1trend;
  resultsT_C3  = fCLIMCAPSL3_night.ptemptrend;
  results_C3   = fCLIMCAPSL3_night.stemptrend;
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
  dWVsurf0_CLIMCAPSL3(iii) = resultsWV_C3(p.nlevs(iii)-1,iii);  
  T2m0_CLIMCAPSL3(iii) = p.ptemp(p.nlevs(iii)-1,iii);

  junkind = 1:p.nlevs(iii)-1;
    junk1 = p.ptemp(:,iii);;
    junk = interp1(log(p.plays(junkind,iii)),junk1(junkind),log(p.spres(iii)),meth,'extrap');
    T2mX_CLIMCAPSL3(iii) = junk;

    junk1 = resultsWV_C3(:,iii);
    junk = interp1(log(p.plays(junkind,iii)),log(junk1(junkind)),log(p.spres(iii)),meth,'extrap');
    dWVsurfX_CLIMCAPSL3(iii) = exp(real(junk));

    junk1 = resultsT_C3(:,iii);
    junk = interp1(log(p.plays(junkind,iii)),junk1(junkind),log(p.spres(iii)),meth,'extrap');
    dT2mX_CLIMCAPSL3(iii) = junk;
end

%dWVsurf_CLIMCAPSL3 = dWVsurf0_CLIMCAPSL3;
dWVsurf_CLIMCAPSL3 = dWVsurfX_CLIMCAPSL3;

%dTx2m_CLIMCAPSL3 = results(:,6)';
dTx2m_CLIMCAPSL3 = dT2mX_CLIMCAPSL3;

%Tx2m_CLIMCAPSL3 = p.stemp;
Tx2m_CLIMCAPSL3 = T2mX_CLIMCAPSL3;

eair0_CLIMCAPSL3 = RHSurf0/100 .* esat2(Tx2m_CLIMCAPSL3)/100;                  %% convert RH from percent to fraction, and esat from Pa to mb
eair_CLIMCAPSL3 = eair0_CLIMCAPSL3;

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

dTx2m_CLIMCAPSL3_orig = dTx2m_CLIMCAPSL3;
dWVsurf_CLIMCAPSL3_orig = dWVsurf_CLIMCAPSL3;

if isfield(fCLIMCAPSL3_night,'frac_e2m')
  load ../FIND_NWP_MODEL_TRENDS/desc_2m.mat
  deltaILR_CLIMCAPSL3_from_actual_monthly_trends = compute_ILR_rate_Bruntsaert(eair_CLIMCAPSL3,Tx2m_CLIMCAPSL3,fCLIMCAPSL3_night.frac_e2m,fCLIMCAPSL3_night.t2m);
  plot(rlat,nanmean(reshape(deltaILR_CLIMCAPSL3_from_actual_monthly_trends,72,64),1),'b',ceres_ilr.ceres_trend.trend_lat,ceres_ilr.ceres_trend.trend_lw_clr,'r')
end

[deltaILR_CLIMCAPSL3,ILR_CLIMCAPSL3] = compute_ILR_rate_Bruntsaert(eair_CLIMCAPSL3,Tx2m_CLIMCAPSL3,dWVsurf_CLIMCAPSL3,dTx2m_CLIMCAPSL3);

aslmap(iFig,rlat65,rlon73,smoothn(reshape(deltaILR_CLIMCAPSL3,72,64)',1), [-90 +90],[-180 +180]); title('CLIMCAPSL3 \deltaILR (W/m2) from Brutsaert'); caxis([-1 +1]*0.25); colormap(llsmap5)
  
clf; plot(rlat,nanmean(reshape(deltaILR_CLIMCAPSL3,72,64),1),'b',ceres_ilr.ceres_trend.trend_lat,ceres_ilr.ceres_trend.trend_lw_clr,'r','linewidth',2);
  title('\deltaILR (W/m2) from Brutsaert'); plotaxis2;
  hl = legend('CLIMCAPSL3 Brutsaert','CERES','location','best','fontsize',10);

if isfield(fCLIMCAPSL3_night,'frac_e2m')
  clf; plot(rlat,nanmean(reshape(deltaILR_CLIMCAPSL3,72,64),1),'b',rlat,nanmean(reshape(deltaILR_CLIMCAPSL3_from_actual_monthly_trends,72,64),1),'co-',...
            ceres_ilr.ceres_trend.trend_lat,ceres_ilr.ceres_trend.trend_lw_clr,'r','linewidth',2);
    title('\deltaILR (W/m2) from Brutsaert'); plotaxis2;
    hl = legend('CLIMCAPSL3 Brutsaert','CLIMCAPSL3 Brutsaert from actual monthly trends','CERES','location','best','fontsize',10);
end
  
iY = -1;
if iY == +1 & isfield(fCLIMCAPSL3,'t2m')
   plot(rlat,nanmean(reshape(dTx2m_CLIMCAPSL3_orig,72,64),1),'m',rlat,nanmean(reshape(fCLIMCAPSL3_night.t2m,72,64),1),'r',rlat,nanmean(reshape(dTx2m_orig,72,64),1),'k',...
        rlat,10*nanmean(reshape(dWVsurf_CLIMCAPSL3_orig,72,64),1),'cx-',rlat,10*nanmean(reshape(fCLIMCAPSL3_night.frac_e2m,72,64),1),'bx-',rlat,10*nanmean(reshape(dWVsurf_orig,72,64),1),'gx-','linewidth',2);
  ylim([0 0.075])
  plotaxis2;
  hl = legend('CLIMCAPSL3 interped T2m rate','CLIMCAPSL3 actual monthly T2m rate','CHIRP\_A interped T2m rate',...
              'CLIMCAPSL3 interped fracWV 2m rate  x10','CLIMCAPSL3 actual monthly fracWV 2m rate x10','CHIRP\_A interped fracWV 2m rate x10','location','best','fontsize',10);
end

 
  
