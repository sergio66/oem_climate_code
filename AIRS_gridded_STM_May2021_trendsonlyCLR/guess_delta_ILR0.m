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
if ~exist('resultsWV')
  resultsWV = fUMBC_night.resultsWV;
  resultsT  = fUMBC_night.resultsT;
  results   = fUMBC_night.results;
end
if ~isfield(p,'plays')
  p.plays = plevs2plays(p.plevs);
end
if ~exist('RHSurf0') & exist('RHSurf')
  RHSurf0 = RHSurf;
end

for iii = 1 : length(p.stemp)
  dWVsurf0(iii) = resultsWV(iii,nlays_straight_from_results(iii));  
  T2m0(iii) = p.ptemp(p.nlevs(iii)-1,iii);

  junkind = 1:p.nlevs(iii)-1;
    junk1 = p.ptemp(:,iii);;
    junk = interp1(log(p.plays(junkind,iii)),junk1(junkind),log(p.spres(iii)),[],'extrap');
    T2mX(iii) = junk;

  junkind = 1:nlays_straight_from_results(iii);
    junk1 = resultsWV(iii,:);
    junk = interp1(log(pavg(junkind)),log(junk1(junkind)),log(p.spres(iii)),[],'extrap');
    dWVsurfX(iii) = exp(real(junk));

  junkind = 1:nlays_straight_from_results(iii);
    junk1 = resultsT(iii,:);
    junk = interp1(log(pavg(junkind)),junk1(junkind),log(p.spres(iii)),[],'extrap');
    dT2mX(iii) = junk;
end

%dWVsurf = dWVsurf0;
dWVsurf = dWVsurfX;

%dTx2m = results(:,6)';
dTx2m = dT2mX;

%Tx2m = p.stemp;
Tx2m = T2mX;

eair0 = RHSurf0/100 .* esat2(Tx2m)/100;                  %% convert RH from percent to fraction, and esat from Pa to mb
eair = eair0;

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

dTx2m_orig = dTx2m;
dWVsurf_orig = dWVsurf;

stefan_boltzmann = 5.67e-8; %% W/m2/K4
deltaILR = stefan_boltzmann*(Tx2m .^4) *1.24/7.*(eair./Tx2m).^(1/7).*(dWVsurf - dTx2m./Tx2m);

iTryOther = -1;
if iTryOther == +1
  
  dTx2m = dTx2m * 0.5;
  deltaILR = stefan_boltzmann*(Tx2m .^4) *1.24/7.*(eair./Tx2m).^(1/7).*(dWVsurf - dTx2m./Tx2m);
  
  dWVsurf = dWVsurf * 2;
  deltaILRpertparam = stefan_boltzmann*(Tx2m .^4) *1.24/7.*(eair./Tx2m).^(1/7).*(dWVsurf - dTx2m./Tx2m);
  
  aslmap(iFig,rlat65,rlon73,smoothn(reshape(deltaILR,72,64)',1), [-90 +90],[-180 +180]); title('\deltaILR (W/m2) from Brutsaert'); caxis([-1 +1]*0.25); colormap(llsmap5)
  
  clf; plot(rlat,nanmean(reshape(deltaILR,72,64),1),'b',rlat,nanmean(reshape(deltaILRpertparam,72,64),1),'c',ceres_ilr.ceres_trend.trend_lat,ceres_ilr.ceres_trend.trend_lw_clr,'r','linewidth',2);
    title('\deltaILR (W/m2) from Brutsaert'); plotaxis2;
    hl = legend('Brutsaert','Brutsaert with WVfrac x2','CERES','location','best');
else
  aslmap(iFig,rlat65,rlon73,smoothn(reshape(deltaILR,72,64)',1), [-90 +90],[-180 +180]); title('\deltaILR (W/m2) from Brutsaert'); caxis([-1 +1]*0.25); colormap(llsmap5)
  
  clf; plot(rlat,nanmean(reshape(deltaILR,72,64),1),'b',ceres_ilr.ceres_trend.trend_lat,ceres_ilr.ceres_trend.trend_lw_clr,'r','linewidth',2);
    title('\deltaILR (W/m2) from Brutsaert'); plotaxis2;
    hl = legend('Brutsaert CHIRP\_A','CERES','location','best');
end
  

iY = -1;
if exist('fERA5_night') & iY > 0
  if isfield(fERA5_night,'ilr_adj')
    wahoo = ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608;
    plot(rlat,nanmean(reshape(deltaILR,72,64),1),'k',...
         rlat,nanmean(reshape(wahoo,72,64),1),'gx-',rlat,nanmean(reshape(fERA5_night.ilr_clr,72,64),1),'c',rlat,nanmean(reshape(fERA5_night.ilr_adj,72,64),1),'b',rlat,nanmean(reshape(fERA5_night.ilr_Rld,72,64),1),'r','linewidth',2);
    plotaxis2;  hl = legend('CHIRP\_A','CERES','ERA5 net?','ERA5 adj','Bruntsaert','location','best','fontsize',10); xlabel('Latitude'); ylabel('dILR/dt W/m2/yr'); xlim([-90 +90])
    title('ILR trends')
  end
end
