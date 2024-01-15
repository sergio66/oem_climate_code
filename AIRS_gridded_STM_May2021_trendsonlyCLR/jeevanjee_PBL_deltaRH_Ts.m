function [dRH,lowest_layer_fracwater_dRH_Held_Jeevanjee,lowest_layer_fracwater_dRHzero] = jeevanjee_PLB_deltaRH_Ts(Ts,RHs,dlnPdT,trend_BT1231,p)

%% input Ts = Tsurf, RHs = RHsurface, dlnPdT = Precipitaion Rate
%%       optional BT1231 trend and "p" structure for plots

%% output dRH                                        in percent (so x100 what Nadir Jeevanjee computes, for 1 K temperature change)
%%        lowest_layer_fracwater_dRHzero             what set_CO2_CH4_N2O_ESRL.m computes till Nov 2023
%%        lowest_layer_fracwater_dRH_Held_Jeevanjee  what Isacc Held/Nadir Jeevanjee would compute

%% https://arxiv.org/pdf/1802.02695.pdf
%% The physics of climate change: simple models in climate science
%% Nadir Jeevanjee
%% February 23, 2018
%% Section 6.1

%% Precipitation P (kg/m2/second) depends on Ts as 2-3% change per Kelvin increase in Ts (Sect 5)
%% in steady state, Evaporation = Precipitation       E = P

%% E = Cd  u  rho_v(Ts) (1 - RH(BL))             where Cd ~ 1e-3 dimensionless drag coeff, u = wind speed, rho_v = saturation vapor pressure
%% But E = P and d(ln P)/dTs = 0.02 K-1
%% this d(ln E)/dTs = Cd u    d ln [rho_v(Ts) (1 - RH)]/dTs  = Cd u d ln[rho_v(Ts)]/dTs + d ln(1-RH)/dTs = 0.02 K-1
%%                                                           = Cd u [L/Rv TsTs   -   1/(1-RH) dRH/dTs]   = 0.02 K-1

%% recall 0.02 K-1 = d(lnP)/dTs
%% dRH/dTs = [L/Rv/TsTs - d(lnP)/dTs/(Cd u)](1-RH(BL))

%% But at 300 K, L/Rv/TsTs = 0.06 K-1 and in tropics RH(surface) ~ 0.75 ~ 3/4

%% this dRH/dTs = [0.06 - 0.02] * (1-RH(PBL)) ~ 0.04 K-1 (1/4) ~ 0.01 K-1

lowest_layer_fracwater_dRH_Held_Jeevanjee = [];
lowest_layer_fracwater_dRHzero            = [];

RHs(RHs < 0) = eps;
junk = RHs(:);
if max(junk) > 2
  error('need 0 < RHs < 1');
end

Lv = 2.5e6;     % J/kg
Rv = 431;       % J/kg/K
Lv_Rv = 6139;
Lv_Rv = Lv/Rv;

dlnP_dTs = 0.02;  %% 2018 arXiv Nadir "Simple Climate Models( says this is how ln(Precipitation) changes with Ts
                  %% so does "Constraining global changes in temperature and precipitation from observable changes in surface radiative heating" by Chirag Dhara
                  %% Dhara, C. (2020), Geophysical Research Letters, 47(9), pp. 1–8. doi: 10.1029/2020GL087576.
if nargin >= 3 
  dlnP_dTs = dlnPdT;
end

dRH_dTs = (Lv_Rv./Ts./Ts - dlnP_dTs) .* (1-RHs);     %% this is how RH changes with Ts, accordin to Held (see GFDL Blog #47) and Jeevanjee, arXiv 2021

%% REMEMBER NADIR does everything in terms of fractions, while I do everything in terms of percent (so evetually multiply by 100) ... and here I assume dTs = 1 K
dRH = dRH_dTs * 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin >= 4
  %% see set_CO2_CH4_N2O_ESRL.m how I set FRACTIONAL WATER VAPOR TREND for lowest layer, assuming d(RH) = 0
  
  trend_BT1231 = reshape(trend_BT1231,1,length(trend_BT1231));

  dBT1231 = 0.0025/2;
  dBT1231 = trend_BT1231;
  dBT1231 = abs(trend_BT1231);  %%% new
  dBT1231_WV = dBT1231 * 0.07; %% remember saturation vapor pressure changes at 0.07/K and we want dRH = 0 BUT THIS 0.07 is for 285 K !!!!!!!!!!!!
  
  %%% see my notes, BK 45
  Lo = 2.5e6;  %%% J/kg
  Rv = 461.52; %%% J/kg/K
  moo = exp(Lo/Rv * dBT1231./Ts./Ts)-1;
  moo = Lo/Rv * dBT1231./Ts./Ts;
  dBT1231_WV = moo; %% remember saturation vapor pressure changes at 0.07/K and we want dRH = 0 BUT THIS 0.07 is for 285 K !!!!!!!!!!!!

  %iAdjLowerAtmWVfracX = iAdjLowerAtmWVfrac;                                                    %% orig, too low
  %iAdjLowerAtmWVfracX = iAdjLowerAtmWVfrac * TfacAdjAtmosphericAmplification * iWVFudgeMult;   %% new
  iAdjLowerAtmWVfracX = 1.0;

  %[dBT1231_WV estimate_fracWV_for_deltaRH_zero(dBT1231,Ts) iAdjLowerAtmWVfracX]

  lowest_layer_fracwater_dRHzero = dBT1231_WV * iAdjLowerAtmWVfracX;

  %%%%%%%%%%%%%%%%%%%%%%%%% 
  dRH_Jeevanjee = dRH_dTs .* dBT1231;                                       %% convert Nadir Jeevanjee/Isaac Held RH change per kelvin change in SKT to actual RH change using dBT1231
  lowest_layer_fracwater_dRH_Held_Jeevanjee = dRH_Jeevanjee ./ RHs;         %% convert this to frac change

  %%%%%%%%%%%%%%%%%%%%%%%%%
   
  %% so now compare lowest_layer_fracwater_dRH_Held_Jeevanjee and lowest_layer_fracwater_dRHzero
  if nargin == 5
    figure(1); clf; scatter_coast(p.rlon,p.rlat,100,lowest_layer_fracwater_dRHzero);              title('Sergio : fracWV(PBL) for dRH=0')
      caxis([0 +1]*10e-3); colormap(jet);
    figure(2); clf; scatter_coast(p.rlon,p.rlat,100,dRH);                                         title('Jeevanjee : dRH (in %) for 1 K')
      caxis([-1 +1]*5); colormap(usa2);
    figure(3); clf; scatter_coast(p.rlon,p.rlat,100,lowest_layer_fracwater_dRH_Held_Jeevanjee);   title('Jeevanjee : fracWV(PBL)')
      caxis([0 +1]*10e-3); colormap(jet);
    figure(4); clf; scatter_coast(p.rlon,p.rlat,100,lowest_layer_fracwater_dRHzero + lowest_layer_fracwater_dRH_Held_Jeevanjee); title('Sum')
      caxis([0 +1]*10e-3); colormap(jet);
    figure(5); clf; scatter_coast(p.rlon,p.rlat,100,lowest_layer_fracwater_dRH_Held_Jeevanjee./lowest_layer_fracwater_dRHzero + 1); title('1 + Held/Sergio')
      colormap(jet); caxis([0 2])
    figure(6); clf; scatter_coast(p.rlon,p.rlat,100,dlnPdT); title('d ln(Precip)/dSKT')
      colormap(usa2); caxis([-1 +1]*0.1)
    
  end
end
