function y = precipitation_vs_skt_changes(lat,landfrac)

%% Nadir Jeevanjee 2018 :" Ciple CLimate Models   https://arxiv.org/pdf/1802.02695.pdf
%% claims global precipitation chhanges with SKT as d lnP /dSKT = 0.02/K

%% see https://archive.ipcc.ch/publications_and_data/ar4/wg1/en/ch10s10-3-2.html, Figure 10.6

%% 10.3.2 Patterns of Change in the 21st Century
%% 10.3.2.1 Warming
% Figure 10.6. Zonal means over land and ocean separately, for annual
% mean surface warming (a, b) and precipitation (c, d), shown as ratios
% scaled with the global mean warming (a, c) and not scaled (b,
% d). Multi-model mean results are shown for two scenarios, A2 and
% Commitment (see Section 10.7), for the period 2080 to 2099 relative to
% the zonal means for 1980 to 1999. Results for individual models can be
% seen in the Supplementary Material for this chapter.

%% Fig 10.6.(c) look at A2 land (red line)
landL = [-90 -45 -30   -10   0  10 20  30  90];
landP = [  9  -3  -0.5   0   2   0  2  -2  +12]; 
landP = landP/100;

%% Fig 10.6.(c) look at A2 ocean (black dashed)
oceanL = [-90 -60 -30  -10  0   10 30  90];
oceanP = [ 3   6   -2  -1   12  0  -2  12];
oceanP = oceanP/100;

if nargin == 0
  plot(landL,landP,'r',oceanL,oceanP,'k--','linewidth',2); plotaxis2;
  xl = -90 : +30 : +90;
  xticks(xl);
  
  moo = -90 : +2 : 90;
  plot(moo,interp1(landL,landP,moo),'r',moo,interp1(oceanL,oceanP,moo),'k--','linewidth',2); plotaxis2;
  xlabel('Latitude [deg]'); ylabel('d ln P/dSKT [1/K]');
  
  mooL = interp1(landL,landP,moo);
  mooO = interp1(oceanL,oceanP,moo);
  
  plot(landL,landP,'r',oceanL,oceanP,'k--','linewidth',2); plotaxis2;
  xlabel('Latitude [deg]'); ylabel('d ln P/dSKT [1/K]');
  
  meanL = nansum(cos(moo*pi/180).*mooL)/nansum(cos(moo*pi/180));
  meanO = nansum(cos(moo*pi/180).*mooO)/nansum(cos(moo*pi/180));
  fprintf(1,'mean change : L %8.4f    O %8.4f    L+O %8.4f \n',meanL,meanO,0.3*meanL + 0.7*meanO)

  y = [];

elseif nargin == 2

  ind = find(landfrac > 0);
  if length(ind) > 0
    y(ind) = interp1(landL,landP,lat(ind));
  end

  ind = find(landfrac == 0);
  if length(ind) > 0
    y(ind) = interp1(oceanL,oceanP,lat(ind));
  end
  
end
