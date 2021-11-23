% Trends in continental temperature and humidity directly linked to ocean warming
% Michael P. Byrnea,b,1 and Paul A. O’Gorman
% www.pnas.org/cgi/doi/10.1073/pnas.1722312115
% PNAS | May 8, 2018 | vol. 115 | no. 19 | 4863–4868

% h = cp*T + Lv*q + g*zg = J/kg/K K + J/kg*kg/kg + m/s2*phi/g = J/kg + J/kg + phi     but phi = integral g(z) gz = energy/unit mass = J/kg ??????
% dh = cp dT + Lv dq

Lv = 2.26e6;    % J/kg
cp = 1005;      % J/kg/K
mdair = 28.966; % molecular mass of dry air
mass_g = 18;    % molecular mass of water

ppmvLAY_WV_0    = layers2ppmv(h,p,1:length(p.stemp),1);
ppmvLAY_WV_pert = layers2ppmv(h,pert,1:length(p.stemp),1);

ppmvLAY_WV_0_wet    = layers2ppmv(h,p,1:length(p.stemp),1,-1);
ppmvLAY_WV_pert_wet = layers2ppmv(h,pert,1:length(p.stemp),1,-1);

%% see toppmv.m
% PPMV = MRgkg*((MDAIR*1E-3)/MASSF)*1E+6
% PPMV = MRgg*((MDAIR)/MASSF)*1E+6
junk = mdair/mass_g *1e6;
mmr_0    = ppmvLAY_WV_0 / junk;
mmr_pert = ppmvLAY_WV_pert / junk;
mmr_0_wet    = ppmvLAY_WV_0_wet / junk;
mmr_pert_wet = ppmvLAY_WV_pert_wet / junk;

for ii = 1 : length(p.stemp)
  nlays = p.nlevs(ii)-1;
  dq(ii)     = mmr_pert(nlays,ii) - mmr_0(nlays,ii);      %% g/g or kg/kg
  dq_wet(ii) = mmr_pert_wet(nlays,ii) - mmr_0_wet(nlays,ii);      %% g/g or kg/kg
end
%dq = dq * 1000;  %% convert to g/kg

dh = cp * results(:,6) + Lv * dq';  %% J/kg          J/kg/K K/yr + J/kg kg/kg/yr = J/kg/yr 
dh = dh/1000;                       %% kJ/kg

%% multiply x10 to go to /decade, and for Fig 50 multiply by 1000 to gto from g/g to g/kg as in Fig 2(B) of paper
iFig = 50;
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,10*smoothn((reshape(dq,72,64)')*1000,1),[-90 +90],[-180 +180]); colormap(llsmap5); title('d(mmr)/dt g/kg/decade')
  caxis([-0.5 +0.5])
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,10*smoothn((reshape(dh,72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5); title('PNAS 2018 dh/dt kJ/kg/decade')
  caxis([-1.5 +1.5])

ocean = find(p.landfrac == 0);
land  = find(p.landfrac >  0);
%% multiply x10 to go to /decade, and for Fig 52 multiply by 1000 to gto from g/g to g/kg as in Fig 2(B) of paper
dq_zonal_ocean = 10*dq*1000; dq_zonal_ocean(land) = NaN; dq_zonal_ocean = nanmean(reshape(dq_zonal_ocean,72,64),1);
dq_zonal_land = 10*dq*1000;  dq_zonal_land(ocean) = NaN; dq_zonal_land = nanmean(reshape(dq_zonal_land,72,64),1);
dh_zonal_ocean = 10*dh; dh_zonal_ocean(land) = NaN; dh_zonal_ocean = nanmean(reshape(dh_zonal_ocean,72,64),1);
dh_zonal_land = 10*dh;  dh_zonal_land(ocean) = NaN; dh_zonal_land = nanmean(reshape(dh_zonal_land,72,64),1);
iFig = iFig + 1; figure(iFig); clf;  plot(rlat,dq_zonal_ocean,'bo-',rlat,dq_zonal_land,'rx-','linewidth',2); xlabel('Latitude'); ylabel('Sp. Humidity d(mmr)/dt \newline g/kg/decade');
  xlim([min(rlat) max(rlat)]); plotaxis2; hl = legend('ocean','land','location','best','fontsize',10); 
iFig = iFig + 1; figure(iFig); clf;  plot(rlat,dh_zonal_ocean,'bo-',rlat,dh_zonal_land,'rx-','linewidth',2); xlabel('Latitude'); ylabel('Most Static Energy d(h)/dt \newline kJ/kg/decade');
  xlim([min(rlat) max(rlat)]); plotaxis2; hl = legend('ocean','land','location','best','fontsize',10); 

%% now look at RHsurf changes
% Understanding Decreases in Land Relative Humidity with Global Warming: Conceptual Model and GCM Simulations
% Michael P. Byrnea,b,1 and Paul A. O’Gorman
% JClim Dec 2016 DOI: 10.1175/JCLI-D-16-0351.1
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,smoothn((reshape((RHSurfpert-RHSurf0),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d/dt (surf RH)');  caxis([-1 +1])

deltaRH_per_K = (RHSurfpert-RHSurf0);
deltaRH_per_K = results(:,6)';
deltaRH_per_K = (RHSurfpert-RHSurf0)./(eps+results(:,6)');
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,smoothn((reshape(deltaRH_per_K,72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d/dt (surf RH)/ dST/dt = percent/K');  caxis([-50 +50])
deltaRH_per_K_zonal_ocean = deltaRH_per_K; deltaRH_per_K_zonal_ocean(land) = NaN; deltaRH_per_K_zonal_ocean = nanmean(reshape(deltaRH_per_K_zonal_ocean,72,64),1);
deltaRH_per_K_zonal_land = deltaRH_per_K;  deltaRH_per_K_zonal_land(ocean) = NaN; deltaRH_per_K_zonal_land = nanmean(reshape(deltaRH_per_K_zonal_land,72,64),1);
iFig = iFig + 1; figure(iFig); clf;  plot(rlat,smooth(deltaRH_per_K_zonal_ocean,10),'bo-',rlat,smooth(deltaRH_per_K_zonal_land,10),'rx-','linewidth',2); xlabel('Latitude'); ylabel('dRH/dST percent/K');
  ylim([-20 +20]); xlim([min(rlat) max(rlat)]); plotaxis2; hl = legend('ocean','land','location','best','fontsize',10); 
