%% see Eqn 9, 10 of Analysis of global trends of total column water vapour from multiple years of OMI observations
%% Christian Borger, Steffen Beirle, and Thomas Wagner
%% Articles  Volume 22, issue 16
%% ACP, 22, 10603–10621, 2022  https://doi.org/10.5194/acp-22-10603-2022

Lv_Rv = 6139;  %% units of kelvin

dRH_RH_umbc       = umbc_mmwrate  ./ umbc_mmw0     - Lv_Rv * umbc.stemprate  ./ p.stemp ./ p.stemp;
dRH_RH_era5       = olr.trend_mmw ./ umbc_mmw0     - Lv_Rv * olr.trend_stemp ./ p.stemp ./ p.stemp;
dRH_RH_merra2     = merra2_mmwtrend ./ umbc_mmw0   - Lv_Rv * merra2_stemptrend ./ p.stemp ./ p.stemp;
dRH_RH_airsL3     = airsL3_mmwtrend ./ umbc_mmw0   - Lv_Rv * airsL3.stemprate ./ p.stemp ./ p.stemp;
dRH_RH_climcapsL3 = climcaps_mmwtrend ./ umbc_mmw0 - Lv_Rv * climcapsL3.stemprate ./ p.stemp ./ p.stemp;

if ~exist('RH')
  addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/
  [RH,RH1km,colwater] = layeramt2RH(h,p);
  RHsurf = RH1km(1,:);
end

clear colormapOWG2
colormapOWG = [000 051 078; 000 097 150; 044 139 179; 075 169 178; 049 189 151; 108 228 168; 155 247 157; 209 255 158; 256 256 256; ...
               255 249 204; 252 239 169; 249 215 149; 250 198 130; 236 168 108; 217 133 078; 203 115 062; 143 075 079];
colormapOWG = colormapOWG/256;
for junk = 1 : 3
  colormapOWG2(:,junk) = interp1((1:17)/17,colormapOWG(:,junk),(1:121)/121);
end
colormapOWG2 = flipud(colormapOWG2);
colormapOWG2 = colormapOWG2(1:117,:);
colormapOWG2 = colormapOWG2(1:111,:);

addpath /home/sergio/MATLABCODE/COLORMAP/COLORBREWER/cbrewer/cbrewer/
aslmap(80,rlat65,rlon73,smoothn(reshape(nanmean(dRH_RH_umbc .* RHsurf,1),72,64)',1),  [-90 +90],[-180 +180]); colormap(colormapOWG2); title('UMBC delta(RH)'); caxis([-1 +1])
aslmap(80,rlat65,rlon73,smoothn(reshape(nanmean(dRH_RH_era5.* RHsurf,1),72,64)',1),  [-90 +90],[-180 +180]); colormap(colormapOWG2); title('ERA5 delta(RH)'); caxis([-1 +1])
aslmap(80,rlat65,rlon73,smoothn(reshape(nanmean(dRH_RH_merra2.* RHsurf,1),72,64)',1),  [-90 +90],[-180 +180]); colormap(colormapOWG2); title('MERRA2 delta(RH)'); caxis([-1 +1])

plotRHoptions.str11 = 'ERA5';  plotRHoptions.str12 = 'MERRA2';  plotRHoptions.str21 = 'AIRS L3';  plotRHoptions.str22 = 'CLIMCAPS'; plotRHoptions.str31 = 'This Work';  plotRHoptions.str32 = ' '; 

plotRHoptions.cx = [-1 +1]*0.40; plotRHoptions.cmap = colormapOWG2; plotRHoptions.maintitle = 'd colwater/year [mm/year]';
plotRHoptions.cx = [-1 +1]*0.20; plotRHoptions.cmap = colormapOWG2; plotRHoptions.barstr = 'd colwater/year [mm/year]';
aslmap_3x2tiledlayout(olr.trend_mmw, merra2_mmwtrend, airsL3_mmwtrend, climcaps_mmwtrend, umbc_mmwrate, 0*umbc_mmwrate,80,plotRHoptions);

plotRHoptions.cx = [-1 +1]; plotRHoptions.cmap = colormapOWG2; plotRHoptions.maintitle = 'dcolwater/dt --> dRH/dt surface';
plotRHoptions.cx = [-1 +1]; plotRHoptions.cmap = colormapOWG2; plotRHoptions.barstr = 'dcolwater/dt --> dRH/dt surface';
aslmap_3x2tiledlayout(dRH_RH_era5.* RHsurf, dRH_RH_merra2.* RHsurf, dRH_RH_airsL3.* RHsurf, dRH_RH_climcapsL3.* RHsurf, dRH_RH_umbc.*RHsurf, 0*RHsurf,81,plotRHoptions);
