jet64 = jet;
choosecolor = jet64;
choosecolor = llsmap5;
addpath /home/sergio/MATLABCODE//COLORMAP/COLORBREWER/cbrewer/cbrewer/
BrBG = cbrewer('div','BrBG',8);
choosecolor = BrBG;

plotoptions.yReverseDir = -1;
if isfield(plotoptions,'yLimits');
  plotoptions = rmfield(plotoptions,'yLimits');
end
if isfield(plotoptions,'yLinearOrLog');
  plotoptions = rmfield(plotoptions,'yLinearOrLog');
end

plotoptions.str11 = 'ERA5';    plotoptions.str12 = 'MERRA2';    
plotoptions.str21 = 'AIRS L3'; plotoptions.str22 = 'CLIMCAPS L3'; 
plotoptions.str31 = 'CMIP6';   plotoptions.str32 = 'AMIP6';     
plotoptions.str41 = 'UMBC';    plotoptions.str42 = 'MLS L3';     
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';

rlat = load('latB64.mat'); rlat65 = rlat.latB2; rlat = 0.5*(rlat.latB2(1:end-1)+rlat.latB2(2:end));
rlon73 = (1:73); rlon73 = -180 + (rlon73-1)*5;  rlon = (1:72); rlon = -177.5 + (rlon-1)*5;
[Y,X] = meshgrid(rlat,rlon); X = X(:); Y = Y(:);
iFig = 31; figure(iFig); clf; aslmapSergio(rlat65,rlon73,smoothn(reshape(era5_geo_rates(50,:),72,64)',1), [-90 +90],[-180 +180]);

plotoptions.cx = [-1 +1]*0.5; plotoptions.maintitle = 'lay98'; plotoptions.plotcolors = choosecolor;
ix = 98; 
  wx1 = era5_geo_rates(ix,:);   wx2 = merra2_geo_rates(ix,:);
  wx3 = airsL3_geo_rates(ix,:); wx4 = climcaps_geo_rates(ix,:);
  wx5 = cmip6_geo_rates(ix,:);  wx6 = amip6_geo_rates(ix,:);
  wx7 = umbc_geo_rates(ix,:);   wx8 = mls_geo_rates(ix,:);
iFig = 31; figure(iFig); clf; aslmap_8tiledlayout(wx1,wx2,wx3,wx4,wx5,wx6,wx7,wx7*0,iFig,plotoptions);

plotoptions.cx = [-1 +1]*0.5; plotoptions.maintitle = 'lay97'; plotoptions.plotcolors = choosecolor;
ix = 97; 
  wx1 = era5_geo_rates(ix,:);   wx2 = merra2_geo_rates(ix,:);
  wx3 = airsL3_geo_rates(ix,:); wx4 = climcaps_geo_rates(ix,:);
  wx5 = cmip6_geo_rates(ix,:);  wx6 = amip6_geo_rates(ix,:);
  wx7 = umbc_geo_rates(ix,:);   wx8 = mls_geo_rates(ix,:);
iFig = 32; figure(iFig); clf; aslmap_8tiledlayout(wx1,wx2,wx3,wx4,wx5,wx6,wx7,wx7*0,iFig,plotoptions);

plotoptions.cx = [-1 +1]*0.5; plotoptions.maintitle = 'lay96'; plotoptions.plotcolors = choosecolor;
ix = 96; 
  wx1 = era5_geo_rates(ix,:);   wx2 = merra2_geo_rates(ix,:);
  wx3 = airsL3_geo_rates(ix,:); wx4 = climcaps_geo_rates(ix,:);
  wx5 = cmip6_geo_rates(ix,:);  wx6 = amip6_geo_rates(ix,:);
  wx7 = umbc_geo_rates(ix,:);   wx8 = mls_geo_rates(ix,:);
iFig = 33; figure(iFig); clf; aslmap_8tiledlayout(wx1,wx2,wx3,wx4,wx5,wx6,wx7,wx7*0,iFig,plotoptions);

plotoptions.cx = [-1 +1]*0.5; plotoptions.maintitle = 'lay95'; plotoptions.plotcolors = choosecolor;
ix = 95; 
  wx1 = era5_geo_rates(ix,:);   wx2 = merra2_geo_rates(ix,:);
  wx3 = airsL3_geo_rates(ix,:); wx4 = climcaps_geo_rates(ix,:);
  wx5 = cmip6_geo_rates(ix,:);  wx6 = amip6_geo_rates(ix,:);
  wx7 = umbc_geo_rates(ix,:);   wx8 = mls_geo_rates(ix,:);
iFig = 34; figure(iFig); clf; aslmap_8tiledlayout(wx1,wx2,wx3,wx4,wx5,wx6,wx7,wx7*0,iFig,plotoptions);

i850 = find(plays >= 850,1); 
plotoptions.cx = [-1 +1]*0.3; plotoptions.maintitle = '850 mb'; plotoptions.plotcolors = choosecolor;
ix = i850; 
  wx1 = era5_geo_rates(ix,:);   wx2 = merra2_geo_rates(ix,:);
  wx3 = airsL3_geo_rates(ix,:); wx4 = climcaps_geo_rates(ix,:);
  wx5 = cmip6_geo_rates(ix,:);  wx6 = amip6_geo_rates(ix,:);
  wx7 = umbc_geo_rates(ix,:);   wx8 = mls_geo_rates(ix,:);
iFig = 35; figure(iFig); clf; aslmap_8tiledlayout(wx1,wx2,wx3,wx4,wx5,wx6,wx7,wx7*0,iFig,plotoptions);

i500 = find(plays >= 500,1); 
plotoptions.cx = [-1 +1]*0.3; plotoptions.maintitle = '500 mb'; plotoptions.plotcolors = choosecolor;
ix = i500; 
  wx1 = era5_geo_rates(ix,:);   wx2 = merra2_geo_rates(ix,:);
  wx3 = airsL3_geo_rates(ix,:); wx4 = climcaps_geo_rates(ix,:);
  wx5 = cmip6_geo_rates(ix,:);  wx6 = amip6_geo_rates(ix,:);
  wx7 = umbc_geo_rates(ix,:);   wx8 = mls_geo_rates(ix,:);
iFig = 36; figure(iFig); clf; aslmap_8tiledlayout(wx1,wx2,wx3,wx4,wx5,wx6,wx7,wx7*0,iFig,plotoptions);

i250 = find(plays >= 250,1); 
plotoptions.cx = [-1 +1]*0.5; plotoptions.maintitle = '250 mb'; plotoptions.plotcolors = choosecolor;
ix = i250; 
  wx1 = era5_geo_rates(ix,:);   wx2 = merra2_geo_rates(ix,:);
  wx3 = airsL3_geo_rates(ix,:); wx4 = climcaps_geo_rates(ix,:);
  wx5 = cmip6_geo_rates(ix,:);  wx6 = amip6_geo_rates(ix,:);
  wx7 = umbc_geo_rates(ix,:);   wx8 = mls_geo_rates(ix,:);
iFig = 37; figure(iFig); clf; aslmap_8tiledlayout(wx1,wx2,wx3,wx4,wx5,wx6,wx7,wx8,iFig,plotoptions);

i100 = find(plays >= 100,1); 
plotoptions.cx = [-1 +1]*0.5; plotoptions.maintitle = '100 mb'; plotoptions.plotcolors = choosecolor;
ix = i100; 
  wx1 = era5_geo_rates(ix,:);   wx2 = merra2_geo_rates(ix,:);
  wx3 = airsL3_geo_rates(ix,:); wx4 = climcaps_geo_rates(ix,:);
  wx5 = cmip6_geo_rates(ix,:);  wx6 = amip6_geo_rates(ix,:);
  wx7 = umbc_geo_rates(ix,:);   wx8 = mls_geo_rates(ix,:);
iFig = 38; figure(iFig); clf; aslmap_8tiledlayout(wx1,wx2,wx3,wx4,wx5,wx6,wx7,wx8,iFig,plotoptions);

i050 = find(plays >= 050,1); 
plotoptions.cx = [-1 +1]*0.5; plotoptions.maintitle = '050 mb'; plotoptions.plotcolors = choosecolor;
ix = i050; 
  wx1 = era5_geo_rates(ix,:);   wx2 = merra2_geo_rates(ix,:);
  wx3 = airsL3_geo_rates(ix,:); wx4 = climcaps_geo_rates(ix,:);
  wx5 = cmip6_geo_rates(ix,:);  wx6 = amip6_geo_rates(ix,:);
  wx7 = umbc_geo_rates(ix,:);   wx8 = mls_geo_rates(ix,:);
iFig = 39; figure(iFig); clf; aslmap_8tiledlayout(wx1,wx2,wx3,wx4,wx5,wx6,wx7,wx8,iFig,plotoptions);

i010 = find(plays >= 010,1); 
plotoptions.cx = [-1 +1]*0.5; plotoptions.maintitle = '010 mb'; plotoptions.plotcolors = choosecolor;
ix = i010; 
  wx1 = era5_geo_rates(ix,:);   wx2 = merra2_geo_rates(ix,:);
  wx3 = airsL3_geo_rates(ix,:); wx4 = climcaps_geo_rates(ix,:);
  wx5 = cmip6_geo_rates(ix,:);  wx6 = amip6_geo_rates(ix,:);
  wx7 = umbc_geo_rates(ix,:);   wx8 = mls_geo_rates(ix,:);
iFig = 40; figure(iFig); clf; aslmap_8tiledlayout(wx1,wx2,wx3,wx4,wx5,wx6,wx7,wx8,iFig,plotoptions);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Fig 31 : lay 98');
display('Fig 32 : lay 97');
display('Fig 33 : lay 96');
display('Fig 34 : lay 95');
display('Fig 35 : 850 mb');
display('Fig 36 : 500 mb');
display('Fig 37 : 250 mb');
display('Fig 38 : 100 mb');
display('Fig 39 : 050 mb');
display('Fig 40 : 010 mb');
