function y = dmmw_dsst_VS_lat(lfmaskA,lfmaskL,lfmaskO,mmw0,dmmw,dsst)

if length(lfmaskA) ~= 72*64  |length(lfmaskO) ~= 72*64 | length(lfmaskL) ~= 72*64 | length(dsst) ~= 72*64 | length(dmmw) ~= 72*64 | length(mmw0) ~= 72*64
  whos lfmask* mmw0 dmmw dsst
  error('need all lfmaskA,lfmaskL,lfmasko,dmmw,dsst,mmw0 to be of length 4608 (72x74)')
end

clear zonal_dmmw_dST zonal_dST
zonal_dmmw_dST = dmmw./mmw0 * 100;
zonal_dmmw_dST = nanmean(reshape(zonal_dmmw_dST,72,64),1);
zonal_dST      = dsst;
zonal_dST      = nanmean(reshape(zonal_dST,72,64),1);
zonal_dmmw_dST = zonal_dmmw_dST ./ zonal_dST;
y.all          = zonal_dmmw_dST;

clear zonal_dmmw_dST zonal_dST
zonal_dmmw_dST = dmmw./mmw0 * 100;
  zonal_dmmw_dST(lfmaskL > 0) = NaN;
zonal_dmmw_dST = nanmean(reshape(zonal_dmmw_dST,72,64),1);
zonal_dST      = dsst;
  zonal_dST(lfmaskL > 0) = NaN;
zonal_dST      = nanmean(reshape(zonal_dST,72,64),1);
zonal_dmmw_dST = zonal_dmmw_dST ./ zonal_dST;
y.ocean        = zonal_dmmw_dST;  %% note we used lfmaskL to get y.ocean

clear zonal_dmmw_dST zonal_dST
zonal_dmmw_dST = dmmw./mmw0 * 100;
  zonal_dmmw_dST(lfmaskO > 0) = NaN;
zonal_dmmw_dST = nanmean(reshape(zonal_dmmw_dST,72,64),1);
zonal_dST      = dsst;
  zonal_dST(lfmaskO > 0) = NaN;
zonal_dST      = nanmean(reshape(zonal_dST,72,64),1);
zonal_dmmw_dST = zonal_dmmw_dST ./ zonal_dST;
y.land         = zonal_dmmw_dST;  %% note we used lfmaskO to get y.land

%%%%%%%%%%%%%%%%%%%%%%%%%

plot(1:64,smooth(y.all,10),'r',1:64,smooth(y.ocean,10),'b',1:64,smooth(y.land,10),'g','linewidth',2); 
plotaxis2; hl = legend('all','ocean','land','location','best','fontsize',10);
ylim([-20 +20])
ylabel('d%mmw/dST'); xlabel('Latitude Bin');
