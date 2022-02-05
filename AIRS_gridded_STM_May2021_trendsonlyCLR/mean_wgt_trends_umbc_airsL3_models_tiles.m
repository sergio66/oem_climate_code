%% see plot_profile_trends2.m

figure(40); clf
ta = tiledlayout(1,3);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

%%%%%%%%%%%%%%%%%%%%%%%%%

tafov(1) = nexttile;
semilogy(mncos100.*nanmean(coswgt100.*deltaT(1:100,xmask),2),plays,'linewidth',2);
hold on

Tlevs = airsL3.Tlevs;
  boo = zeros(72,64,24); for ijunk = 1 : 24; boo(:,:,ijunk) = xmaskLFmatr'; end
  junk = airsL3.thestats64x72.ptemprate.*boo; junk = reshape(junk,72*64,24);  
  semilogy(mncos024.*nanmean(coswgt024'.*junk(xmask,:),1)',Tlevs,'linewidth',2);

Tlevs = plays;
  boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
  junk = cmip6.trend_ptemp; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  semilogy(mncos100'.*nanmean(coswgt100.*junk(:,xmask),2)',Tlevs,'linewidth',2);   

Tlevs = plays;
  boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
  junk = era5.trend_ptemp; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  semilogy(mncos100'.*nanmean(coswgt100.*junk(:,xmask),2)',Tlevs,'linewidth',2);   

plotaxis2; 
hl = legend(ocbstr,'AIRSL3',mip6str,'ERA5','location','north','fontsize',8);
xlim([-1 +1]*0.05)
hold off; 

%%%%%%%%%%%%%%%%%%%%%%%%%

tafov(2) = nexttile;
semilogy(mncos100.*nanmean(coswgt100.*fracWV(1:100,xmask),2),plays,'linewidth',2); 
hold on

Qlevs = airsL3.Qlevs;
  boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = xmaskLFmatr'; end
  junk = airsL3.thestats64x72.waterrate.*boo; junk = reshape(junk,72*64,12);  
  semilogy(mncos012.*nanmean(coswgt012'.*junk(xmask,:),1)',Qlevs,'linewidth',2); 

Qlevs = plays;
  boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
  junk = cmip6.trend_gas_1; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  semilogy(mncos100'.*nanmean(coswgt100.*junk(:,xmask),2)',Qlevs,'linewidth',2);   

Qlevs = plays;
  boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
  junk = era5.trend_gas_1; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  semilogy(mncos100'.*nanmean(coswgt100.*junk(:,xmask),2)',Qlevs,'linewidth',2);         

plotaxis2; 
%hl = legend(ocbstr,'AIRSL3',mip6str,'ERA5','location','north','fontsize',8);
xlim([-1 +1]*0.01)
hold off; 

%%%%%%%%%%%%%%%%%%%%%%%%%

tafov(3) = nexttile;
semilogy(mncos100.*nanmean(coswgt100.*deltaRH(1:100,xmask),2),plays,'linewidth',2);
hold on

Qlevs = airsL3.Qlevs;
  boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = xmaskLFmatr'; end
  junk = airsL3.thestats64x72.RHrate.*boo; junk = reshape(junk,72*64,12);  
  semilogy(mncos012.*nanmean(coswgt012'.*junk(xmask,:),1)',Qlevs,'linewidth',2); 

Qlevs = plays;
  boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
  junk = cmip6.trend_RH; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  semilogy(mncos100'.*nanmean(coswgt100.*junk(:,xmask),2)',Qlevs,'linewidth',2);   

Qlevs = plays;
  boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
  junk = era5.trend_RH; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  semilogy(mncos100'.*nanmean(coswgt100.*junk(:,xmask),2)',Qlevs,'linewidth',2);   

plotaxis2; 
%hl = legend(ocbstr,'AIRSL3',mip6str,'ERA5','location','north','fontsize',8);
xlim([-1 +1]*0.25)
hold off; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;

set(tafov,'ydir','reverse');   

set(tafov(1),'yscale','log');  
set(tafov(1),'Ylim',[T_Ylim 1000]);

set(tafov(2),'yscale','linear');  
set(tafov(2),'Ylim',[WV_Ylim 1000]);
set(tafov(3),'yscale','linear');  
set(tafov(3),'Ylim',[WV_Ylim 1000]);

set(tafov(2),'yscale','linear');  
set(tafov(2),'Ylim',[WV_Ylim 1000]);
set(tafov(3),'yscale','linear');  
set(tafov(3),'Ylim',[WV_Ylim 1000]);

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';
ta.TileSpacing = 'compact';

% Remove all ytick labels except for 1st column
%for ii = [2 3]
%   tafov(ii).YTickLabel = '';
%end
% Remove all ytick labels from 3rd column
for ii = [3]
   tafov(ii).YTickLabel = '';
end

% Put in xlabel and ylable in the “middle”
tafov(1).XLabel.String = 'dT/dt \newline Kelvin/yr';    tafov(3).XLabel.FontSize = 10;
tafov(2).XLabel.String = 'dWVfrac/dt \newline frac/yr'; tafov(1).XLabel.FontSize = 10;
tafov(3).XLabel.String = 'dRH/dt \newline percent/yr';  tafov(2).XLabel.FontSize = 10;

tafov(1).YLabel.String = 'P(mb)'; tafov(1).YLabel.FontSize = 10;
%tafov(2).YLabel.String = 'P(mb)'; tafov(3).YLabel.FontSize = 10;
%tafov(3).YLabel.String = 'P(mb)'; tafov(3).YLabel.FontSize = 10;

%% put titles
title(tafov(1),'Tz',     'Units', 'normalized', 'Position', [0.5, +1.025, 0]);
title(tafov(2),'WVfrac', 'Units', 'normalized', 'Position', [0.5, +1.025, 0]);
title(tafov(3),'RH',     'Units', 'normalized', 'Position', [0.5, +1.025, 0]);

