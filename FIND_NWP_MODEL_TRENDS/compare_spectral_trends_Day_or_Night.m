disp('plotting spectral chisqr between D and N spectra')

iFrac = input('Show Day (-1) or Night (+1/default) or average (0) or difference (-10): ');
if length(iFrac) == 0
  iFrac = 1;
end
if iFrac == +1
  nightfrac = 1;
  dayfrac = 0;
  ghgfrac = 1.0;
elseif iFrac == -1
  nightfrac = 0;
  dayfrac = 1;
  ghgfrac = 1.0;
elseif iFrac == 0
  nightfrac = 0.5;
  dayfrac = 0.5;
  ghgfrac = 1.0;
elseif iFrac == -10
  nightfrac = -1.0;
  dayfrac = +1.0;
  ghgfrac = 0.0;
end

junk = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/h2645structure.mat');
f = junk.h.vchan;

if ~exist('raaBadFov')
  
  junk = load(umbc_day_file,'rates','nedt');   spectral_rates_day = junk.rates;   nedt_day = junk.nedt;
  junk = load(umbc_night_file,'rates','nedt'); spectral_rates_night = junk.rates; nedt_day = junk.nedt;
  chanset = 1 : 2000;
  %settings.iIgnoreChans_CH4 = -1;
  %plot_spectral_region_chisqr(spectral_rates_day(chanset,:),0*spectral_rates_day(chanset,:),0*spectral_rates_day(chanset,:),spectral_rates_night(chanset,:),f(chanset,:),nedt_day(chanset,:),-1,settings,[]);
  settings.iIgnoreChans_CH4 = -1;
  [raaBadFov,indBadFov,chisqrX,chisqrR] = plot_spectral_region_chisqr(spectral_rates_day(chanset,:),0*spectral_rates_day(chanset,:),0*spectral_rates_day(chanset,:),spectral_rates_night(chanset,:),f(chanset,:),nedt_day(chanset,:));
end
  
ind = 1 : 4608;                           allbias = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
ind = find(p.landfrac == 1);              landbias = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
ind = find(p.landfrac == 0);              oceanbias = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
figure(95); clf; plot(f,allbias,'k',f,landbias,'g',f,oceanbias,'b');
xlim([640 1620]); plotaxis2; hl = legend('All','Land','Ocean','location','best','fontsize',10); title('Day')

ind = find(p.rlat < -60);                 spolarbias_lo = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
ind = find(p.rlat > +60);                 npolarbias_lo = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
ind = find(p.rlat >= -60 & p.rlat < -30); smidlatbias_lo = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
ind = find(p.rlat >= +30 & p.rlat < +60); nmidlatbias_lo = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
ind = find(abs(p.rlat) < 30);             tropicalbias_lo = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
figure(96); clf; plot(f,allbias,'k',f,tropicalbias_lo,'g',f,smidlatbias_lo,'m--',f,nmidlatbias_lo,'r',f,spolarbias_lo,'c--',f,npolarbias_lo,'b')
xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','best','fontsize',10); title('Day OCEAN + LAND')

%figure(97); clf; 
%  plot(f,0.33*landbias+0.67*oceanbias,f,0.33*tropicalbias_l+0.67*tropicalbias_o,'g',f,0.33*smidlatbias_l+0.67*smidlatbias_o,'m--',f,0.33*nmidlatbias_l+0.67*nmidlatbias_o,'r',...
%       f,0.33*spolarbias_l+0.67*spolarbias_o,'c--',f,0.33*npolarbias_l+0.67*npolarbias_o,'b')
%xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','best','fontsize',10); title('Day 0.33LAND+0.67OCEAN')

jett = jet; jett(1,:) = 1;
maskLF = ones(size(chisqrX));
figure(98); clf;
aslmap(98,rlat65,rlon73,smoothn((reshape(maskLF.*chisqrX,72,64)') ,1), [-90 +90],[-180 +180]); title('Chisqr');     caxis([0 +1]*10); colormap(jett)
clear plotoptions;
plotoptions.maintitle = '\chi^2'; plotoptions.cmap = jet(64); plotoptions.cmap = jett;
plotoptions.str11 = 'ALL';
plotoptions.str12 = 'T(z)';
plotoptions.str21 = 'Window';
plotoptions.str22 = '5*WV(z)';
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
plotoptions.yLinearOrLog = +1;

plotoptions.cx = [0 +1]*10; 
z11 = maskLF.*chisqrR.iAll;
z12 = maskLF.*chisqrR.i15um;
z21 = maskLF.*chisqrR.iWindow;
z22 = maskLF.*chisqrR.iWV;

plotoptions.cx = log10(plotoptions.cx);
%plotoptions.cx = [-2 +1];
z11 = log10(z11);  z12 = log10(z12);  z21 = log10(z21);  z22 = log10(z22); 

iFig = 98; figure(iFig); clf;
aslmap_2x2tiledlayout(z11,z12,z21,z22*5,iFig,plotoptions);
disp('so looks like the chisqr is largest over land because B/N spectra differ in (b) T(z) and also(c) window region. Also there are differences in WV band in zones 35-50 N and 30-50 S');

if abs(iFrac) <= 1
  for junk = 95:97; figure(junk); axis([640 1620 -0.1 +0.1]); end
else
  for junk = 95:97; figure(junk); axis([640 1620 -0.02 +0.01]); end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = find(p.landfrac == 0 & p.rlat < -60);                 spolarbias_o = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
ind = find(p.landfrac == 0 & p.rlat > +60);                 npolarbias_o = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
ind = find(p.landfrac == 0 & p.rlat >= -60 & p.rlat < -30); smidlatbias_o = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
ind = find(p.landfrac == 0 & p.rlat >= +30 & p.rlat < +60); nmidlatbias_o = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
ind = find(p.landfrac == 0 & abs(p.rlat) < 30);             tropicalbias_o = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
  signal_o = [tropicalbias_o smidlatbias_o nmidlatbias_o spolarbias_o npolarbias_o]; 

ind = find(p.landfrac > 0 & p.rlat < -60);                  spolarbias_l = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
ind = find(p.landfrac > 0 & p.rlat > +60);                  npolarbias_l = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
ind = find(p.landfrac > 0 & p.rlat >= -60 & p.rlat < -30);  smidlatbias_l = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
ind = find(p.landfrac > 0 & p.rlat >= +30 & p.rlat < +60);  nmidlatbias_l = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
ind = find(p.landfrac > 0 & abs(p.rlat) < 30);              tropicalbias_l = nanmean(dayfrac*spectral_rates_day(:,ind) + nightfrac*spectral_rates_night(:,ind),2);
  signal_l = [tropicalbias_l smidlatbias_l nmidlatbias_l spolarbias_l npolarbias_l];

figure(99); clf; 
  ta = tiledlayout(2,1);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  
  tafov(1) = nexttile; 
  plot(f,oceanbias,'k',f,tropicalbias_o,'g',f,smidlatbias_o,'m--',f,nmidlatbias_o,'r',f,spolarbias_o,'c--',f,npolarbias_o,'b')
  xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','best','fontsize',6); ylabel('OCEAN')
  title('Obervations');

  if abs(iFrac) <= 1
    axis([640 1620 -0.1 +0.1]);
  else
    axis([640 1620 -0.02 +0.01]);
  end

  tafov(2) = nexttile; plot(f,landbias,'k',f,tropicalbias_l,'g',f,smidlatbias_l,'m--',f,nmidlatbias_l,'r',f,spolarbias_l,'c--',f,npolarbias_l,'b')
  xlim([640 1620]); plotaxis2; hl = legend('All','tropical','S Midlat','N Midlat','S Polar','N Polar','location','best','fontsize',6); ylabel('LAND')

  if abs(iFrac) <= 1
    axis([640 1620 -0.1 +0.1]);
  else
    axis([640 1620 -0.02 +0.01]);
  end

  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  % Remove all xtick labels except for 2nd row
  for ii = [1]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
  end
  for ii = [1 2]
   tafov(ii).FontSize = 8;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('now SIMULATING average spectra baased on AVERAGE geophysicsl trends .. should really be reaading in the SAVED spectral trends')

get_avg_profile

simulate_spectral_trends_Day_or_Night_profiles
simulate_spectral_trends_Day_or_Night_profiles_X

