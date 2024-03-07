%% ocean
ppX = pprad0_ocean;

junk = find(p.landfrac == 0 & abs(p.rlat) < 30);
  boo = nanmean(dayfrac*fUMBC_day.results(junk,6)    + nightfrac*fUMBC_night.results(junk,6));        ppX.stemp(1) = ppX.stemp(1) + boo;
  boo = nanmean(dayfrac*fUMBC_day.ptemptrend(:,junk) + nightfrac*fUMBC_night.ptemptrend(:,junk),2);   ppX.ptemp(1:100,1) = ppX.ptemp(1:100,1) + reshape(boo,100,1);
  boo = nanmean(dayfrac*fUMBC_day.gas_1trend(:,junk) + nightfrac*fUMBC_night.gas_1trend(:,junk),2);   ppX.gas_1(1:100,1) = ppX.gas_1(1:100,1) .* (1 + reshape(boo,100,1));
  ppX.gas_2(1:100,1) = ppX.gas_2(1:100,1) .* (1 + 2.2/400 * ghgfrac);
  ppX.gas_4(1:100,1) = ppX.gas_4(1:100,1) .* (1 + 1/300 * ghgfrac);
  ppX.gas_6(1:100,1) = ppX.gas_6(1:100,1) .* (1 + 6/1800 * ghgfrac);
junk = find(p.landfrac == 0 & p.rlat > -60 & p.rlat < -30);
  boo = nanmean(dayfrac*fUMBC_day.results(junk,6)    + nightfrac*fUMBC_night.results(junk,6));        ppX.stemp(2) = ppX.stemp(2) + boo;
  boo = nanmean(dayfrac*fUMBC_day.ptemptrend(:,junk) + nightfrac*fUMBC_night.ptemptrend(:,junk),2);   ppX.ptemp(1:100,2) = ppX.ptemp(1:100,2) + reshape(boo,100,1);
  boo = nanmean(dayfrac*fUMBC_day.gas_1trend(:,junk) + nightfrac*fUMBC_night.gas_1trend(:,junk),2);   ppX.gas_1(1:100,2) = ppX.gas_1(1:100,2) .* (1 + reshape(boo,100,1));
  ppX.gas_2(1:100,2) = ppX.gas_2(1:100,2) .* (1 + 2.2/400 * ghgfrac);
  ppX.gas_4(1:100,2) = ppX.gas_4(1:100,2) .* (1 + 1/300 * ghgfrac);
  ppX.gas_6(1:100,2) = ppX.gas_6(1:100,2) .* (1 + 6/1800 * ghgfrac);
junk = find(p.landfrac == 0 & p.rlat > +30 & p.rlat < +60);
  boo = nanmean(dayfrac*fUMBC_day.results(junk,6)    + nightfrac*fUMBC_night.results(junk,6));        ppX.stemp(3) = ppX.stemp(3) + boo;
  boo = nanmean(dayfrac*fUMBC_day.ptemptrend(:,junk) + nightfrac*fUMBC_night.ptemptrend(:,junk),2);   ppX.ptemp(1:100,3) = ppX.ptemp(1:100,3) + reshape(boo,100,1);
  boo = nanmean(dayfrac*fUMBC_day.gas_1trend(:,junk) + nightfrac*fUMBC_night.gas_1trend(:,junk),2);   ppX.gas_1(1:100,3) = ppX.gas_1(1:100,3) .* (1 + reshape(boo,100,1));
  ppX.gas_2(1:100,3) = ppX.gas_2(1:100,3) .* (1 + 2.2/400 * ghgfrac);
  ppX.gas_4(1:100,3) = ppX.gas_4(1:100,3) .* (1 + 1/300 * ghgfrac);
  ppX.gas_6(1:100,3) = ppX.gas_6(1:100,3) .* (1 + 6/1800 * ghgfrac);
junk = find(p.landfrac == 0 & p.rlat < -60);
  boo = nanmean(dayfrac*fUMBC_day.results(junk,6)    + nightfrac*fUMBC_night.results(junk,6));        ppX.stemp(4) = ppX.stemp(4) + boo;
  boo = nanmean(dayfrac*fUMBC_day.ptemptrend(:,junk) + nightfrac*fUMBC_night.ptemptrend(:,junk),2);   ppX.ptemp(1:100,4) = ppX.ptemp(1:100,4) + reshape(boo,100,1);
  boo = nanmean(dayfrac*fUMBC_day.gas_1trend(:,junk) + nightfrac*fUMBC_night.gas_1trend(:,junk),2);   ppX.gas_1(1:100,4) = ppX.gas_1(1:100,4) .* (1 + reshape(boo,100,1));
  ppX.gas_2(1:100,4) = ppX.gas_2(1:100,4) .* (1 + 2.2/400 * ghgfrac);
  ppX.gas_4(1:100,4) = ppX.gas_4(1:100,4) .* (1 + 1/300 * ghgfrac);
  ppX.gas_6(1:100,4) = ppX.gas_6(1:100,4) .* (1 + 6/1800 * ghgfrac);
junk = find(p.landfrac == 0 & p.rlat > +60);
  boo = nanmean(dayfrac*fUMBC_day.results(junk,6)    + nightfrac*fUMBC_night.results(junk,6));        ppX.stemp(5) = ppX.stemp(5) + boo;
  boo = nanmean(dayfrac*fUMBC_day.ptemptrend(:,junk) + nightfrac*fUMBC_night.ptemptrend(:,junk),2);   ppX.ptemp(1:100,5) = ppX.ptemp(1:100,5) + reshape(boo,100,1);
  boo = nanmean(dayfrac*fUMBC_day.gas_1trend(:,junk) + nightfrac*fUMBC_night.gas_1trend(:,junk),2);   ppX.gas_1(1:100,5) = ppX.gas_1(1:100,5) .* (1 + reshape(boo,100,1));
  ppX.gas_2(1:100,5) = ppX.gas_2(1:100,5) .* (1 + 2.2/400 * ghgfrac);
  ppX.gas_4(1:100,5) = ppX.gas_4(1:100,5) .* (1 + 1/300 * ghgfrac);
  ppX.gas_6(1:100,5) = ppX.gas_6(1:100,5) .* (1 + 6/1800 * ghgfrac);

rtpwrite(fop,hh2,hha2,ppX,ppa2);
eval(sartaer);

[~,~,ppradX,~] = rtpread(frp);
simbiasO = rad2bt(hh2.vchan,ppradX.rcalc)-rad2bt(hh2.vchan,pprad0_ocean.rcalc); 

pprad0_ocean.plays = plevs2plays(pprad0_ocean.plevs);
figure(100); clf; 
  ta = tiledlayout(2,2);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  
  tafov(1) = nexttile; 
  plot(ppradX.ptemp(1:97,1)-pprad0_ocean.ptemp(1:97,1),pprad0_ocean.plays(1:97,1),'g',ppradX.ptemp(1:97,2)-pprad0_ocean.ptemp(1:97,2),pprad0_ocean.plays(1:97,2),'m--',ppradX.ptemp(1:97,3)-pprad0_ocean.ptemp(1:97,3),pprad0_ocean.plays(1:97,3),'r',...
       ppradX.ptemp(1:97,4)-pprad0_ocean.ptemp(1:97,4),pprad0_ocean.plays(1:97,4),'c--',ppradX.ptemp(1:97,5)-pprad0_ocean.ptemp(1:97,5),pprad0_ocean.plays(1:97,5),'b','linewidth',2)
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); plotaxis2; title('T'); axis([-0.1 +0.1 10 1000])
  if iFrac == -10;  axis([-0.1/2 +0.1/2 10 1000]); end
  hl = legend('tropical','S Midlat','N Midlat','S Polar','N Polar','location','best','fontsize',6); 

  tafov(2) = nexttile; 
  plot(ppradX.gas_1(1:97,1)./pprad0_ocean.gas_1(1:97,1),pprad0_ocean.plays(1:97,1),'g',ppradX.gas_1(1:97,2)./pprad0_ocean.gas_1(1:97,2),pprad0_ocean.plays(1:97,2),'m--',ppradX.gas_1(1:97,3)./pprad0_ocean.gas_1(1:97,3),pprad0_ocean.plays(1:97,3),'r',...
       ppradX.gas_1(1:97,4)./pprad0_ocean.gas_1(1:97,4),pprad0_ocean.plays(1:97,4),'c--',ppradX.gas_1(1:97,5)./pprad0_ocean.gas_1(1:97,5),pprad0_ocean.plays(1:97,5),'b','linewidth',2); axis([0.99 1.01 10 1000])
  if iFrac == -10;  axis([1-0.01/2 1+0.01/2 10 1000]); end
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); plotaxis2([1]); title('WV')

disp('OCEAN : dSTEMP and dMMW for TRP, NML,SMP, NP,SP')
fprintf(1,' %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',[ppradX.stemp - pprad0_ocean.stemp; mmwater_rtp(hh2,ppradX)-mmwater_rtp(hh2,pprad0_ocean)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% land
ppX = pprad0_land;

junk = find(p.landfrac == 1 & abs(p.rlat) < 30);
  boo = nanmean(dayfrac*fUMBC_day.results(junk,6)    + nightfrac*fUMBC_night.results(junk,6));        ppX.stemp(1) = ppX.stemp(1) + boo;
  boo = nanmean(dayfrac*fUMBC_day.ptemptrend(:,junk) + nightfrac*fUMBC_night.ptemptrend(:,junk),2);   ppX.ptemp(1:100,1) = ppX.ptemp(1:100,1) + reshape(boo,100,1);
  boo = nanmean(dayfrac*fUMBC_day.gas_1trend(:,junk) + nightfrac*fUMBC_night.gas_1trend(:,junk),2);   ppX.gas_1(1:100,1) = ppX.gas_1(1:100,1) .* (1 + reshape(boo,100,1));
  ppX.gas_2(1:100,1) = ppX.gas_2(1:100,1) .* (1 + 2.2/400 * ghgfrac);
  ppX.gas_4(1:100,1) = ppX.gas_4(1:100,1) .* (1 + 1/300 * ghgfrac);
  ppX.gas_6(1:100,1) = ppX.gas_6(1:100,1) .* (1 + 6/1800 * ghgfrac);
junk = find(p.landfrac == 1 & p.rlat > -60 & p.rlat < -30);
  boo = nanmean(dayfrac*fUMBC_day.results(junk,6)    + nightfrac*fUMBC_night.results(junk,6));        ppX.stemp(2) = ppX.stemp(2) + boo;
  boo = nanmean(dayfrac*fUMBC_day.ptemptrend(:,junk) + nightfrac*fUMBC_night.ptemptrend(:,junk),2);   ppX.ptemp(1:100,2) = ppX.ptemp(1:100,2) + reshape(boo,100,1);
  boo = nanmean(dayfrac*fUMBC_day.gas_1trend(:,junk) + nightfrac*fUMBC_night.gas_1trend(:,junk),2);   ppX.gas_1(1:100,2) = ppX.gas_1(1:100,2) .* (1 + reshape(boo,100,1));
  ppX.gas_2(1:100,2) = ppX.gas_2(1:100,2) .* (1 + 2.2/400 * ghgfrac);
  ppX.gas_4(1:100,2) = ppX.gas_4(1:100,2) .* (1 + 1/300 * ghgfrac);
  ppX.gas_6(1:100,2) = ppX.gas_6(1:100,2) .* (1 + 6/1800 * ghgfrac);
junk = find(p.landfrac == 1 & p.rlat > +30 & p.rlat < +60);
  boo = nanmean(dayfrac*fUMBC_day.results(junk,6)    + nightfrac*fUMBC_night.results(junk,6));        ppX.stemp(3) = ppX.stemp(3) + boo;
  boo = nanmean(dayfrac*fUMBC_day.ptemptrend(:,junk) + nightfrac*fUMBC_night.ptemptrend(:,junk),2);   ppX.ptemp(1:100,3) = ppX.ptemp(1:100,3) + reshape(boo,100,1);
  boo = nanmean(dayfrac*fUMBC_day.gas_1trend(:,junk) + nightfrac*fUMBC_night.gas_1trend(:,junk),2);   ppX.gas_1(1:100,3) = ppX.gas_1(1:100,3) .* (1 + reshape(boo,100,1));
  ppX.gas_2(1:100,3) = ppX.gas_2(1:100,3) .* (1 + 2.2/400 * ghgfrac);
  ppX.gas_4(1:100,3) = ppX.gas_4(1:100,3) .* (1 + 1/300 * ghgfrac);
  ppX.gas_6(1:100,3) = ppX.gas_6(1:100,3) .* (1 + 6/1800 * ghgfrac);
junk = find(p.landfrac == 1 & p.rlat < -60);
  boo = nanmean(dayfrac*fUMBC_day.results(junk,6)    + nightfrac*fUMBC_night.results(junk,6));        ppX.stemp(4) = ppX.stemp(4) + boo;
  boo = nanmean(dayfrac*fUMBC_day.ptemptrend(:,junk) + nightfrac*fUMBC_night.ptemptrend(:,junk),2);   ppX.ptemp(1:100,4) = ppX.ptemp(1:100,4) + reshape(boo,100,1);
  boo = nanmean(dayfrac*fUMBC_day.gas_1trend(:,junk) + nightfrac*fUMBC_night.gas_1trend(:,junk),2);   ppX.gas_1(1:100,4) = ppX.gas_1(1:100,4) .* (1 + reshape(boo,100,1));
  ppX.gas_2(1:100,4) = ppX.gas_2(1:100,4) .* (1 + 2.2/400 * ghgfrac);
  ppX.gas_4(1:100,4) = ppX.gas_4(1:100,4) .* (1 + 1/300 * ghgfrac);
  ppX.gas_6(1:100,4) = ppX.gas_6(1:100,4) .* (1 + 6/1800 * ghgfrac);
junk = find(p.landfrac == 1 & p.rlat > +60);
  boo = nanmean(dayfrac*fUMBC_day.results(junk,6)    + nightfrac*fUMBC_night.results(junk,6));        ppX.stemp(5) = ppX.stemp(5) + boo;
  boo = nanmean(dayfrac*fUMBC_day.ptemptrend(:,junk) + nightfrac*fUMBC_night.ptemptrend(:,junk),2);   ppX.ptemp(1:100,5) = ppX.ptemp(1:100,5) + reshape(boo,100,1);
  boo = nanmean(dayfrac*fUMBC_day.gas_1trend(:,junk) + nightfrac*fUMBC_night.gas_1trend(:,junk),2);   ppX.gas_1(1:100,5) = ppX.gas_1(1:100,5) .* (1 + reshape(boo,100,1));
  ppX.gas_2(1:100,5) = ppX.gas_2(1:100,5) .* (1 + 2.2/400 * ghgfrac);
  ppX.gas_4(1:100,5) = ppX.gas_4(1:100,5) .* (1 + 1/300 * ghgfrac);
  ppX.gas_6(1:100,5) = ppX.gas_6(1:100,5) .* (1 + 6/1800 * ghgfrac);

rtpwrite(fop,hh2,hha2,ppX,ppa2);
eval(sartaer);

[~,~,ppradX,~] = rtpread(frp);
simbiasL = rad2bt(hh2.vchan,ppradX.rcalc)-rad2bt(hh2.vchan,pprad0_land.rcalc); 

pprad0_land.plays = plevs2plays(pprad0_land.plevs);
  tafov(3) = nexttile;
  plot(ppradX.ptemp(1:97,1)-pprad0_land.ptemp(1:97,1),pprad0_land.plays(1:97,1),'g',ppradX.ptemp(1:97,2)-pprad0_land.ptemp(1:97,2),pprad0_land.plays(1:97,2),'m--',ppradX.ptemp(1:97,3)-pprad0_land.ptemp(1:97,3),pprad0_land.plays(1:97,3),'r',...
       ppradX.ptemp(1:97,4)-pprad0_land.ptemp(1:97,4),pprad0_land.plays(1:97,4),'c--',ppradX.ptemp(1:97,5)-pprad0_land.ptemp(1:97,5),pprad0_land.plays(1:97,5),'b','linewidth',2); axis([-0.1 +0.1 10 1000])
  if iFrac == -10;  axis([-0.1/2 +0.1/2 10 1000]); end
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); plotaxis2; 
  hl = legend('tropical','S Midlat','N Midlat','S Polar','N Polar','location','best','fontsize',6); 

  tafov(4) = nexttile;
  plot(ppradX.gas_1(1:97,1)./pprad0_land.gas_1(1:97,1),pprad0_land.plays(1:97,1),'g',ppradX.gas_1(1:97,2)./pprad0_land.gas_1(1:97,2),pprad0_land.plays(1:97,2),'m--',ppradX.gas_1(1:97,3)./pprad0_land.gas_1(1:97,3),pprad0_land.plays(1:97,3),'r',...
       ppradX.gas_1(1:97,4)./pprad0_land.gas_1(1:97,4),pprad0_land.plays(1:97,4),'c--',ppradX.gas_1(1:97,5)./pprad0_land.gas_1(1:97,5),pprad0_land.plays(1:97,5),'b','linewidth',2); axis([0.99 1.01 10 1000])
  if iFrac == -10;  axis([1-0.01/2 1+0.01/2 10 1000]); end
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); plotaxis2([1]); 

  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  % Remove all xtick labels except for 2nd row
  for ii = [1 2]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
  end
  for ii = [1 2 3]
   tafov(ii).FontSize = 8;
  end
  ax = gca; ax.XAxis.FontSize = 8; ax.YAxis.FontSize = 8;

disp('LAND : dSTEMP and dMMW for TRP, NML,SML, NP,SP')
fprintf(1,' %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',[ppradX.stemp - pprad0_land.stemp; mmwater_rtp(hh2,ppradX)-mmwater_rtp(hh2,pprad0_land)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(101); clf; 
  ta = tiledlayout(2,1);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  
  tafov(1) = nexttile; 
  plot(hh2.vchan,simbiasO(:,1),'g',hh2.vchan,simbiasO(:,2),'m--',hh2.vchan,simbiasO(:,3),'r',hh2.vchan,simbiasO(:,4),'c--',hh2.vchan,simbiasO(:,5),'b');
  plotaxis2; hl = legend('tropical','S Midlat','N Midlat','S Polar','N Polar','location','best','fontsize',6); axis([640 1620 -0.1 +0.1]); ylabel('OCEAN');
  title('AIRS\_RT')

  tafov(2) = nexttile; 
  plot(hh2.vchan,simbiasL(:,1),'g',hh2.vchan,simbiasL(:,2),'m--',hh2.vchan,simbiasL(:,3),'r',hh2.vchan,simbiasL(:,4),'c--',hh2.vchan,simbiasL(:,5),'b');
  plotaxis2; hl = legend('tropical','S Midlat','N Midlat','S Polar','N Polar','location','best','fontsize',6); axis([640 1620 -0.1 +0.1]); ylabel('LAND')

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

figure(102); 
  ta = tiledlayout(2,1);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  
  tafov(1) = nexttile; 
  plot(hh2.vchan,signal_o(:,1)-simbiasO(:,1),'g',hh2.vchan,signal_o(:,2)-simbiasO(:,2),'m--',hh2.vchan,signal_o(:,3)-simbiasO(:,3),'r',...
       hh2.vchan,signal_o(:,4)-simbiasO(:,4),'c--',hh2.vchan,signal_o(:,5)-simbiasO(:,5),'b'); 
  plotaxis2; hl = legend('tropical','S Midlat','N Midlat','S Polar','N Polar','location','best','fontsize',6); ylabel('OCEAN'); axis([640 1620 -0.03 +0.03]); 
  title('OBS-AIRS\_RT')

  tafov(2) = nexttile; 
  plot(hh2.vchan,signal_l(:,1)-simbiasL(:,1),'g',hh2.vchan,signal_l(:,2)-simbiasL(:,2),'m--',hh2.vchan,signal_l(:,3)-simbiasL(:,3),'r',...
       hh2.vchan,signal_l(:,4)-simbiasL(:,4),'c--',hh2.vchan,signal_l(:,5)-simbiasL(:,5),'b'); 
  plotaxis2; hl = legend('tropical','S Midlat','N Midlat','S Polar','N Polar','location','best','fontsize',6); ylabel('LAND');  axis([640 1620 -0.03 +0.03])

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

%if abs(iFrac) > 1
%  figure(101); axis([640 1620 -0.02 +0.01]);
%  figure(103); axis([640 1620 -0.02 +0.01]);
%end
