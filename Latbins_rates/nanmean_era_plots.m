
figure(6); clf
  subplot(121)
  shadedErrorBarY(nanmean(water(lat_index,:)),plays,nanstd(water(lat_index,:)),'bo-',1);
  hold on
  shadedErrorBarY(nanmean(waterrate(lat_index,:)),plays,nanstd(waterrate(lat_index,:)),'rx-',1);
  hold off; title('AIRS(b) ERA(r) Water frac/yr'); %set(hl,'fontsize',10); grid
  set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 0 1000]); 

  subplot(122)
  shadedErrorBarY(nanmean(temp(lat_index,:)),plays,nanstd(temp(lat_index,:)),'bo-',1);
  hold on
  shadedErrorBarY(nanmean(ptemprate(lat_index,:)),plays,nanstd(ptemprate(lat_index,:)),'rx-',1);

  hl=errorbar_x(nanmean(nTMT.trend(oink)),xPRESS_tmt,nanstd(nTMT.trend(oink)),'ko'); set(hl,'linewidth',2)
  hl=errorbar_x(nanmean(nTUT.trend(oink)),xPRESS_ttt,nanstd(nTUT.trend(oink)),'ko'); set(hl,'linewidth',2)
  hl=errorbar_x(nanmean(nTLS.trend(oink)),xPRESS_tls,nanstd(nTLS.trend(oink)),'ko'); set(hl,'linewidth',2)

  %these were 30 yr trends
  %oink = find(abs(rssTLT(:,1)) < 30);   errorbar_x(nanmean(rssTLT(oink,2))/10,xPRESS_tlt,nanstd(rssTLT(oink,2)),'k*')
  %oink = find(abs(rssTMT(:,1)) < 30);   errorbar_x(nanmean(rssTMT(oink,2))/10,xPRESS_tmt,nanstd(rssTMT(oink,2)),'k*')
  %oink = find(abs(rssTTT(:,1)) < 30);   errorbar_x(nanmean(rssTTT(oink,2))/10,xPRESS_ttt,nanstd(rssTTT(oink,2)),'k*')
  %oink = find(abs(rssTLS(:,1)) < 30);   errorbar_x(nanmean(rssTLS(oink,2))/10,xPRESS_tls,nanstd(rssTLS(oink,2)),'k*')
  %oink = find(abs(rssTTS(:,1)) < 30);   errorbar_x(nanmean(rssTTS(oink,2))/10,xPRESS_tts,nanstd(rssTTS(oink,2)),'k*')

  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(1),xPRESS_tls,rssTRP.amsu_rss_error_robust(1),'g*'); set(hl,'linewidth',2)
  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(2),xPRESS_tts,rssTRP.amsu_rss_error_robust(2),'g*'); set(hl,'linewidth',2)
  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(3),xPRESS_tmt,rssTRP.amsu_rss_error_robust(3),'g*'); set(hl,'linewidth',2)
  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(4),xPRESS_tlt,rssTRP.amsu_rss_error_robust(4),'g*'); set(hl,'linewidth',2)

  hold off; title('AIRS(b) ERA(r) Temp K/yr'); %set(hl,'fontsize',10); grid
  set(gca,'ydir','reverse'); grid; axis([-0.10 +0.10 0 1000]);

figure(4); clf
  pcolor(save_lat,plays,double(waterrate')); shading flat; colorbar; 
  set(gca,'ydir','reverse'); grid; title('ERA WV frac/yr ')
figure(5); clf
  pcolor(save_lat,plays,double(ptemprate')); shading flat; colorbar; 
  set(gca,'ydir','reverse'); grid; title('ERA T K/yr ')

figure(1); caxis([-0.02 +0.02]); shading flat; colorbar
figure(4); caxis([-0.02 +0.02]); shading flat; colorbar

figure(2); caxis([-0.1 +0.1]); colorbar  
figure(5); caxis([-0.1 +0.1]); colorbar

figure(8);
  shadedErrorBar(save_lat,params(:,6),params_sigs(:,6),'bo-',1);
  hold on
  shadedErrorBar(save_lat,stemprate,stempratestd,'ro-',1);
  hold off
grid
title('AIRS(b) ERA(r) STemp K/yr'); %set(hl,'fontsize',10); grid

%figure(4)
%hold off
%  set(gca,'ydir','reverse'); grid; title('T K/yr')
%  hl = legend('OEM','ERA','location','east'); set(hl,'fontsize',10)

figure(9); clf
  shadedErrorBarY(nanmean(temp(lat_index,:)),log10(plays),nanstd(temp(lat_index,:)),'bo-',1);
  hold on
  shadedErrorBarY(nanmean(ptemprate(lat_index,:)),log10(plays),nanstd(ptemprate(lat_index,:)),'rx-',1);

  lxPRESS_tmt = log10(xPRESS_tmt);
  lxPRESS_ttt = log10(xPRESS_ttt);
  lxPRESS_tlt = log10(xPRESS_tlt);
  lxPRESS_tts = log10(xPRESS_tts);
  lxPRESS_tls = log10(xPRESS_tls);

  hl=errorbar_x(nanmean(nTMT.trend(oink)),lxPRESS_tmt,nanstd(nTMT.trend(oink)),'ko'); set(hl,'linewidth',2)
  hl=errorbar_x(nanmean(nTUT.trend(oink)),lxPRESS_ttt,nanstd(nTUT.trend(oink)),'ko'); set(hl,'linewidth',2)
  hl=errorbar_x(nanmean(nTLS.trend(oink)),lxPRESS_tls,nanstd(nTLS.trend(oink)),'ko'); set(hl,'linewidth',2)

  %these were 30 yr trends
  %oink = find(abs(rssTLT(:,1)) < 30);   errorbar_x(nanmean(rssTLT(oink,2))/10,lxPRESS_tlt,nanstd(rssTLT(oink,2)),'k*')
  %oink = find(abs(rssTMT(:,1)) < 30);   errorbar_x(nanmean(rssTMT(oink,2))/10,lxPRESS_tmt,nanstd(rssTMT(oink,2)),'k*')
  %oink = find(abs(rssTTT(:,1)) < 30);   errorbar_x(nanmean(rssTTT(oink,2))/10,lxPRESS_ttt,nanstd(rssTTT(oink,2)),'k*')
  %oink = find(abs(rssTLS(:,1)) < 30);   errorbar_x(nanmean(rssTLS(oink,2))/10,lxPRESS_tls,nanstd(rssTLS(oink,2)),'k*')
  %oink = find(abs(rssTTS(:,1)) < 30);   errorbar_x(nanmean(rssTTS(oink,2))/10,lxPRESS_tts,nanstd(rssTTS(oink,2)),'k*')

  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(1),lxPRESS_tls,rssTRP.amsu_rss_error_robust(1),'g*'); set(hl,'linewidth',2)
  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(2),lxPRESS_tts,rssTRP.amsu_rss_error_robust(2),'g*'); set(hl,'linewidth',2)
  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(3),lxPRESS_tmt,rssTRP.amsu_rss_error_robust(3),'g*'); set(hl,'linewidth',2)
  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(4),lxPRESS_tlt,rssTRP.amsu_rss_error_robust(4),'g*'); set(hl,'linewidth',2)

  hold off; title('AIRS(b) ERA(r) Temp K/yr'); %set(hl,'fontsize',10); grid
  set(gca,'ydir','reverse'); grid; axis([-0.10 +0.10 1 3]);

figure(9); clf
  shadedErrorBarY(nanmean(ozonerate(lat_index,:)),(plays),nanstd(ozonerate(lat_index,:)),'rx-',1);
  hold on
  orate = nanmean(ozonerate(lat_index,:)); 
  prate = plays;  
    prate = sum(abs(orate).*prate')/sum(abs(orate)); orate = nanmean(orate);
  hl=errorbar_x(orate,prate,orate*2,'rx');
  set(hl,'linewidth',2)
  hl=errorbar_x(nanmean(tracegases(lat_index,2)),200,nanstd(tracegases(lat_index,2)),'bo');
  set(hl,'linewidth',2)
  hold off; title('AIRS(b) ERA(r) Ozone frac/yr'); %set(hl,'fontsize',10); 
  set(gca,'ydir','reverse'); grid; axis([-0.05 +0.05 0 1000]); grid

figure(7); clf
  hold on
  plot(ff(g),input_rates(lat_index,g),'b',ff(g),fitted_rates(lat_index,g),'r');   
  plot(ff(chanset),input_rates(lat_index,chanset),'bo',ff(chanset),fitted_rates(lat_index,chanset),'ro');   
  hold off
grid

%{
chanset = jacobian.chanset;
%g = dogoodchan;
figure(5);
plot(f,input_rates(18,:),'b',f,fitted_rates(18,:),'r',...
      f,input_rates(18,:)-fitted_rates(18,:),'k',...
      f(chanset),input_rates(18,chanset)-fitted_rates(18,chanset),'ko')
  hl = legend('input','fit','bias=input-fit','location','north'); set(hl,'fontsize',10)
  axis([640 2780 -0.10 +0.10]); grid

wvrates.oem = wv(18,:);
wvrates.oem_d = dwv(18,:);
wvrates.era = nanmean(waterrate(lat_index,:));
wvrates.plays = plays;

trates.oem = t(18,:);
trates.oem_d = dt(18,:);
trates.era = nanmean(ptemprate(lat_index,:));
trates.plays = plays;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'rateset = %s \n',rateset.datafile);
fprintf(1,'jacset  = %s \n',jacobian.filename);

%}
