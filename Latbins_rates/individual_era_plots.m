iLinearOrLog = +1;   %% linear plot
iLinearOrLog = -1;   %% log plot

figure(6); clf
if iLinearOrLog == +1 
  subplot(121)
  shadedErrorBarY(water(lat_index,:),plays,water_sigs(lat_index,:),'bo-',1);
  hold on
  shadedErrorBarY(waterrate(lat_index,:),plays,waterratestd(lat_index,:),'rx-',1);
  hold off; hl = title('AIRS(b) ERA(r) Water frac/yr'); set(hl,'fontsize',10); 
  set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 0 1000]); 

  subplot(122)
  shadedErrorBarY(temp(lat_index,:),plays,temp_sigs(lat_index,:),'bo-',1);
  hold on
  shadedErrorBarY(ptemprate(lat_index,:),plays,ptempratestd(lat_index,:),'rx-',1);

  hl=errorbar_x(nTMT.trend(oink),xPRESS_tmt,nTMT.err_trend(oink),'ko'); set(hl,'linewidth',2)
  hl=errorbar_x(nTUT.trend(oink),xPRESS_ttt,nTUT.err_trend(oink),'ko'); set(hl,'linewidth',2)
  hl=errorbar_x(nTLS.trend(oink),xPRESS_tls,nTLS.err_trend(oink),'ko'); set(hl,'linewidth',2)

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

  hold off; hl = title('AIRS(b) ERA(r) Temp K/yr'); set(hl,'fontsize',10);
  set(gca,'ydir','reverse'); grid; axis([-0.25 +0.25 0 1000]);

else
  subplot(121)
  shadedErrorBarYLog10(water(lat_index,:),plays,water_sigs(lat_index,:),'bo-');
  hold on
  shadedErrorBarYLog10(waterrate(lat_index,:),plays,waterratestd(lat_index,:),'rx-')
  hold off; hl = title('AIRS(b) ERA(r) Water frac/yr'); set(hl,'fontsize',10); 
  set(gca,'ydir','reverse'); grid; axis([-0.025 +0.025 0 1000]); 

  subplot(122)
  shadedErrorBarYLog10(temp(lat_index,:),plays,temp_sigs(lat_index,:),'bo-');
  hold on
  shadedErrorBarYLog10(ptemprate(lat_index,:),plays,ptempratestd(lat_index,:),'rx-');

%  shadedErrorBarY(nTMT.trend(oink),(xPRESS_tmt),nTMT.err_trend(oink),'ko',1,-1);
%  hlxx=errorbar_x(nTMT.trend(oink),(xPRESS_tmt),nTMT.err_trend(oink),'ko'); set(hlxx,'linewidth',2); 
%  set(hlxx,'YScale','log')
%  hl=errorbar_x(nTUT.trend(oink),log10(xPRESS_ttt),nTUT.err_trend(oink),'ko'); set(hl,'linewidth',2)
%  hl=errorbar_x(nTLS.trend(oink),log10(xPRESS_tls),nTLS.err_trend(oink),'ko'); set(hl,'linewidth',2)

  %these were 30 yr trends
  %oink = find(abs(rssTLT(:,1)) < 30);   errorbar_x(nanmean(rssTLT(oink,2))/10,xPRESS_tlt,nanstd(rssTLT(oink,2)),'k*')
  %oink = find(abs(rssTMT(:,1)) < 30);   errorbar_x(nanmean(rssTMT(oink,2))/10,xPRESS_tmt,nanstd(rssTMT(oink,2)),'k*')
  %oink = find(abs(rssTTT(:,1)) < 30);   errorbar_x(nanmean(rssTTT(oink,2))/10,xPRESS_ttt,nanstd(rssTTT(oink,2)),'k*')
  %oink = find(abs(rssTLS(:,1)) < 30);   errorbar_x(nanmean(rssTLS(oink,2))/10,xPRESS_tls,nanstd(rssTLS(oink,2)),'k*')
  %oink = find(abs(rssTTS(:,1)) < 30);   errorbar_x(nanmean(rssTTS(oink,2))/10,xPRESS_tts,nanstd(rssTTS(oink,2)),'k*')

%  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(1),xPRESS_tls,rssTRP.amsu_rss_error_robust(1),'g*'); set(hl,'linewidth',2)
%  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(2),xPRESS_tts,rssTRP.amsu_rss_error_robust(2),'g*'); set(hl,'linewidth',2)
%  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(3),xPRESS_tmt,rssTRP.amsu_rss_error_robust(3),'g*'); set(hl,'linewidth',2)
%  hl=errorbar_x(rssTRP.amsu_rss_rate_robust(4),xPRESS_tlt,rssTRP.amsu_rss_error_robust(4),'g*'); set(hl,'linewidth',2)

  hold off; hl = title('AIRS(b) ERA(r) Temp K/yr'); set(hl,'fontsize',10);
  set(gca,'ydir','reverse'); grid; axis([-0.25 +0.25 0 1000]);
end

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
  shadedErrorBarY(ozonerate(lat_index,:),(plays),ozoneratestd(lat_index,:),'rx-',1);
  hold on
  orate = ozoneratestd(lat_index,:); 
  prate = plays;  
    prate = sum(abs(orate).*prate')/sum(abs(orate)); orate = nanmean(orate);
  hl=errorbar_x(orate,prate,orate*2,'rx');
  set(hl,'linewidth',2)
  hl=errorbar_x(tracegases(lat_index,2),200,tracegases_sigs(lat_index,2),'bo');
  set(hl,'linewidth',2)
  hold off; hl = title('AIRS(b) ERA(r) Ozone frac/yr'); set(hl,'fontsize',10); 
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
