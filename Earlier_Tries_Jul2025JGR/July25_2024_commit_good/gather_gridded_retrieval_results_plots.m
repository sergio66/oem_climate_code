maskLFchan = ones(2645,1) * maskLF;
mask = find(maskLF == 1);

redo_fig8_spectralrates_with_mask

i730 = find(f >= 730,1);   figure(9);  clf; scatter_coast(X(:),Y(:),50,maskLF.*(rates(i730,:)-fits(i730,:))); title('BT730 bias')
i900 = find(f >= 900,1);   figure(10); clf; scatter_coast(X(:),Y(:),50,maskLF.*(rates(i900,:)-fits(i900,:))); title('BT900 bias')
i1419 = find(f >= 1419,1); figure(11); clf; scatter_coast(X(:),Y(:),50,maskLF.*(rates(i1419,:)-fits(i1419,:))); title('BT1419 bias')
i667 = find(f >= 667,1);   figure(12); clf; scatter_coast(X(:),Y(:),50,maskLF.*(rates(i667,:)-fits(i667,:))); title('BT667 bias')
for ii = 9:12; figure(ii); caxis([-1 +1]*0.01); colormap(llsmap5); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% have added 4 new figures, so these subsequent figures displaced by 4 eg fig 11 --> fig 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chisqr = rates(jacobian.chanset,:)-fits(jacobian.chanset,:); chisqr = sum(chisqr.*chisqr,1)/length(jacobian.chanset); figure(13); clf; scatter_coast(X(:),Y(:),50,log10(chisqr.*maskLF)); title('log10(\chi^2) oem chans')

junkBT  = rates(i900,:).*maskLF;
junkT20 = resultsT(:,iNumLay)'.*maskLF;
junkW20 = resultsWV(:,iNumLay)'.*maskLF;
figure(14); clf; plot(rlat,nanmean(reshape(junkBT,72,64),1),'g',rlat,10*nanmean(reshape(junkW20,72,64),1),'bx-',...
                      rlat,nanmean(reshape(junkT20,72,64),1),'k',rlat,nanmean(reshape(results(:,6).*maskLF',72,64),1),'r','linewidth',2); title('zonal rates LAND+OCEAN'); plotaxis2; ylim([-0.04 +0.04]); xlim([-60 +60])
hl = legend('BT900 obs','10*W20','T20','SST','location','best'); 

land = (landfrac > 0); 
fprintf(1,'Out of %5i grid points, the fraction over land is %8.6f \n',length(landfrac),sum(land)/length(landfrac))

resultsNAN = results;               resultsNAN(land,:) = NaN;
junkBTNAN  = rates(i900,:);         junkBTNAN(land) = NaN;
junkT20NAN = resultsT(:,iNumLay)';  junkT20NAN(land) = NaN;
junkW20NAN = resultsWV(:,iNumLay)'; junkT20NAN(land) = NaN;
figure(15); clf; plot(rlat,nanmean(reshape(junkBTNAN,72,64),1),'g',rlat,10*nanmean(reshape(junkW20NAN,72,64),1),'bx-',...
                      rlat,nanmean(reshape(junkT20NAN,72,64),1),'k',rlat,nanmean(reshape(resultsNAN(:,6),72,64),1),'r','linewidth',2); title('zonal rates OCEAN'); plotaxis2; ylim([-0.04 +0.04]); xlim([-60 +60])
hl = legend('BT900 obs','10*W20','T20','SST','location','best'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% these are NOT masked for L/O
figure(16); clf; pcolor(rlon,pavg,resultsWV((1:72)+30*72,:)'); xlabel('longitude'); title('Equator WV'); caxis([-1e-2 +1e-2]);      
figure(17); clf; pcolor(rlon,pavg,resultsT((1:72)+30*72,:)');  xlabel('longitude'); title('Equator T');  caxis([-0.1 +0.1]);        
figure(18); clf; pcolor(rlon,pavg,resultsO3((1:72)+30*72,:)'); xlabel('longitude'); title('Equator O3'); caxis([-0.01 +0.01]*0.25); 
for ii = 16:18; figure(ii); colorbar; shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); end

%% these are NOT masked for L/O
ind = 1:72:4608;
ind = ind + 31;
figure(19); clf; pcolor(rlat,pavg,resultsWV(ind,:)'); xlabel('latitude'); title('GMT line WV'); caxis([-1e-2 +1e-2]);      
figure(20); clf; pcolor(rlat,pavg,resultsT(ind,:)');  xlabel('latitude'); title('GMT line T');  caxis([-0.1 +0.1]);        
figure(21); clf; pcolor(rlat,pavg,resultsO3(ind,:)'); xlabel('latitude'); title('GMT line O3'); caxis([-0.01 +0.01]*0.25); 
for ii = 19:21; figure(ii); colorbar; shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); end

boo = zeros(iNumLay,72,64); for ijunk = 1 : iNumLay; boo(ijunk,:,:) = maskLFmatr'; end
xresultsWV = resultsWV'; xresultsWV = boo .* reshape(xresultsWV,iNumLay,72,64);
xresultsT  = resultsT';   xresultsT = boo .* reshape(xresultsT,iNumLay,72,64);
xresultsO3 = resultsO3'; xresultsO3 = boo .* reshape(xresultsO3,iNumLay,72,64);
figure(22); clf; pcolor(rlat,pavg,squeeze(nanmean(xresultsWV,2))); xlabel('latitude'); title('zonal mean WV'); caxis([-1e-2 +1e-2]);      
figure(23); clf; pcolor(rlat,pavg,squeeze(nanmean(xresultsT,2)));  xlabel('latitude'); title('zonal mean T');  caxis([-0.1 +0.1]);        
figure(24); clf; pcolor(rlat,pavg,squeeze(nanmean(xresultsO3,2))); xlabel('latitude'); title('zonal mean O3'); caxis([-0.01 +0.01]*0.25);
for ii = 22:24; figure(ii); ylim([10 1000]); xlim([-90 +90]); colorbar; shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iSavePrintBasics_gather_gridded_retrieval_results_plots

