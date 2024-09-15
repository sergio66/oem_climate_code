cos64lat = cos(rlat*pi/180) * ones(1,iNumAnomTimeSteps);
scalebands = [0 10 10/2 10/5 10/2 10/2 10 10 25 25 100 100 100 300];
smn = 23; %% this is 1 year,   23 * 16 = 368 days
smn = 06; %% this is 1/4 year, 06 * 16 =  96 days 
smn = 03; %% this is 1/8 year, 03 * 16 =  48 days 

for ii = 2 : length(RRTM_bands)
  figure(1); plot(yymm,smooth(squeeze(nanmean(cos64lat.*fluxanomD(:,:,1),1)),01),'b',yymm,smooth(squeeze(nanmean(cos64lat.*fluxanomD(:,:,ii),1)),01)*scalebands(ii),'r');
  figure(1); plot(yymm,smooth(squeeze(nanmean(cos64lat.*fluxanomA(:,:,1),1)),smn),'bx-',yymm,smooth(squeeze(nanmean(cos64lat.*fluxanomA(:,:,ii),1)),smn)*scalebands(ii),'c',...
                  yymm,smooth(squeeze(nanmean(cos64lat.*fluxanomD(:,:,1),1)),smn),'rx-',yymm,smooth(squeeze(nanmean(cos64lat.*fluxanomD(:,:,ii),1)),smn)*scalebands(ii),'m','linewidth',2);
             plotaxis2; title(['ii=' num2str(ii) ' : ' num2str(RRTM_bands(ii-1)) '-' num2str(RRTM_bands(ii)) ' cm-1'])
             ax = axis; line([2022+15/30/12  2022+15/30/12],[ax(3) ax(4)],'color','k','linewidth',2);
             hl = legend('A all','A band','D all','D band','location','best','fontsize',8);
             xlim([2020 2025])
  pause
end
