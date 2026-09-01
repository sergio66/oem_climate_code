%if sum(maskLF) == 4608
%  figure(8); clf; plot(f,nanmean(rates,2)-nanmean(rates.*maskLFchan,2),'b',f,nanmean(rates,2)-nanmean(rates'.*maskLFchan')','c',f,nanmean(rates,2)-nanmean(rates(:,mask),2),'g',f,nanmean(rates,2)-nanmean(rates(:,maskLF),2),'r')
%  title('RED is BAD way to do things (maskLF)')
%end

figure(8); plot(f,nanmean(rates'.*maskLFchan'),'r',f,nanstd(rates'.*maskLFchan'),'m',f,nanmean(rates'.*maskLFchan'-fits'.*maskLFchan'),'b',f,nanstd(rates'.*maskLFchan'-fits'.*maskLFchan'),'c'); 
  grid; hl = legend('mean(signal)','std(signal)','mean(signal-fit)','std(signal-fit)','location','best','fontsize',10);
  xlim([min(f) max(f)])
figure(8); plot(f,nanmean(rates'.*maskLFchan'),'r',f,nanstd(rates'.*maskLFchan'),'m',f,nanmean(fits'.*maskLFchan'),'b',f,nanstd(fits'.*maskLFchan'),'c'); 
  grid; hl = legend('mean(signal)','std(signal)','mean(fit)','std(fit)','location','best','fontsize',10);
  xlim([min(f) max(f)])

  ta = tiledlayout(2,1);
  % nexttile says go to the obvious, save it’s handle

  taHandle(1)= nexttile;
  plot(f,nanmean(rates'.*maskLFchan'),'r',f,nanstd(rates'.*maskLFchan'),'m',f,nanmean(fits'.*maskLFchan'),'b',f,nanstd(fits'.*maskLFchan'),'c'); 
  grid; hl = legend('mean(signal)','std(signal)','mean(fit)','std(fit)','location','best','fontsize',10);
  xlim([min(f) max(f)])

  taHandle(2)= nexttile;
  plot(f,nanmean(rates'.*maskLFchan'-fits'.*maskLFchan'),'b',f,nanstd(rates'.*maskLFchan'-fits'.*maskLFchan'),'c'); 
  grid; hl = legend('mean(signal-fit)','std(signal-fit)','location','best','fontsize',10);
  xlim([min(f) max(f)])

  % Get rid of all extra space I can
  ta.Padding = 'none';
  ta.TileSpacing = 'none';  %% or compact

  %taHandle(1).XTickLabel = '';
  %taHandle(1).YLabel.String = 'BT(K)';
  %taHandle(2).YLabel.String = 'BT(K)';
  %taHandle(2).XLabel.String = '\nu cm^{-1}';

  %% https://www.mathworks.com/help/matlab/ref/linkaxes.html
  %% synchronize
  linkaxes(taHandle,'xy')
  %% unsynchronize
  linkaxes(taHandle,'off')
  
  linkaxes(taHandle,'x');
  xlabel(ta,'\nu cm^{-1}');
  ylabel(ta,'BT(K)');
  xticklabels(taHandle(1),{})
  xlim([640 1640])
