disp('Fig 21 has to be handstrtched so it is as wide as Fig 4 + Fig 20')
disp('Fig 21 has to be handstrtched so it is as wide as Fig 4 + Fig 20')
disp('Fig 21 has to be handstrtched so it is as wide as Fig 4 + Fig 20')

figure(21); clf; 
ta = tiledlayout(2,2);

% Tile 1
% Span across two rows and columns, this is Fig 4
nexttile([2 1]);
%plot(h.vchan,q_mean_BT); xlabel('Wavenumber [cm^{-1}]','Interpreter', 'latex'); ylabel('BT [K]');
plot(h.vchan,q_mean_BT); xlabel('Wavenumber [cm^{-1}]'); ylabel('BT [K]'); 
xlim([640 1620]); plotaxis2; hl = legend('Q50','Q80','Q90','Q95','Q97','location','best','fontsize',12);
ylim([200 300]);
set(gca,'fontsize',14)

%% this is Fig 20
tafov(1) = nexttile; plot(h.vchan,czoo01,h.vchan,czoo02,h.vchan,czoo03,h.vchan,czoo04,h.vchan,czoo05); plotaxis2;
  %xlim([640 1620]); ylabel('d(BT)/dt [K yr^{-1}]','Interpreter', 'latex');
  xlim([640 1620]); ylabel('d(BT)/dt [K yr^{-1}]');  
  %legend('Q50','Q80','Q90','Q95','Q97','fontsize',10,'location','south'); 
  %legend('Q50','Q80','Q90','Q95','Q97','fontsize',10,'Position',[0.45 0.45 0.25 0.25])
  set(gca,'fontsize',14)
tafov(2) = nexttile; plot(h.vchan,cuoo01,h.vchan,cuoo02,h.vchan,cuoo03,h.vchan,cuoo04,h.vchan,cuoo05); plotaxis2;
  xlim([640 1620]); 
  %legend('Q50','Q80','Q90','Q95','Q97','fontsize',8,'location','southeast'); 
  %legend('Q50','Q80','Q90','Q95','Q97','fontsize',10,'Position',[0.45 0.4625 0.2 0.2]) %% july 2025
  legend('Q50','Q80','Q90','Q95','Q97','fontsize',12,'Position',[0.60 0.4625 0.2 0.2]) %% july 2025  
  %xlabel('Wavenumber [cm^{-1}]','Interpreter', 'latex'); ylabel('d(BT)/dt [K yr^{-1}]','Interpreter', 'latex');
  xlabel('Wavenumber [cm^{-1}]'); ylabel('d(BT)/dt [K yr^{-1}]');  
  ylim([0 0.031])
  set(gca,'fontsize',14)

ta.Padding = 'none';
ta.TileSpacing = 'compact';
ta.Padding = 'compact';
ta.TileSpacing = 'tight';
tafov(1).XTickLabel = '';

%% sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/trend_dBTdt_desc_20_years_Q50_Q80_Q90_Q95_Q97_bigfont_Aug2025')
