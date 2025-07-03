%% see /home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/make_fig4_matfile.m

clear all
load fig4.mat

iiMax = 4;

figure(1); clf
for ii = 1 : iiMax
  plot(data.x{ii},data.y{ii},'linewidth',2); hold on
end

figure(1); clf
ii = 4;  plot(data.x{ii},data.y{ii},'linewidth',3); hold on
ii = 3;  plot(data.x{ii},data.y{ii},'linewidth',2); hold on
ii = 2;  plot(data.x{ii},data.y{ii},'linewidth',2); hold on
ii = 1;  plot(data.x{ii},data.y{ii},'linewidth',2); hold on

hold off
legend('AIRS L1C (observed)','AIRS v7 (closure)','CLIMCAPS (closure)','ERA5 (closure)','location','best'); 
grid on
xlabel('Wavenumber (cm^{-1}'); ylabel('dBT/dt (K/yr)');
set(gca,'fontsize',14)

axis(ax)
ylim([-0.025 +0.02])

%pos = [2700 500 700 500];
%set(gcf,'position',pos)
