airsL1Crates = load('/asl/s1/sergio/JUNK/test8junk3.mat','rates','nedt');
era5rates    = load('/asl/s1/sergio/JUNK/test7_retrieve_syntheticERA5_march14_2023v2.mat','rates','nedt');

if ~exist('f','var')
  junk = load('h2645structure.mat');
  f    = junk.h.vchan;
end

figure(70); plot(f,nanmean(airsL1Crates.rates'),'b',f,+nanmean(airsL1Crates.nedt'),'c',f,-nanmean(airsL1Crates.nedt'),'c',...
                 f,nanmean(era5rates.rates'),'r',f,+nanmean(era5rates.nedt'),'m',f,-nanmean(era5rates.nedt'),'m')
plotaxis2; xlim([645 1645]); hl = legend('AIRS L1C rates','AIRS L1C noise','AIRS L1C noise','ERA5 rates','ERA5 noise','ERA5 noise','location','best');

g4 = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g4_jac.mat');
g6 = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g6_jac.mat');

hold on; plot(g4.fout,sum(g4.jout')/100,'k',g6.fout,sum(g6.jout')/100,'g'); hold off

figure(71);
i820 = find(f >= 822.4,1);
i961 = find(f >= 961.1,1);
slope961_820_obs = airsL1Crates.rates(i961,:)-airsL1Crates.rates(i820,:);
slope961_820_era = era5rates.rates(i961,:)-era5rates.rates(i820,:);
plot(1:64,nanmean(reshape(slope961_820_obs,72,64),1),1:64,nanmean(reshape(slope961_820_era,72,64),1))

aslmap(71,rlat65,rlon73,smoothn((reshape(slope961_820_obs,72,64)') ,1), [-90 +90],[-180 +180]); title('Obs Slope 820 - 961 cm-1'); caxis([-1 +1]*1.5e-2); colormap(llsmap5)
aslmap(72,rlat65,rlon73,smoothn((reshape(slope961_820_era,72,64)') ,1), [-90 +90],[-180 +180]); title('Era Slope 820 - 961 cm-1'); caxis([-1 +1]*1.5e-2); colormap(llsmap5)
