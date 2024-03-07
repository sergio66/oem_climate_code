addpath /asl/matlib/h4tools
addpath /home/sergio/KCARTA/MATLAB

%{
DAN ZHOU
/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/set_driver_jacfile.m tells me to look at 
/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Mar2023_startSept2002_endAug2022_trendsonly/clust_put_together_jacs_clrERA5.m which uses
/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.rp.rtp
which is symbolic link to 
/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.rp.rtp

CAMEL
/home/sergio/MATLABCODE_Git/CAMEL_emissivity/Trends_camelV003_emis4608tiles/emistrendsV003.mat
%}

sarta = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';
sartaer = ['!time ' sarta ' fin=junk.op.rtp fout=junk.rp.rtp listj=100'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,ha,p0,pa] = rtpread('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.rp.rtp');
p = p0;
rtpwrite('junk.op.rtp',h,ha,p,pa);
eval(sartaer);
[h,ha,pZhou,pa] = rtpread('junk.rp.rtp');
[w,jacZhou,iaProf,iaNumLay] = readsarta_jac('junk.rp.rtp_jacTZ',100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exists('camel')
  camel = load('/home/sergio/MATLABCODE_Git/CAMEL_emissivity/Trends_camelV003_emis4608tiles/emistrendsV003.mat');
end
p = p0;
p.nemis = nanmean(camel.new_nemisave,1);
p.efreq = camel.new_efreqsave;
p.emis  = squeeze(nanmean(camel.new_emissave,1));
p.rho = (1-p.emis)/pi;
compare_two_structures(p0,p)

for ii = 1 : 4608
  zhou_emiss_100(:,ii) = interp1(pZhou.efreq(1:pZhou.nemis(ii),ii),pZhou.emis(1:pZhou.nemis(ii),ii),p.efreq(:,ii));
end
plot(nanmean(p0.efreq'),nanmean(p0.emis'),nanmean(p.efreq'),nanmean(p.emis'),'r'); hl = legend('Zhou','Camel','location','best'); title('blue is wierd because of different nemis')
plot(nanmean(p.efreq'),nanmean(zhou_emiss_100'),nanmean(p.efreq'),nanmean(p.emis'),'r'); hl = legend('Zhou','Camel','location','best');

rtpwrite('junk.op.rtp',h,ha,p,pa);
eval(sartaer);
[h,ha,pCamel,pa] = rtpread('junk.rp.rtp');
[w,jacCamel,iaProf,iaNumLay] = readsarta_jac('junk.rp.rtp_jacTZ',100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tZhou = rad2bt(h.vchan,pZhou.rcalc);
tCamel = rad2bt(h.vchan,pCamel.rcalc);
plot(h.vchan,nanmean(tZhou'-tCamel'),h.vchan,nanstd(tZhou'-tCamel'))
plotaxis2;

jacSTCamel = squeeze(jacCamel(:,:,101))';
jacSTZhou  = squeeze(jacZhou(:,:,101))';
for ii = 1 : 4608
  jacSTCamel(:,ii) = squeeze(jacCamel(ii,:,p.nlevs(ii)))';
  jacSTZhou(:,ii)  = squeeze(jacZhou(ii,:,p.nlevs(ii)))';
end
wahZ = squeeze(jacZhou(2300,:,p.nlevs(2300)));
wahC = squeeze(jacCamel(2300,:,p.nlevs(2300)));

figure(1); clf
plot(h.vchan,nanmean(jacSTZhou'-jacSTCamel'),'b',h.vchan,nanstd(jacSTZhou'-jacSTCamel'),'c--',h.vchan,nanmean(camel.t00'-camel.t0X'),'r',h.vchan,nanstd(camel.t00'-camel.t0X'),'m--'); plotaxis2;
hl = legend('bias jac Zhou-Camel','std jac Zhou-Camel','mean BT trend Camel emiss changes','std BT trend Camel emiss changes','location','best','fontsize',8); title('Global')
xlim([640 1640])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%bad asc 20 year trends
plot_72x64_tiles([1434 2419 2491 2716 2787 2788 2859 2860 2392 3004])

%bad asc 20 year trends over India
plot_72x64_tiles([2716 2787 2788 2859 2860])
plot_72x64_tiles([[2571:2573] [2643:2645] [2715:2717] [2787:2789]])
plot_72x64_tiles([[2571:2573] [2643:2645] [2715:2717] [2787:2789] [2859:2862] [2929:2935] [3000:3008]])
  colormap(usa2); colorbar('off')

miaow = load('/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_day_removeemisstrends.mat','results');
miaow = load('/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_day.mat','results');
plot_72x64_tiles([[2571:2573] [2643:2645] [2715:2717] [2787:2789] [2859:2862] [2929:2935] [3000:3008]],miaow.results(:,6))
plot_72x64_tiles([[2571:2573] [2643:2645] [2715:2717] [2787:2789]],miaow.results(:,6))
%}

figure(2); clf
ind = [[2571:2573] [2643:2645] [2715:2717] [2787:2789]];
plot(h.vchan,nanmean(jacSTZhou(:,ind)'-jacSTCamel(:,ind)'),'b',h.vchan,nanstd(jacSTZhou(:,ind)'-jacSTCamel(:,ind)'),'c--',...
     h.vchan,nanmean(camel.t00(:,ind)'-camel.t0X(:,ind)'),'r',h.vchan,nanstd(camel.t00(:,ind)'-camel.t0X(:,ind)'),'m--'); plotaxis2;
hl = legend('bias jac Zhou-Camel','std jac Zhou-Camel','mean BT trend Camel emiss changes','std BT trend Camel emiss changes','location','best','fontsize',8); title('India')
xlim([640 1640])

addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load llsmap5
i900 = find(h.vchan >= 900,1);
i1231 = find(h.vchan >= 1231,1);
plot_72x64_tiles([],(camel.t00(i900,:)-camel.t0X(i900,:))); colormap(llsmap5); caxis([-1 +1]*0.15)
figure(3); clf; plot_72x64_tiles([[2571:2573] [2643:2645] [2715:2717] [2787:2789]],camel.t00(i900,:)-camel.t0X(i900,:)); colormap(llsmap5); caxis([-1 +1]*0.05)
figure(4); clf; plot_72x64_tiles([[2571:2573] [2643:2645] [2715:2717] [2787:2789]],jacSTZhou(i900,:)-jacSTCamel(i900,:)); colormap(llsmap5); caxis([-1 +1]*0.05)

figure(5); clf
plot(nanmean(p.efreq(:,ind)'),nanmean(zhou_emiss_100(:,ind)'),nanmean(p.efreq(:,ind)'),nanmean(p.emis(:,ind)'),'r'); hl = legend('Zhou','Camel','location','best');
plot(nanmean(p.efreq(:,ind)'),20*nanmean(zhou_emiss_100(:,ind)'-p.emis(:,ind)'),'b',nanmean(p.efreq(:,ind)'),namean(camel.trendall(:,ind)'),'r'); plotaxis2; hl = legend('20*emss diff','BT diff','location','best'); title('India : Zhou - Camel')
plot(nanmean(p.efreq(:,ind)'),nanmean(zhou_emiss_100(:,ind)'-p.emis(:,ind)'),'b',nanmean(p.efreq(:,ind)'),20*nanmean(camel.trendall(:,ind)'),'r'); plotaxis2; hl = legend('emss diff','20 * BT trend','location','best'); title('India : Zhou - Camel')

figure(6); clf
plot(nanmean(p.efreq(:,ind)'),20*nanmean(zhou_emiss_100(:,ind)'-p.emis(:,ind)'),'b',h.vchan,nanmean(tZhou(:,ind)'-tCamel(:,ind)'),'r'); plotaxis2; hl = legend('20*emss diff','BT diff','location','best'); title('India : Zhou - Camel')
plot(nanmean(p.efreq(:,ind)'),01*nanmean(zhou_emiss_100(:,ind)'-p.emis(:,ind)'),'b',h.vchan,nanmean(jacSTZhou(:,ind)'-jacSTCamel(:,ind)'),'r'); plotaxis2; hl = legend('emss diff','SKT jac diff','location','best'); title('India : Zhou - Camel')
