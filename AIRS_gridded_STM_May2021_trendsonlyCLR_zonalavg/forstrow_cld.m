load sarta_chans_for_l1c.mat
for ii = 1 : 64
  kcrad = ['KCARTAJACS_CLD_tracegas/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  loader = ['load ' kcrad];
  eval(loader)
  rcld(ii,:) = rKc;

  kcjac = ['KCARTAJACS_CLD_tracegas/individual_prof_convolved_kcarta_airs_' num2str(ii) '_coljac.mat'];
  loader = ['load ' kcjac];
  eval(loader)
  cldjacKc(ii,:,:) = rKc; %% gids 2,3,4,5,6 .. TA,ST

  %%%%%%%%%%%%%%%%%%%%%%%%%
  kcrad = ['KCARTAJACS_CLR_tracegas/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  loader = ['load ' kcrad];
  eval(loader)
  rclr(ii,:) = rKc;

  kcjac = ['KCARTAJACS_CLR_tracegas/individual_prof_convolved_kcarta_airs_' num2str(ii) '_coljac.mat'];
  loader = ['load ' kcjac];
  eval(loader)
  clrjacKc(ii,:,:) = rKc; %% gids 2,4,5,6,51,52 .. TA,ST
end
subplot(221); plot(fKc(ichan),nanmean(cldjacKc(:,ichan,1),1),'r',fKc(ichan),nanmean(clrjacKc(:,ichan,1),1),'b'); xlim([640 1640]); title('CO2')
subplot(222); plot(fKc(ichan),nanmean(cldjacKc(:,ichan,3),1),'r',fKc(ichan),nanmean(clrjacKc(:,ichan,2),1),'b'); xlim([640 1640]); title('N2O')
subplot(223); plot(fKc(ichan),nanmean(cldjacKc(:,ichan,4),1),'r',fKc(ichan),nanmean(clrjacKc(:,ichan,4),1),'b'); xlim([640 1640]); title('CH4')
subplot(224); plot(fKc(ichan),nanmean(cldjacKc(:,ichan,6),1),'r',fKc(ichan),nanmean(clrjacKc(:,ichan,8),1),'b'); xlim([640 1640]); title('Stemp')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load forsergio.mat

co2jac = zeros(2645,1);
n2ojac = zeros(2645,1);
ch4jac = zeros(2645,1);

for ii = 1 : 64
  jacname = ['/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR_zonalavg/'];
  jacname = [jacname '/AVGJAC64/Latbin' num2str(ii,'%02d') '/avg_72bins_m_ts_jac_latbin' num2str(ii,'%02d') '.mat'];
  jacx = load(jacname);   % read in a jac file
   co2jac = co2jac + jacx.m_ts_jac_fast(:,1);
   n2ojac = n2ojac + jacx.m_ts_jac_fast(:,2);
   ch4jac = ch4jac + jacx.m_ts_jac_fast(:,3);
end

co2jac = co2jac/64;
n2ojac = n2ojac/64;
ch4jac = ch4jac/64;

plot(fKc(ichan),nanmean(cldjacKc(:,ichan,1),1),fKc(ichan),co2jac)
plot(fKc(ichan),nanmean(cldjacKc(:,ichan,3),1),fKc(ichan),n2ojac)
plot(fKc(ichan),nanmean(cldjacKc(:,ichan,4),1),fKc(ichan),ch4jac)

co2jackc = nanmean(cldjacKc(:,ichan,1),1)';
n2ojackc = nanmean(cldjacKc(:,ichan,3),1)';
ch4jackc = nanmean(cldjacKc(:,ichan,4),1)';

% remove CO2, N2O,CH4
dbt_desc_allsky_noTG = dbt_desc_allsky;
dbt_desc_allsky_noTG = dbt_desc_allsky_noTG - co2jackc*2.2/400;
dbt_desc_allsky_noTG = dbt_desc_allsky_noTG - n2ojackc*0.8/300;
dbt_desc_allsky_noTG = dbt_desc_allsky_noTG - ch4jackc*5.0/1860;

plot(jacx.freq2645,dbt_desc_allsky_noTG,jacx.freq2645,dbt_desc_allsky)
%save forsergio_allsky_noTG.mat dbt_desc_allsky_noTG
