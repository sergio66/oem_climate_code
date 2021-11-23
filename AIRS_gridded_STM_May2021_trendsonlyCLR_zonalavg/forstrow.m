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

% remove CO2, N2O,CH4
dbt_desc_hot_noTG = dbt_desc_hot;
dbt_desc_hot_noTG = dbt_desc_hot_noTG - co2jac*2.2/400;
dbt_desc_hot_noTG = dbt_desc_hot_noTG - n2ojac*0.8/300;
dbt_desc_hot_noTG = dbt_desc_hot_noTG - ch4jac*5.0/1860;

plot(jacx.freq2645,dbt_desc_hot_noTG,jacx.freq2645,dbt_desc_hot)
%save forsergio_noTG.mat dbt_desc_hot_noTG
