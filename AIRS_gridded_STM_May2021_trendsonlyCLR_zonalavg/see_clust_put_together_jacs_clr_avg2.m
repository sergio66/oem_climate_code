%% see clust_put_together_jacs_clr.m

JOB = iLatBin;

thedir0 = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Oct2020_startSept2002_trendsonly/';
thedir0 = ['AVGJAC64/Latbin' num2str(JOB,'%02d')];

favgmat = [thedir0 '/avg_72bins_m_ts_jac_latbin' num2str(JOB,'%02d') '.mat'];
if ~exist(favgmat)

  [h,ha,p,pa] = rtpread('pavg64.op.rtp');
  miaow = load('sarta_chans_for_l1c.mat');
  ind2834to2645 = miaow.ichan;
  
  for lon = 1 : 72
    fjac  = [thedir0 '/m_ts_jac_latbin' num2str(JOB,'%02d') '_lonbin' num2str(lon,'%02d') '.mat'];
    ajac = load(fjac,'nlays');
    ianumlay(lon) = ajac.nlays;
  end
  
  kcarta.subjac.ppmv2 = 0;
  kcarta.subjac.ppmv4 = 0;
  kcarta.subjac.ppmv6 = 0;
  m_ts_jac_fast = zeros(2645,297);
  iaSum = zeros(1,297);
  for lon = 1 : 72
    fjac  = [thedir0 '/m_ts_jac_latbin' num2str(JOB,'%02d') '_lonbin' num2str(lon,'%02d') '.mat'];
    ajac = load(fjac);
  
    kcarta.subjac.ppmv2 = kcarta.subjac.ppmv2 + ajac.theppmv.ppmv2;
    kcarta.subjac.ppmv4 = kcarta.subjac.ppmv4 + ajac.theppmv.ppmv4;
    kcarta.subjac.ppmv6 = kcarta.subjac.ppmv6 + ajac.theppmv.ppmv6;
  
    [mmf,nnf] = size(ajac.m_ts_jac_fast);
    xboo(lon) = (nnf-6)/3;
    n97   = 97;
    indTG = 1:6;                                                   iaSum(indTG)  = iaSum(indTG) + 1;
    indWV = xboo(lon)*0 + (1:xboo(lon)) + 6;    indWVx = n97*0 + (1:xboo(lon)) + 6;  iaSum(indWVx) = iaSum(indWVx) + 1;
    indTZ = xboo(lon)*1 + (1:xboo(lon)) + 6;    indTZx = n97*1 + (1:xboo(lon)) + 6;  iaSum(indTZx) = iaSum(indTZx) + 1;
    indO3 = xboo(lon)*2 + (1:xboo(lon)) + 6;    indO3x = n97*2 + (1:xboo(lon)) + 6;  iaSum(indO3x) = iaSum(indO3x) + 1;
  
    m_ts_jac_fast(:,indTG)  = m_ts_jac_fast(:,indTG)  + ajac.m_ts_jac_fast(:,indTG);
    m_ts_jac_fast(:,indWVx) = m_ts_jac_fast(:,indWVx) + ajac.m_ts_jac_fast(:,indWV);  
    m_ts_jac_fast(:,indTZx) = m_ts_jac_fast(:,indTZx) + ajac.m_ts_jac_fast(:,indTZ);  
    m_ts_jac_fast(:,indO3x) = m_ts_jac_fast(:,indO3x) + ajac.m_ts_jac_fast(:,indO3);  
  end
  
  nlays = 97;
  freq2645 = ajac.freq2645;
  m_ts_jac_fast = m_ts_jac_fast/72;
  kcarta.subjac.ppmv2 = kcarta.subjac.ppmv2/72;
  kcarta.subjac.ppmv4 = kcarta.subjac.ppmv4/72;
  kcarta.subjac.ppmv6 = kcarta.subjac.ppmv6/72;

  comment = 'code in see_clust_put_together_jacs_clr_avg2.m';
  saver = ['save ' favgmat ' nlays freq2645 m_ts_jac_fast kcarta comment'];
  fprintf(1,'saving zonal avg jacs for latbin %2i %s \n',JOB,favgmat);
  eval(saver)
  
else
  fprintf(1,'reading in zonal avg jacs for latbin %2i %s \n',JOB,favgmat);
  loader = ['load ' favgmat];
  eval(loader);
end

