fprintf(1,'loading in %s for iVersJac = %4i would be sloooooooow\n',fname,iVersJac);

loader = ['load ' fname];
eval(loader);

%% no need to add in extra column for CFC12 >>>>>>
m_ts_jac_slow = [];

%boo = squeeze(subjac.jacCO2z(:,:,iLonBin)); boo = sum(boo,1);  m_ts_jac_slow = [m_ts_jac_slow boo'/log(10)];
boo = squeeze(subjac.coljacCO2(:,iLonBin));                    m_ts_jac_slow = [m_ts_jac_slow boo/log(10)];
boo = subjac.coljacN2O(:,iLonBin);                             m_ts_jac_slow = [m_ts_jac_slow boo/log(10)];
boo = subjac.coljacCH4(:,iLonBin);                             m_ts_jac_slow = [m_ts_jac_slow boo/log(10)];
boo = subjac.coljacCFC11(:,iLonBin);                           m_ts_jac_slow = [m_ts_jac_slow boo/log(10)];
boo = subjac.coljacCFC12(:,iLonBin);                           m_ts_jac_slow = [m_ts_jac_slow boo/log(10)];
boo = subjac.jacST(:,iLonBin);                                 m_ts_jac_slow = [m_ts_jac_slow boo];

nlays = subjac.nlevs(iLonBin)-1
boo = squeeze(subjac.jacWV(1:nlays,:,iLonBin));                 m_ts_jac_slow = [m_ts_jac_slow boo'];
boo = squeeze(subjac.jacT(1:nlays,:,iLonBin));                  m_ts_jac_slow = [m_ts_jac_slow boo'];
boo = squeeze(subjac.jacO3(1:nlays,:,iLonBin));                 m_ts_jac_slow = [m_ts_jac_slow boo'];
