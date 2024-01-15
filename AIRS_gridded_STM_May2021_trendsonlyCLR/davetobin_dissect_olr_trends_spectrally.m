disp('see "dissect_olr_trends_spectrally.m" ')

addpath /asl/matlib/h4tools/
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP

if ~exist('olr0_bandsr')
  fUMBC = '/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_20_defaultGULP_STM_Oct2023_but_inconsistent.mat'; %% as the name says, WHO KNOWS, supposedly done around time of SOunder STM
  fUMBC = '/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_20.mat';                                          %% default, but things get written this all the time .. so could be dangerous/inconsistent
  fUMBC = '/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_20_defaultGULP_STM_Oct2023.mat';                  %% this is more consistent around time of Sounder Oct 2023 STM
  fUMBC = '/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_20_SEQN_vers10.mat';                              %% trends paper

  junk = load(fUMBC,'umbc_spectral_olr','strUMBC','iNumYears');  %% this is in compare_OLR_trend.m

  fprintf(1,'UMBC : %2i year dataset = %s \n',junk.iNumYears,junk.strUMBC)
  spectral_olr = junk.umbc_spectral_olr;
  dataset_name = 'UMBC';
end

olr0_bands    = spectral_olr.olr0_ecRad.bands;
skt_bands     = spectral_olr.skt_ecRad.bands;
planck_bands  = spectral_olr.planck_ecRad.bands;
lapse_bands   = spectral_olr.lapse_ecRad.bands;
o3_bands      = spectral_olr.o3_ecRad.bands;
wv_bands      = spectral_olr.wv_ecRad.bands;
t_co2_bands   = spectral_olr.ptemp_co2_ecRad.bands;               %% compute_feedbacks_generic_ecRad.m, line 352 : T + CO2 only
t_bands       = spectral_olr.perts9999.atmT_only_ecRad.bands;     %% compute_feedbacks_generic_ecRad.m, line 352 : T only
ghg_bands     = spectral_olr.perts9999.ghg_only_ecRad.bands;      %% compute_feedbacks_generic_ecRad.m, line 352 : T only
allpert_bands = spectral_olr.perts9999.atm_skt_ghg_ecRad.bands;   %% compute_feedbacks_generic_ecRad.m, line 156 : the whole shebang : SKT,T,g1,g2,g3,g4,g6

bands = [10 250 500 630 700 820 980 1080 1180 1390 1480 1800 2080 2250 2380 2600 3000];
bandcenter = 0.5*(bands(1:end-1) + bands(2:end));
olr0_clr = olr0_bands.clr;

% save dave_tobin_4608_olr0.mat olr0_clr bandcenter bands
