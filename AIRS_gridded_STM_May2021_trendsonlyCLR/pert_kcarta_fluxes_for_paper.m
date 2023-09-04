%% from driver_compute_feedbacks_from_smallFATfile.m

addpath /asl/matlib/h4tools

fprintf(1,'iNumYears = %2i .. reading in %s \n',iNumYears,strUMBC);
loader = ['load ' strUMBC];
eval(loader);

read_fileMean17years
h = hMean17years;
p = pMean17years;

get_nan_bottom_layer
p = make_rtp_plays(p);
p.plays(p.plays <= eps) = NaN;
plays = nanmean(p.plays,2);
plays = plays(1:100);

resultsTunc  = 0 * resultsT;
resultsWVunc = 0 * resultsT;
resultsO3unc = 0 * resultsT;
resultsunc   = 0 * results;

iMakeUnc = +1;
if ~exist('fracO3') | iMakeUnc > 0
  disp('computing fracO3')
  if ~exist('maskLF')
    disp('  warning : setting maskLF = ones(1,4608)')
    maskLF = ones(1,4608);
    maskLFmatr = reshape(maskLF,72,64)';
  end
  if ~exist('iNumLay')
    disp('  warning : setting iNumLay = 49')
    iNumLay = 49; 
  end
  if ~exist('xb')
    disp('  warning : setting xb = 0')
    xb = zeros(iNumLay,4608);
  end
  if ~exist('pavg')
    disp('  warning : setting pavg = pjunk20')
    pavg = pjunk20;
  end
  if ~exist('resultsunc')
    disp('  warning : setting resultsunc = results/100')
    resultsunc = results/100;
  end
  interp_resultsT_WV_O3_to_p
  get_deltaO3
end
