function ind = get_climateQA_LW_noCFC(f)

%{
%clust_run_retrieval_setlatbin_AIRS_loop_lonbin.m
  iChSet = topts.iChSet;
  iChSet = 2; %% new chans
  iChSet = 4; %% new chans + Tonga (high alt) + more LW
  iChSet = 1; %% old chans (default)
  iChSet = 3; %% new chans, but no CFC11   STROW PRISTINE SET, AMT 2019, also used for JPL April 2021 Sounder Science Meeting
  %iChSet = 5; %% new chans, + CO2 laser lines (window region, low altitude T)
  %iChSet = 6; %% SW T(z) chans + 800 - 1600 cm- 1 lines (window region, low altitude T)
  if topts.dataset == 4
    iChSet = 3; %% new chans, but no CFC11   STROW PRISTINE SET, AMT 2019, also used for JPL April 2021 Sounder Science Meeting
  end
  topts.iChSet = iChSet;

find_the_oem_channels.m
%}


iNum = length(f);

if iNum == 2378
  ch = choose_goodchans_from_2378;                      %% Old 2378 chans
else
   %% change_important_topts_settings.m:41:topts.chan_LW_SW =  0;  %% just LW/MW DEFAULT, keep strat T chans                       470 chans
  settings_chan_LW_SW = 0;
  ch = choose_goodchans_from_2645(settings_chan_LW_SW); %% New 2645 chans
end

    iless700 = find(f(ch) < 700);
      iless700 = iless700(1:3:length(iless700));
    imore700 = find(f(ch) >= 700 & f(ch) < 730);
      imore700 = imore700(1:3:length(imore700));
    imore730 = find(f(ch) >= 730 & f(ch) < 800);
      imore730 = imore730(1:2:length(imore730));
    chx = ch;
    chx = chx([iless700; imore700; imore730]);
    chx = setdiff(chx,76);   %% 668.0348 cm-1 is a bad chan

    moo = load('climateQA_LW_noCFC.txt');
    ch = moo(:,2);
    %ch = union(ch,chx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp('see Mitchell_Goldberg-Dissertation.pdf for list of chans : wget https://aosc.umd.edu/sites/default/files/dissertations-theses/Mitchell%20Goldberg-Dissertation.pdf')
  disp('Fig 5.8        667.766 cm-1       1 mb')
  disp('Fig 5.9        667.715 cm-1       2 mb')
  disp('Fig 5.9        667.270 cm-1      15 mb')
  disp('Fig 5.8        667.018 cm-1      25 mb')
  disp('Fig 4.4        681.457 cm-1      90 mb')
  disp('Fig 4.5        704.436 cm-1     350 mb')
  disp('Fig 4.6        723.029 cm-1     700 mb')
  disp('Fig 4.7        801.099 cm-1     850 mb')
  disp(' ')
  disp('Fig 4.8        1519.07 cm-1     315 mb')
  disp('Fig 4.9        1598.49 cm-1     490 mb  MidTropWV')
  disp('Fig 5.10       1519.07 cm-1     315 mb  UTWV')
  disp('Fig 5.1        1520.87 cm-1     315 mb???')
  disp(' ')
  disp('Fig 5.4        1040.03 cm-1     80 mb???')
  disp(' ')

ii = [];
junk = find(f >= 0667.781,1); ii = [ii junk];
junk = find(f >= 0667.276,1); ii = [ii junk];
junk = find(f >= 0667.024,1); ii = [ii junk];
junk = find(f >= 0681.467,1); ii = [ii junk];
junk = find(f >= 0704.436,1); ii = [ii junk];
junk = find(f >= 0723.029,1); ii = [ii junk];
junk = find(f >= 0801.099,1); ii = [ii junk];
junk = find(f >= 1419.07,1); ii = [ii junk];
junk = find(f >= 1519.07,1); ii = [ii junk];
junk = find(f >= 1598.49,1); ii = [ii junk];
junk = find(f >= 1520.87,1); ii = [ii junk];
junk = find(f >= 1040.03,1); ii = [ii junk];

addpath /home/sergio/MATLABCODE/DUSTFLAG
freqsNchans
junk = find(f >= ff(1),1); ii = [ii junk];
junk = find(f >= ff(2),1); ii = [ii junk];
junk = find(f >= ff(3),1); ii = [ii junk];
junk = find(f >= ff(4),1); ii = [ii junk];
junk = find(f >= ff(5),1); ii = [ii junk];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ch = union(ch,ii);
ch = sort(ch);
fprintf(1,'found %3i chans \n',length(ch));

ind = ch;
