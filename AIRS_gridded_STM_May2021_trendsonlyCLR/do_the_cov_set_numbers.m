%%            1  |      2        3       4     |       5         6           7  |     8         9        10     |    11       12        13   
%%               |   sigT_t    sigT_s          |    sigWV_t   sigWV_s           |   sigO3_t   sigO3_s           |
%%            lc |  ct.lev1  ct.lev2   ct_wide |    cw.lev1  cw.lev2    cw_wide |  coz.lev1  coz.lev2  coz_wide |   alpha_T  alpha_w  alpha_oz

if iCovSetNumber == 4.16
  cov_set = [1.0  0.05*1        0.05*3          1/2       0.15/50*1         0.15/50*3         1/2      0.15/50*1       0.15/50*3           1/2        20*1E-7     20*1E-7  20*1E-7];  %% try strat only x3   unc
  cov_set = [1.0  0.05*1        0.05*3          1/2       0.15/50*1         0.15/50*3         1/2      0.15/50*1       0.15/50*3           1/2        20*1E+2     20*1E+2  20*1E+2];  %% try strat only x3   unc

  cov_set = [1.0  0.025         0.05            1/2       0.15/50           0.15/50           1/2      0.15/50         0.15/50             1/2        20*1E-7     20*1E-7  20*1E-7];  %% ok   excellent simulated ERA5 spectral rates Feb 4, 2022
  cov_set = [1.0  0.05*1        0.05*1          1/2       0.15/50*1         0.15/50*1         1/2      0.15/50*1       0.15/50*1           1/2        20*1E-7     20*1E-7  20*1E-7];  %% very excellent simulated ERA5 spectral rates Feb 4, 2022 till Feb 15, 2022

  cov_set = [1.0  0.05*3        0.05*3          1/2       0.02              0.02              1/2      0.02            0.02                1/2        20*1E-7     20*1E-7  20*1E-7];  %% try x100 unc, Feb 16 2022-Apr7,2022 : great JPLMay 2022 talk!!! huge retr uncertainty
  cov_set = [1.0  0.05*3        0.05*3          1/2       0.15/50*3         0.15/50*3         1/2      0.15/50*3       0.15/50*3           1/2        20*1E-7     20*1E-7  20*1E-7];  %% try x3   unc, Feb 16 2022-Apr7,2022 : great but maybe still constricts WV/O3
  cov_set = [1.0  0.05*1        0.05*3          1/2       0.15/50*1         0.15/50*3         1/2      0.15/50*1       0.15/50*3           1/2        20*1E-7     20*1E-7  20*1E-7];  %% try strat only x3   unc  04/23/2022 commit 30d2e554a97b34b0923ad58346d183a3c10d6bcb

elseif iCovSetNumber == 12
  cov_set = [1.0  0.05*3        0.05*3          1/2       0.02              0.02              1/2      0.02            0.02                1/2        20*1E-7     20*1E-7  20*1E-7];  %% 12 year rates, init try 2002/09-2014/08
  cov_set = [1.0  0.05*3        0.05*3          1/2       0.02              0.02              1/2      0.02            0.02                1/2        20*1E-4     20*1E-4  20*1E-4];  %% 2002/09-2014/08, * used this for Princeton iQuant=50, 
                                                                                                                                                                                      %% and GOOD ERA5 retr dataset4,Quant16 **
  cov_set = [1.0  0.05*1        0.05*1          1/2       0.02/5            0.02/5            1/2      0.02/5          0.02/5              1/2        20*1E-4     20*1E-4  20*1E-4];  %% 2002/09-2014/08 12 years AMIP6/CMIP6 for Princeton
  cov_set = [1.0  0.05*1        0.05*1          1/2       0.08/20           0.08/20           1/2      0.08/20         0.08/20             1/2        20*1E-4     20*1E-4  20*1E-4];  %% 2002/09-2014/08 12 years AMIP6/CMIP6 for Princeton

elseif iCovSetNumber == 18 | iCovSetNumber == 19
  cov_set = [1.0  0.05*3        0.05*3          1/2       0.02              0.02              1/2      0.02            0.02                1/2        05*1E-4     05*1E-4  05*1E-4];  %% 2002/09-2020/08, * reproduces ERA5 20 year gophysical rates dataset9,Quant16 **
  cov_set = [1.0  0.05*5        0.05*5          1/2       0.02              0.02              1/2      0.02            0.02                1/2        01*1E-4     05*1E-4  05*1E-4];  %% 2002/09-2020/08, * reproduces ERA5 20 year gophysical rates dataset9,Quant16 **
  cov_set = [1.0  0.05*5        0.09*5          1/2       0.04              0.02              1/2      0.02            0.02                1/2        01*1E-4     05*1E-4  05*1E-4];  %% 2002/09-2020/08, * reproduces ERA5 20 year gophysical rates dataset9,Quant16 **

elseif iCovSetNumber == 20.0
  cov_set = [1.0  0.05*1       0.09*1        1/2        0.01              0.01              1/2      0.01            0.01                1/2        01*1E-3     05*1E-3  05*1E-2];  %% Nov 2022 -- 2002/09-2022/08, * too loosy goosy * 
  cov_set = [1.0  0.05*1       0.09*1        1/2        0.01              0.01              1/2      0.01            0.01                1/2        05*1E-2     05*1E-2  05*1E-2];  %% Nov 2022 -- 2002/09-2022/08, * pretty good but might be little too tight *
  cov_set = [1.0  0.05*1       0.09*1        1/2        0.01              0.01              1/2      0.01            0.01                1/2        01*1E-2     05*1E-2  05*1E-2];  %% Nov 2022 -- 2002/09-2022/08, * pretty good but might be little too loose *

  %%%%% March 2023  
  cov_set = [1.0  0.05*5        0.09*5          1/2       0.04              0.02              1/2      0.02            0.02                1/2        02*1E-0     05*1E+2  05*1E+1];  %% 2002/09-2022/08, * reproduces ERA5 20 year gophysical rates dataset9,Quant16 *
  cov_set0 = cov_set;
  cov_set = cov_set0;
  cov_set(11:13) = cov_set0(11:13) .* [1e2 1 1e3]/20;   %% loosen tikonov params
    cov_set(2:3) = cov_set0(2:3)*2;                     %% loosen cov
    cov_set(5:6) = cov_set0(5:6)*2;                     %% loosen cov

elseif iCovSetNumber == 20.1
  %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%%
  cov_set = [1.0  0.05*1       0.09*1        1/2        0.01              0.01              1/2      0.01            0.01                1/2        05*1E-2     02*1E-2  05*1E-2];  %% Nov 2022 -- 2002/09-2022/08, *** bloody good ocb_set=1 testing, Q=16,dataset=9 
                                                                                                                                                                                    %% NP is iffy BUT obs T may be too tight, obs WV maybe too loose***
                                                                                                                                                                                    %% see eg Output_CAL/Quantile16_20years and Output_CAL/Quantile16_20yearsV2
  %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%%
  if topts.ocb_set == 1
    %% this is the Nov 16, 2022 commit commit 41a6282ca1902035330b78d4e378c6d9aba23491
    cov_set = [1.0  0.05*5        0.09*5          1/2       0.04              0.02              1/2      0.02            0.02                1/2        01*1E-4     05*1E-4  05*1E-4];  %% 2002/09-2020/08, * reproduces ERA5 20 year gophysical rates dataset9,Quant16 *
    %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%%
    cov_set = [1.0  0.05*5        0.09*5          1/2       0.04              0.02              1/2      0.02            0.02                1/2        02*1E-0     05*1E+2  05*1E+1];  %% 2002/09-2022/08, * reproduces ERA5 20 year gophysical rates dataset9,Quant16 *
    cov_set(11:13) = cov_set(11:13) *1e3;  %% Feb 4,2022 commit
    %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%% %%%% YES YES YES for ERA5 cal %%%%%
  end

elseif iCovSetNumber == 20.2 | iCovSetNumber == 23
  do_covset_20p2
   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif iCovSetNumber == 0
  %% testing and trying

  %  cov_set = [1.0  0.05*1        0.09*1        1/2        0.01/4               0.01/4              1/2      0.01              0.01                  1/2        03*1E-2     09*1E-2  05*1E-2];  %% Nov 2022 -- 2002/09-2022/08, *** WV too flexible?
  %  cov_set = [1.0  0.05*1        0.09*1        1/2        0.01/4               0.01/4              1/2      0.01              0.01                  1/2        03*1E-2     10*1E-2  05*1E-2];  %% Nov 2022 -- 2002/09-2022/08, *** WV too flexible?
  %  cov_set = [1.0  0.05*1/3      0.05*1/3        1/3       0.02/5/3            0.02/5/3            1/3      0.02/5/3          0.02/5/3              1/2        20*1E-4     20*1E-4  20*1E-4];  %% 20 years, iQAX=3
  %  cov_set = [1.0  0.05*1/2      0.05*1/2        1/2       0.02/5/2            0.02/5/2            1/2      0.02/5/2          0.02/5/2              1/2        20*1E-4     20*1E-4  20*1E-4];  %% 20 years, iQAX=3
  %  cov_set = [1.0  0.05*3        0.05*3          1/2       0.02/5/2            0.02/5/2            1/2      0.02/5/2          0.02/5/2              1/2        20*1E-4     20*1E-4  20*1E-4];  %% 20 years trying SHTUFF SHTUFF, T rates good, WV too constrained
  %    cov_set(12:13) = cov_set(12:13) *1e2;
  %    cov_set(11:11) = cov_set(11:11) *1e1;

  cov_set = [1.0  0.05*1       0.09*1        1/2        0.01              0.01              1/2      0.01            0.01                1/2        03*1E-2     08*1E-2  05*1E-2];  %% Nov 2022 -- 2002/09-2022/08, *** TRY THIS, WV still too loose
  cov_set = [1.0  0.05*1       0.09*1        1/2        0.01/8            0.01/8            1/2      0.01            0.01                1/2        03*1E-2     08*1E-2  05*1E-2];  %% Nov 2022 -- 2002/09-2022/08, *** very bad bad WV
  cov_set = [1.0  0.05*1       0.09*1        1/2        0.01/4            0.01/4            1/2      0.01            0.01                1/2        03*1E-2     08*1E-1  05*1E-2];  %% Nov 2022 -- 2002/09-2022/08, *** WV too tight, commmitted Dec 6 
  cov_set = [1.0  0.05*1       0.09*1        1/2        0.01/4            0.01/4            1/2      0.01            0.01                1/2        03*1E-2     08*1E-2  05*1E-2];  %% Dec 2022 -- 2002/09-2022/08, *** WV too flexible, bias great, std dev little large

  cov_set = [1.0  0.05*1       0.09*1        1/2        0.01              0.01              1/2      0.01            0.01                1/2        05*1E-2     02*1E-2  05*1E-2];  %% Nov 2022 -- 2002/09-2022/08, *** bloody good ocb_set=1 testing, Q=16,dataset=9 
  cov_set = [1.0  0.05*1       0.05*1        1/2        0.01/4            0.01/4            1/2      0.01            0.01                1/2        03*1E-2     10*1E-2  05*1E-2]; 

  %  cov_set = [1.0  0.05*1       0.09*1        1/2        0.01/5            0.01/5            1/2      0.01            0.01                1/2        03*1E-2     08*1E-1  05*1E-2];  %% try this, bias good std dev sucks, WV too tight
  %  cov_set = [1.0  0.05*1       0.09*1        1/2        0.01/3            0.01/3            1/2      0.01            0.01                1/2        03*1E-2     08*1E-1  05*1E-2];  %% try this, bias good std dev sucks, WV too tight

else
  iCovSetNumber
  error('unknown iCovSetNumber')
end

%%%%%%%%%% <<<<<<<<<< %%%%%%%%%% >>>>>>>>>> %%%%%%%%%% <<<<<<<<<< %%%%%%%%%% >>>>>>>>>> %%%%%%%%%%
%%%%%%%%%% <<<<<<<<<< %%%%%%%%%% >>>>>>>>>> %%%%%%%%%% <<<<<<<<<< %%%%%%%%%% >>>>>>>>>> %%%%%%%%%%

%{
%%% this complicates things tooooooo much to remember
cov_setX = cov_set; 
if (topts.iChSet == 4 | topts.iChSet == 5) & topts.dataset >= 8
  if iCovSetNumber == 20.0 | iCovSetNumber == 20.1
    cov_set(11:13) = cov_setX(11:13) *1e2;  %% default but I think a little tooooo loosey goosey till Nov 2022
    cov_set(11:13) = cov_setX(11:13) .* [1e0 1e0 1e0];  %% THIS EXTRA MULT WAS TOO COMPLICATED! and it was because I screwed up mat_od with scale lenghts below! so just set to 1
    %%   cov_set(11:13) = cov_setX(11:13) .* [1e2 5e5 1e2];       %% NOT bad at all, but could relax it a little
    %%   cov_set(11:13) = cov_setX(11:13) .* [1e2 5e5 1e2] * 1/2;
  elseif iCovSetNumber == 20.2
    %%% THIS IS ALSO PRETTY DARN GOOD FOR OBS and ERA CAL, though WV is a but tooooo tight %%% THIS IS ALSO PRETTY DARN GOOD FOR OBS and ERA CAL, though WV is a but tooooo tight %%% THIS IS ALSO PRETTY DARN GOOD FOR OBS and ERA CAL, though WV is a but tooooo tight 
    cov_set(11:13) = cov_setX(11:13) .* [1e2 1e4 5e1];  %% NOT bad at all, but could relax it a little ... used in the Dec 6, 2022 commit where for OBS, dT(z)/dt from obs are mostly pretty good, WV too overdamped, great bias/std dev ****************
    cov_set(11:13) = cov_setX(11:13) .* [1e0 1e0 1e0];  %% THIS EXTRA MULT WAS TOO COMPLICATED! and it was because I screwed up mat_od with scale lenghts below! so just set to 1
    %%% THIS IS ALSO PRETTY DARN GOOD FOR OBS and ERA CAL, though WV is a but tooooo tight %%% THIS IS ALSO PRETTY DARN GOOD FOR OBS and ERA CAL, though WV is a but tooooo tight %%% THIS IS ALSO PRETTY DARN GOOD FOR OBS and ERA CAL, though WV is a but tooooo tight 
  end
endd
%}

%%%%%%%%%% <<<<<<<<<< %%%%%%%%%% >>>>>>>>>> %%%%%%%%%% <<<<<<<<<< %%%%%%%%%% >>>>>>>>>> %%%%%%%%%%
%%%%%%%%%% <<<<<<<<<< %%%%%%%%%% >>>>>>>>>> %%%%%%%%%% <<<<<<<<<< %%%%%%%%%% >>>>>>>>>> %%%%%%%%%% 
%%   if sig_q -> 0   then you say you are VERY sure about a-priori ==> do not change ==> delta(param) --> 0
%%      sig_q -> INF then you say you are DO NOT TRUST    a-priori ==>        change ==> delta(param) --> bigly wigly
%%   if alpha -> 0   then you say you are DO NOT TRUST    a-priori ==>        change ==> delta(param) --> bigly wigly
%%      alpha -> INF then you say you are VERY sure about a-priori ==> do not change ==> delta(param) --> 0

