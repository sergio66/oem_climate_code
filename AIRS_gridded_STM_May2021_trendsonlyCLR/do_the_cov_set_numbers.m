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

elseif iCovSetNumber == 20.2
  %%% THIS IS ALSO PRETTY DARN GOOD FOR OBS and ERA CAL, though WV is a but tooooo tight %%% THIS IS ALSO PRETTY DARN GOOD FOR OBS and ERA CAL, though WV is a but tooooo tight %%% THIS IS ALSO PRETTY DARN GOOD FOR OBS and ERA CAL, though WV is a but tooooo tight 
  cov_set = [1.0  0.05*1       0.09*1        1/2        0.01/4            0.01/4            1/2      0.01            0.01                1/2        03*1E-2     08*1E-1  05*1E-2];  %% Nov 2022 -- 2002/09-2022/08, *** Dec 6, 2022 commit 
                                                                                                                                                                                    %% GREAT T(z),WV too overdamped, awesome biases/std dev     
  %%% THIS IS ALSO PRETTY DARN GOOD FOR OBS and ERA CAL, though WV is a but tooooo tight %%% THIS IS ALSO PRETTY DARN GOOD FOR OBS and ERA CAL, though WV is a but tooooo tight %%% THIS IS ALSO PRETTY DARN GOOD FOR OBS and ERA CAL, though WV is a but tooooo tight 

  cov_set = [1.0  0.05*3       0.05*3        1/2        0.15/50*3         0.15/50*3         1/2      0.15/50*3       0.15/50*3           1/2        20*1E-7     20*1E-7  20*1E-7];  %% copied from iCovSetNumber == 4.16!!!!!! should work here, eh???
  cov_set = [1.0  0.05*1       0.10*1        1/2        0.15/50*2         0.15/50*2         1/2      0.15/50*3       0.15/50*3           1/2        20*1E-5     20*1E-5  20*1E-5];  %% adjst iCovSetNumber == 4.16!!!!!! oscillates
  cov_set = [1.0  0.05*1       0.10*1        1/2        0.01/4            0.01/4            1/2      0.025           0.025               1/2        10*1E+1     10*1E+1  10*1E+1];  %% try, not bad!!! T(z) too constrained, WV(z) could be loosed
  cov_set = [1.0  0.05*1       0.10*1        1/2        0.01/4            0.01/4            1/2      0.025           0.025               1/2        10*1E+0     10*1E+1  10*1E+1];  %% try, good!!! T(z) and WV(z) could still be loosened
  cov_set = [1.0  0.05*1       0.10*1        1/2        0.01/4            0.01/4            1/2      0.025           0.025               1/2        05*1E+0     05*1E+1  10*1E+1];  %% try, good!!! T(z) and WV(z) could still be loosened
  cov_set = [1.0  0.05*1       0.10*1        1/2        0.01/4            0.01/4            1/2      0.025           0.025               1/2        01*1E+0     01*1E+1  01*1E+1];  %% try, good!!! T(z) good but WV(z) oscillates
  cov_set = [1.0  0.05*1       0.10*1        1/2        0.01/4            0.01/4            1/2      0.025           0.025               1/2        01*1E+0     05*1E+1  01*1E+1];  %% try, good!!! QUITE GOOD AT TROPICS/MIDLATs, bad at N/S. Pole!!!! SAVE THIS!!!
  cov_set = [1.0  0.05*1       0.10*1        1/2        0.01/4            0.01/4            1/2      0.025           0.025               1/2        02*1E-1     10*1E+1  01*1E+1];  %% 

  %% cov_setA = tropics/midlats; cov_setB  = poles
  cov_setA = [1.0  0.05*1       0.10*1        1/2        0.01/4            0.01/4            1/2      0.025           0.025               1/2        01*1E+0     05*1E+1  01*1E+1];  %% try, good!!! QUITE GOOD AT TROPICS/MIDLATs, bad at N/S. Pole!!!! SAVE THIS!!!
  cov_setA = [1.0  0.05*10      0.09*10       1/2        0.01/1            0.01/1            1/2      0.01            0.01                1/2        03*1E-2     10*1E+0  05*1E-1];  %% combine the above two for POLES and use for TROPICS/MIDLATS iaSequential = -1
  cov_setA = [1.0  0.05*1       0.09*1        1/2        0.01/1            0.01/1            1/2      0.01            0.01                1/2        03*1E+1     20*1E-1  05*1E-2];  %% use for TROPICS/MIDLATS iaSequential = [150 60 100 -1 150 60 -1]
  cov_setA = [1.0  0.05*10      0.09*10       1/2        0.01/1            0.01/1            1/2      0.01            0.01                1/2        05*1E0     10*1E+0  05*1E-1];  %% combine the above two for POLES and use for TROPICS/MIDLATS iaSequential = -1
  cov_setA = [1.0  0.05*0.5      0.09*0.5       1/2        0.01/2            0.01/2            1/2      0.01            0.01                1/2        09*1E2     9*1E+3  05*1E-1];  %% combine the above two for POLES and use for TROPICS/MIDLATS iaSequential = -1
  cov_setA = [1.0  0.05*0.5      0.09*0.5       1/8        0.01/2            0.01/2            1/8      0.01/10        0.01                1/8        04*1E3     75*1E+4  05*1E-1];  %% TROPICS/MIDLATS iaSequential = -1 iia_OorC_DataSet_Quantile = [+0 09 05]; 2/4/23 git commit

  cov_setB = [1.0  0.05*1       0.09*1        1/2        0.01/4            0.01/4            1/2      0.01            0.01                1/2        03*1E-2     08*1E-1  05*1E-2];  %% Nov 2022 -- 2002/09-2022/08, *** Dec 6, 2022 commit 
      cov_setB(11:13) = cov_setB(11:13) .* [1e2 1e4 5e1];  %% NOT bad at all, but could relax it a little ... used in the Dec 6, 2022 commit where for OBS, dT(z)/dt from obs are mostly pretty good, WV too overdamped, great bias/std dev ****************
  cov_setB = [1.0  0.05*10      0.09*10       1/2        0.01/4            0.01/4            1/2      0.01            0.01                1/2        03*1E-2      08*1E+3  05*1E-1];  %% combine the above two for POLES
  cov_setB = [1.0  0.05*10      0.09*10       1/2        0.01/1            0.01/1            1/2      0.01            0.01                1/2        03*1E-2      10*1E+0  05*1E-1];  %% combine the above two for POLES
  cov_setB = [1.0  0.05*10      0.09*10       1/2        0.01/1            0.01/1            1/2      0.01            0.01                1/2        03*1E-2      10*1E+0  05*1E-1];  %% combine the above two for POLES *** NOT BAD 

  cov_setB = [1.0  0.50*10      0.10*10       1/2        0.01/1            0.01/1            1/2      0.01            0.01                1/2        01*1E-2      10*1E+0  05*1E-2];  %% combine the above two for POLES, pretty nice but T too wiggly
  cov_setB = [1.0  0.50*10      0.10*06       1/2        0.01/1            0.01/1            1/2      0.01            0.01                1/2        07*1E-2      10*1E+0  01*1E-2];  %% combine the above two for POLES, WV good, T too wiggly
  cov_setB = [1.0  0.10         0.10          1/2        0.01/1            0.01/1            1/2      0.01            0.01                1/2        07*1E-2      10*1E+0  01*1E-2];  %% works nicely for iaSequential = [150 60 100 -1 150 60 -1]

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%               sigT_t    sigT_s                sigWV_t   sigWV_s             sigO3_t   sigO3_s
  %%         lc   ct.lev1  ct.lev2   ct_wide     cw.lev1  cw.lev2    cw_wide  coz.lev1  coz.lev2  coz_wide    alpha_T  alpha_w  alpha_oz
  fprintf(1,'  in build_cov_matrices.m, iLatX = %2i \n',iLatX)

  %% COPY THE SETTINGS FOR ia_OorC_DataSet_Quantile = [01 09 16] synthetic calc
  cov_set = [1.0  0.05*5        0.09*5          1/2       0.04              0.02              1/2      0.02            0.02                1/2        02*1E-0     05*1E+2  05*1E+1];  %% 2002/09-2022/08, * reproduces ERA5 20 year gophysical rates dataset9,Quant16 *
  cov_set0 = cov_set;

  iSimpleLife = +1;
  iSimpleLife = -1;
  if iSimpleLife == +1
    %%% WHEN LIFE WAS SIMPLE, FEB 09, 2023 commit
    cov_setA = cov_set0;
    cov_setA(11:13) = cov_set0(11:13) *1e3;               %% Feb 4,2022 synthetic calc commit, works great for tropics obs as well

    cov_setB = cov_set0;
    if driver.iLat >= 64-iLatX | driver.iLat <= iLatX
      cov_setB(11:13) = cov_set0(11:13) .* [1.0 0.1 1.0] *1e4;     %% trying to figure out polar
      cov_setB(11:13) = cov_set0(11:13) .* [1.0 0.5 1.0] *1e3;     %% trying to figure out polar --- ok, not great, tried 2/5/23
      cov_setB(11:13) = cov_set0(11:13) .* [1.0 1.0 1.0] *1e3;     %% trying to figure out polar

      cov_setB(11:13) = cov_set0(11:13) .* [1.0   1.0e-1 1.0] *5e2;  %% trying to figure out polar
      cov_setB(11:13) = cov_set0(11:13) .* [1.0e4 1.0e-0 1.0] *1e-1; %% trying to figure out polar
      cov_setB([4 7 10]) = 1/8;
      cov_setB([2]) = 1.0;
      cov_setB([3]) = 1.0;
    end

  elseif iSimpleLife == -1
    %% MAKING LIFE MORE N MORE COMPLICATED
    cov_setA = cov_set0;
    cov_setA(11:13) = cov_set0(11:13) *1e3;               %% Feb 4,2022 synthetic calc commit, works great for tropics obs as well  
    cov_setB = cov_set0;

    if topts.resetnorm2one == +1

      if driver.iLat >= 64-iLatX | driver.iLat <= iLatX
        cov_setB(11:13) = cov_set0(11:13) .* [1.0 0.1 1.0] *1e4;     %% trying to figure out polar
        cov_setB(11:13) = cov_set0(11:13) .* [1.0 0.5 1.0] *1e3;     %% trying to figure out polar --- ok, not great, tried 2/5/23
        cov_setB(11:13) = cov_set0(11:13) .* [1.0 1.0 1.0] *1e3;     %% trying to figure out polar
  
        cov_setB(11:13) = cov_set0(11:13) .* [1.0   1.0e-1 1.0]   *5e2;  %% trying to figure out polar
  %      cov_setB(11:13) = cov_set0(11:13) .* [1.0e4 1.0e-0 1.0]   *1e-1; %% trying to figure out polar %%% THIS WAS COMMITED 8 am 2/9/23
  %      cov_setB(11:13) = cov_set0(11:13) .* [1.0e4 1.0e+4 1.0e2] *1e-1; %% trying to figure out polar %% mid Feb 2023

        cov_setB([4 7 10]) = 1/8;
        cov_setB([2 3]) = 1.0;

  %      cov_setB(11:13) = cov_set0(11:13) .* [1.0e4 1.0e+5 1.0e2] *1e-1; %% trying to figure out polar %% Mar 2023
  %      cov_setB(11:13) = cov_set0(11:13) .* [1.0e4 1.0e+2 1.0e2] *1e-1; %% trying to figure out polar %% Mar 2023
%  %        cov_setB([5 6]) = cov_set0([5 6])/5;
      end
    else
      if driver.iLat >= 64-iLatX | driver.iLat <= iLatX
        cov_setB(11:13) = cov_set0(11:13) .* [1.0e4 1.0e+4 1.0e2]  *1e-1;  %% trying to figure out polar  Feb 2023 worse
        cov_setB(11:13) = cov_set0(11:13) .* [2.0e3 1.0e+0 1.0e+1] *5e-2;  %% trying to figure out polar  Feb 2023 better ***
        cov_setB(11:13) = cov_set0(11:13) .* [5.0e2 5.0e-1 1.0e+1] *5e-2;  %% trying to figure out polar  Feb 2023 better
        cov_setB(11:13) = cov_set0(11:13) .* [2.0e2 2.0e-1 5.0e+0] *5e-2;  %% trying to figure out polar  Feb 2023 better, used in /asl/s1/sergio/JUNK/test7_guessstartWV_Vers1_march11_2023.mat, git commit Sun Mar 12 09:51:40 2023, Sun Mar 12 10:01:58 2023
        cov_setB([4 7 10]) = 1/8;
        cov_setB([2]) = 1.0;
        cov_setB([3]) = 1.0;

      end
    end
  end
  
  [wgtA,wgtB] = find_wgtA_wgtB(driver,iLatX);
  cov_set = wgtA * cov_setA + wgtB * cov_setB;

  %%%%%%%%%%%%%%%%%%%%%%%%% TESTING TO GET OBS COV MATRICES SAME AS SYNTHETIC ERA5 COV MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%    
  %%%%%%%%%%%%%%%%%%%%%%%%% TESTING TO GET OBS COV MATRICES SAME AS SYNTHETIC ERA5 COV MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%    

  %{
  cov_set = cov_setA;

  cov_set = [1.0  0.05*5        0.09*5          1/2       0.04              0.02              1/2      0.02            0.02                1/2        02*1E-0     05*1E+2  05*1E+1];  %% 2002/09-2022/08, * reproduces ERA5 20 year gophysical rates dataset9,Quant16 *
  cov_set0 = cov_set;
     cov_set(11:13) = cov_set0(11:13) *1e0;            %% this loosy goosy setting actually made dWVfrac/dt > 0 at surface in all latitudes but has low S/N : about 2-3
     cov_set(11:13) = cov_set0(11:13) *1e3;            %% to duplicate settings of synthetic rate retrieval, Feb 4,2022 commit, and has higher S/N, about 6, but ERA5/UMBC dWVfrac/dt is about x2
     cov_set(11:13) = cov_set0(11:13) *1e1;            %% loosen things, not bad the ERA5/UMBC troposphere WVfrac ratio now 1.14 but things too wiggly at bottom
     
     %% this is not bad, but have not really improved dcolW/dt 
     cov_set = cov_set0;
     cov_set(11:13) = cov_set0(11:13) .* [1e2 1 1e3];  %% loosen tikonov params
        cov_set(2:3) = cov_set0(2:3)/5;                %% tighten cov
        cov_set(5:6) = cov_set0(5:6) .* [1/4 1/2];     %% tighten cov

     cov_set = cov_set0;
     cov_set(11:13) = cov_set0(11:13) .* [1e2 1 1e3]/10;  %% loosen tikonov params
        cov_set(2:3) = cov_set0(2:3)/1;                   %% leave cov alone
        cov_set(5:6) = cov_set0(5:6) .* [1 1];            %% leave cov alone

     %%% way toooooo loooosey goooosey
     cov_set = cov_set0;
     cov_set(11:13) = cov_set0(11:13) .* [1e2 1 1e3]/100; %% loosen tikonov params
        cov_set(2:3) = cov_set0(2:3)*10;                   %% loosen cov
        cov_set(5:6) = cov_set0(5:6)*10;                   %% loosen cov

     %%% way toooooo loooosey goooosey
     cov_set = cov_set0;
     cov_set(11:13) = cov_set0(11:13) .* [1e2 1 1e3]/20;   %% loosen tikonov params
        cov_set(2:3) = cov_set0(2:3)*2;                    %% loosen cov
        cov_set(5:6) = cov_set0(5:6)*2;                    %% loosen cov
     cov_set(11:13) = cov_set0(11:13) .* [1e2 1 1e3]/50;   %% loosen tikonov params
        cov_set(2:3) = cov_set0(2:3)*5;                    %% loosen cov
        cov_set(5:6) = cov_set0(5:6)*5;                    %% loosen cov
     cov_set(11:13) = cov_set0(11:13) .* [1e2 1 1e3]/50;   %% loosen tikonov params
        cov_set(2:3) = cov_set0(2:3)*2.5;                  %% loosen cov
        cov_set(5:6) = cov_set0(5:6)*2.5;                  %% loosen cov
  %}

  %%%%%%%%%%%%%%%%%%%%%%%%% TESTING TO GET OBS COV MATRICES SAME AS SYNTHETIC ERA5 COV MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%    
  %%%%%%%%%%%%%%%%%%%%%%%%% TESTING TO GET OBS COV MATRICES SAME AS SYNTHETIC ERA5 COV MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%    

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

