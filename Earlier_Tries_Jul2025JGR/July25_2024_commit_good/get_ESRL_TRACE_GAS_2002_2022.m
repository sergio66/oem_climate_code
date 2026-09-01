function   co2x = get_ESRL_TRACE_GAS_2002_2022(JOBJOBJOB,iNumYears);

%{
plot(2002 + (1:20),get_ESRL_TRACE_GAS_2002_2022(2300,+(1:20)),'bx-');         title('trends')    %% 2300 is roughly half of the 608 tiles
hold
plot(2002 + (1:20),fliplr(get_ESRL_TRACE_GAS_2002_2022(2300,-(1:20))),'ro-'); title('trends')    %% 2300 is roughly half of the 608 tiles
hold off
hl = legend('forwards from 2002 t0 20XY','backwards from 2022 to 20XY','location','best','fontsize',10); grid
line([2000 2025],[2.22 2.22],'color','k')
title('CO2 trends; \newline notice blue end = red begin')
%}

if iNumYears > 0
  %% going forwards since 2002/09
  boo05 = load('/home/sergio/MATLABCODE_Git/ESRL_TRACE_GAS/carbon_tracker_500mb_2002_09_2007_08.mat');
  boo10 = load('/home/sergio/MATLABCODE_Git/ESRL_TRACE_GAS/carbon_tracker_500mb_2002_09_2012_08.mat');
  boo15 = load('/home/sergio/MATLABCODE_Git/ESRL_TRACE_GAS/carbon_tracker_500mb_2002_09_2017_08.mat');
  boo20 = load('/home/sergio/MATLABCODE_Git/ESRL_TRACE_GAS/carbon_tracker_500mb_2002_09_2022_08.mat');
else
  %% going backwards from  2022/08
  boo05 = load('/home/sergio/MATLABCODE_Git/ESRL_TRACE_GAS/carbon_tracker_500mb_2017_09_2022_08.mat');
  boo10 = load('/home/sergio/MATLABCODE_Git/ESRL_TRACE_GAS/carbon_tracker_500mb_2012_09_2022_08.mat');
  boo15 = load('/home/sergio/MATLABCODE_Git/ESRL_TRACE_GAS/carbon_tracker_500mb_2007_09_2022_08.mat');
  boo20 = load('/home/sergio/MATLABCODE_Git/ESRL_TRACE_GAS/carbon_tracker_500mb_2002_09_2022_08.mat');
end

co2x_05 = boo05.trend(JOBJOBJOB);
co2x_10 = boo10.trend(JOBJOBJOB);
co2x_15 = boo15.trend(JOBJOBJOB);
co2x_20 = boo20.trend(JOBJOBJOB);

if iNumYears > 0
  co2x = interp1([05 10 15 20],[co2x_05 co2x_10 co2x_15 co2x_20],iNumYears,[],'extrap');
else
  co2x = interp1(-[05 10 15 20],[co2x_05 co2x_10 co2x_15 co2x_20],iNumYears,[],'extrap');
end

if length(co2x) == 1
  fprintf(1,'set_CO2_CH4_N2O_ESRL.m : co2 --> co2 CarbonTracker, iNumYears = %2i\n',iNumYears)
end
