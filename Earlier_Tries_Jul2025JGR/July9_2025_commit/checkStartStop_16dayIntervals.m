load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/timestepsStartEnd_2002_09_to_2024_09.mat

iStep = input('Enter timestep (1:504) .. if you give two, then it prints out the list ... ');
if length(iStep) == 1
  [thedateS(iStep,:); thedateE(iStep,:)]
else
  [thedateS(iStep(1):iStemp(2),:) thedateE(iStep():iStep(2),:)]
end

