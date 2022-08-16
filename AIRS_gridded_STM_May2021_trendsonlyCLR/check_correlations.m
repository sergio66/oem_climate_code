addpath /home/sergio/MATLABCODE/COLORMAP
for jj = 01 : 64
  for ii = 1 : 72
    istr = num2str(ii,'%02d');
    jstr = num2str(jj,'%02d');
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/'];
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/LatBin' jstr '/LonBin' istr '/fits_LonBin' istr '_LatBin' jstr ];
      %fname = [fname '_V1_200300010001_201200120031_TimeStepsX228.mat'];
      %fname = [fname '_V1_200500010001_201400120031_TimeStepsX228.mat'];      
      fname = [fname '_V1_TimeSteps433.mat'];
    x = load(fname,'dbt_desc');
    dbt(:,ii) = x.dbt_desc(:,16);
  end
  correlations = corr(dbt',dbt');
  imagesc(correlations); colorbar; caxis([-1 +1]); title(num2str(jj)); colormap(usa2);
  pause(0.1)
    
  correlationsJJ(jj,:,:) = correlations;

end

comment = 'see check_correlations.m';
save -v7.3 /asl/s1/sergio/JUNK/latitude_correlations correlationsJJ comment
