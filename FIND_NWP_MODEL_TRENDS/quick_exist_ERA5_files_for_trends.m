for ii = iaMin : iaMax
  % iForwardFrom2002_or_BackwardFrom2022 = +1;
  %if mod(ii,100) == 0 & iForwardFrom2002_or_BackwardFrom2022 > 0
  %  fprintf(1,'+ \n')
  %elseif mod(ii,10) == 0 & iForwardFrom2002_or_BackwardFrom2022 > 0
  %  fprintf(1,'x')
  %elseif iForwardFrom2002_or_BackwardFrom2022 > 0
  %  fprintf(1,'.')
  %end

  if iOLR < 0
    if iDorA == 1
      fin = ['/asl/s1/sergio/alldata/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == -1
      fin = ['/asl/s1/sergio/alldata/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == 10
      fin = ['/asl/s1/sergio/alldata/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC/randomptera5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == -10
      fin = ['/asl/s1/sergio/alldata/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC/randomptera5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    end
  else
    if iDorA == +1
      fin = ['/asl/s1/sergio/alldata/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == -1
      fin = ['/asl/s1/sergio/alldata/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC_WithOLR/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == +10
      fin = ['/asl/s1/sergio/alldata/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/randomptera5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == -10
      fin = ['/asl/s1/sergio/alldata/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC_WithOLR/randomptera5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    end
  end

  if exist(fin)
    fprintf(1,'%3i : %s exists \n',ii,fin)
    iaFound(ii) = 1;
    if iExistOnly_or_NumProfsAlso < 0
      moo = load(fin,'pnew_op');
      iaNumProf(ii) = length(moo.pnew_op.stemp);
      if length(moo.pnew_op.stemp) ~= 4608
        fprintf(1,'OH NO %s %4i has %4i profiles instead of 4608 \n',fin,ii,iaNumProf(ii))
      end
    end
  else
    iaFound(ii) = 0;
    iaNumProf(ii) = -1;
    fprintf(1,'%3i : %s DNE \n',ii,fin)
  end
end
