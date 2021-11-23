addpath /home/sergio/MATLABCODE/TIME

for yy = 2003 : 2017
  for ii = 1 : length(latbins)-1
    indx = (yy-2003+1)*40 + ii;
    pbefore.rlat(indx) = (latbins(ii) + latbins(ii+1))*0.5;
    pbefore.rlon(indx) = 0.0;
    pbefore.rtime(indx) = utc2taiSergio(yy,12,12,12.00);

    pafter.rlat(indx) = (latbins(ii) + latbins(ii+1))*0.5;
    pafter.rlon(indx) = 0.0;
    pafter.rtime(indx) = utc2taiSergio(yy,12,12,12.00);
  end
end
