function T = reader(fin,iTorR);

% This program reads binary data for 365/366 days and writes in ascii file.

% IMD High resolution 1By1 degree gridded daily temperature data
% (1951-2018)*. This data is arranged in 31x31 grid points. Lat 7.5N,
% 8.5N ... 36.5, 37.5 (31 Values). Long 67.5E, 68.5E ... 96.5, 97.5 (31
% Values). For leap years, data for 366 days are included. The unit of
% tempereture is in Celcious.

% Rainfall:
% 
% IMD New High Spatial Resolution (0.25X0.25 degree) Long Period
% (1901-2022) Daily Gridded Rainfall Data Set Over India. This data
% product is a very high spatial resolution daily gridded rainfall data
% (0.25 x 0.25 degree). The unit of rainfall is in millimeter (mm). Data
% available for 122 years, 1901 to 2022. Data is arranged in 135x129
% grid points. The first data in the record is at 6.5N & 66.5E, the
% second is at 6.5N & 66.75E, and so on. The last data record
% corresponds to 38.5N & 100.0E. The yearly data file consists of
% 365/366 records corresponding to non leap/ leap years.

if nargin == 0
  fin = '/asl/s1/sergio/NSD_India_Met_Data/Maxtemp_MaxT_2022.GRD';
  iTorR = +1;
end

if abs(iTorR) == 1
  ISIZ = 31;
  JSIZ = 31;
  NDAY = 366;
else
  ISIZ = 135;
  JSIZ = 129;
  NDAY = 366;
end

T = nan(NDAY,JSIZ,ISIZ);

fid = fopen(fin);
flen1 = fread(fid, 1, 'integer*4');

%for IDAY = 1:NDAY
%  %fprintf(1,'IDAY = %4i \n',IDAY);
%  for I=1:ISIZ
%    boo = fread(fid,[JSIZ 1],'real*4');
%    T(IDAY,I,1:length(boo)) = boo;
%  end
%end

for IDAY = 1:NDAY
  %fprintf(1,'IDAY = %4i \n',IDAY);
  for J=1:JSIZ
    boo = fread(fid,[ISIZ 1],'real*4');
    T(IDAY,J,1:length(boo)) = boo;
  end
end
flen1 = fread(fid, 1, 'integer*4');

fclose(fid);

if abs(iTorR) == 1
  T(T > 60) = nan;
  pcolor(squeeze(nanmean(T,1))); colorbar; shading interp; caxis([20 40])
else
  T(T < 0) = nan;
  pcolor(squeeze(nanmean(T,1))); colorbar; shading interp;
end
