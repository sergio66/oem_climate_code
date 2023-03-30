function A = read_hilo(fdir,fin0,p0);

if nargin == 0
  fdir = 'gml.noaa.gov/aftp/data/ozwv/WaterVapor/Hilo_New/';
  fin0 = [];
  p0 = 50 : 50 : 1000;
elseif nargin == 1
  fin0 = [];
  p0 = 50 : 50 : 1000;
elseif nargin == 2
  p0 = 50 : 50 : 1000;
end

p0 = sort(p0,'desc');
if length(fin0) == 0
  fin0 = 'HIH_H2O_20191209.txt';
end
fin = [fdir fin0];

fid = fopen(fin,'r');
tline = fgetl(fid);
%fprintf(1,'%s \n',tline)
while strcmp(tline(1:9),'Number km') == 0
  tline = fgetl(fid);  
  %fprintf(1,'%s \n',tline)
end

sizeA = [13 inf];
formatSpec = '%d %f %f %f %f %f %f %f %f %f %f %f %f';
A = fscanf(fid,formatSpec,sizeA);
A = A';
good = find(A(:,3) > 0 & A(:,6) > 0 & A(:,8) > 0);
A = A(good,:);
fclose(fid);

figure(1); clf

% Level  Alt    Press    Temp   Theta  RH_FPH RH_RS  H2Omr      H2Osd     H2Omr_orig  TFP_FPH TFP_RS  O3mr
% Number km     hPa      degC   degK    %RH    %RH   ppmv       ppmv      ppmv         degC    degC   ppbv

moo = find(A(:,3) == min(A(:,3)));
A = A(1:moo,:);

moo = find(A(:,3) <= 50);
A = A(1:moo,:);

A0 = A;
clear A;

[mm,nn] = size(A0);

if mm < 30
  A = [];
else
  A(:,1) = (1 : length(p0))';
  A(:,3) = p0;
  for ii = 2
    A(:,2) = interp1(log(A0(:,3)),A0(:,2),log(A(:,3)),[],'extrap');
  end
  
  for ii = 4 : 13
    A(:,ii) = max(interp1(log(A0(:,3)),A0(:,ii),log(A(:,3)),[],'extrap'),0);
  end
   
  subplot(131); semilogy(A(:,4)+273.13,A(:,3)); set(gca,'ydir','reverse'); title('T(z)')
  subplot(132); semilogy(A(:,6),A(:,3)); set(gca,'ydir','reverse'); title('RH')
  subplot(133); semilogy(A(:,8),A(:,3)); set(gca,'ydir','reverse'); title('ppmv')
end
