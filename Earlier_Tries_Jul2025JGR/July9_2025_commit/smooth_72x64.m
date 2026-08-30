function y = smooth_72x64(y0,N)

%% assumes y0 = 4608 x 1
%% can do a 3x3 smooth or a 5x5 smooth or 2x2 or 4x4

if nargin == 1
 N = 3;
end

if N == 1
  y = y0;
  return
end

if N < 2 | N > 5
  N = 3;
end

%{
iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);

load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

imagesc(reshape(junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2_ecRad,72,64)); colorbar; colormap(usa2); caxis([-1 +1]*2)
imagesc(reshape(smooth_72x64(junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2_ecRad),72,64)); colorbar; colormap(usa2); caxis([-1 +1]*2)

boo0 = junk.umbc_spectral_olr.feedback_ecRad.planck_ecRad + junk.umbc_spectral_olr.feedback_ecRad.lapse_ecRad + junk.umbc_spectral_olr.feedback_ecRad.o3_ecRad + junk.umbc_spectral_olr.feedback_ecRad.wv_ecRad;

N=1; figure(1); clf; imagesc(reshape(smooth_72x64(boo0,N),72,64)'); colorbar; colormap(usa2); caxis([-1 +1]*4)
N=2; figure(2); clf; imagesc(reshape(smooth_72x64(boo0,N),72,64)'); colorbar; colormap(usa2); caxis([-1 +1]*4)
N=3; figure(2); clf; imagesc(reshape(smooth_72x64(boo0,N),72,64)'); colorbar; colormap(usa2); caxis([-1 +1]*4)
N=4; figure(2); clf; imagesc(reshape(smooth_72x64(boo0,N),72,64)'); colorbar; colormap(usa2); caxis([-1 +1]*4)
N=5; figure(2); clf; imagesc(reshape(smooth_72x64(boo0,N),72,64)'); colorbar; colormap(usa2); caxis([-1 +1]*4)

booX = boo0;
booX(abs(booX) > 6 | abs(junk.results(:,6)') < 1e-4) = NaN;
N=1; figure(3); clf; imagesc(reshape(smooth_72x64(booX,N),72,64)'); colorbar; colormap(usa2); caxis([-1 +1]*4)
N=2; figure(4); clf; imagesc(reshape(smooth_72x64(booX,N),72,64)'); colorbar; colormap(usa2); caxis([-1 +1]*4)
N=3; figure(4); clf; imagesc(reshape(smooth_72x64(booX,N),72,64)'); colorbar; colormap(usa2); caxis([-1 +1]*4)
N=4; figure(4); clf; imagesc(reshape(smooth_72x64(booX,N),72,64)'); colorbar; colormap(usa2); caxis([-1 +1]*4)
N=5; figure(4); clf; imagesc(reshape(smooth_72x64(booX,N),72,64)'); colorbar; colormap(usa2); caxis([-1 +1]*4)

boo = boo0;
boo(abs(boo) > 6 | abs(junk.results(:,6)') < 1e-4) = NaN;
for ii = 1 : 10
  boo = smooth_72x64(boo);
end

booX = boo0;
booX(abs(booX) > 6 | abs(junk.results(:,6)') < 1e-4) = NaN;
for ii = 1 : 5
  booX = smooth_72x64(booX,5);
end

figure(1); clf; aslmapSergio(rlat65,rlon73,reshape(boo0,72,64)',[-90 +90],[-180 +180]);                 colormap(usa2); caxis([-1 +1]*3); colorbar
figure(2); clf; aslmapSergio(rlat65,rlon73,reshape(smooth_72x64(boo0),72,64)',[-90 +90],[-180 +180]);   colormap(usa2); caxis([-1 +1]*3); colorbar
figure(3); clf; aslmapSergio(rlat65,rlon73,reshape(boo,72,64)',[-90 +90],[-180 +180]);                  colormap(usa2); caxis([-1 +1]*3); colorbar
figure(4); clf; aslmapSergio(rlat65,rlon73,reshape(smooth_72x64(boo0,5),72,64)',[-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*3); colorbar
figure(5); clf; aslmapSergio(rlat65,rlon73,reshape(booX,72,64)',[-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*3); colorbar
%}

y = nan(size(y0));

z0 = 1 : 4608; z0 = z0';
z0 = reshape(z0,72,64);

y0 = y0(:);
[mm,nn] = size(y0);
if mm == 1
  y0 = y0';
end

y0 = reshape(y0,72,64);

if N == 3 | N == 5
  %% eg if N = 3 can easily loop over lon from 2 to 71, loop over lat from 2 to 63 and put the square there
  iS = 1 + (N-1)/2; iE = 72 - (N-1)/2;
  jS = 1 + (N-1)/2; jE = 64 - (N-1)/2;
elseif N == 2 | N == 4
  %% eg if N = 2 can easily loop over lon from 1 to 71, loop over lat from 1 to 63 and put the square there
  iS = 1; iE = 72 - (N)/2;
  jS = 1; jE = 64 - (N)/2;
end

%% and these are the indices you easily start from
ii0 = iS; jj0 = jS;

if N == 2
  indBlock = [001 002; 073 074];                        %% this is to do a 2x2 square around any point
elseif N == 3
  indBlock = [001 002 003; 073 074 075; 145 146 147];   %% this is to do a 3x3 square around any point
elseif N == 4
  indBlock = [001 002 003 004; 073 074 075 076; 145 146 147 148; 217 218 219 220];   %% this is to do a 4x4 square around any point
elseif N == 5
  indBlock = [001 002 003 004 005; 073 074 075 076 077; 145 146 147 148 149; 217 218 219 220 221; 289 290 291 292 293];   %% this is to do a 5x5 square around any point
end

ind0 = (jj0-1)*72 + ii0;

for ii = iS : iE
  for jj = jS : jE
    indx = (jj-1)*72 + ii;
    ind = indx - ind0 + indBlock;
    ind = ind(ind <= 4608);
    %[mm,nn]  = size(ind)
    %[ii jj indx mm nn min(ind(:)) max(ind(:))]
    y(indx) = nanmean(y0(ind(:)));
  end
end


    
