str = 'ERA5';
str = 'CMIP6';
str = 'AIRSL3';

for ii = 1 : 64
  ii

  f1 = ['simulate64bins' str     num2str(ii) '.ip.rtp'];
  f2 = ['simulate64bins' str '_' num2str(ii) '.ip.rtp'];
  mver = ['!/bin/mv ' f1 ' ' f2];
  eval(mver);

  f1 = ['simulate64bins' str     num2str(ii) '.op.rtp'];
  f2 = ['simulate64bins' str '_' num2str(ii) '.op.rtp'];
  mver = ['!/bin/mv ' f1 ' ' f2];
  eval(mver);

  f1 = ['simulate64bins' str     num2str(ii) '.rp.rtp'];
  f2 = ['simulate64bins' str '_' num2str(ii) '.rp.rtp'];
  mver = ['!/bin/mv ' f1 ' ' f2];
  eval(mver);
end
