%Year Mo Globe  Land Ocean   NH   Land Ocean   SH   Land Ocean Trpcs  Land Ocean NoExt  Land Ocean SoExt  Land Ocean NoPol  Land Ocean SoPol  Land Ocean USA48
data = load('uah.txt');

dope = find(data(:,1) < 0);
trend = data(dope,[8 11 14])

dope = find(data(:,1) >= 0);
data = data(dope,:);


data_oceantrp = data(:,[1 2 14]);
data_oceanNH  = data(:,[1 2 8]);
data_oceanSH  = data(:,[1 2 11]);

lala = 1:length(dope);
figure(1); plot(1979 + (lala-1)/12,[data_oceanSH(:,3) data_oceantrp(:,3) data_oceanNH(:,3)]);
  grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = load('uah.txt');
dope = find(data(:,1) >= 2002);
data = data(dope,:);

data_oceantrp = data(:,[1 2 14]);
data_oceanNH  = data(:,[1 2 8]);
data_oceanSH  = data(:,[1 2 11]);

lala = 1:length(dope);
lala = lala(9:9+10*12 - 1); %% Sept 2002 - Aug 2012
data_oceantrp = data_oceantrp(lala,:);
data_oceanNH  = data_oceanNH(lala,:);
data_oceanSH  = data_oceanSH(lala,:);

figure(2); plot(2002 + 8/12 + (1:length(lala))/12,[data_oceanSH(:,3) data_oceantrp(:,3) data_oceanNH(:,3)]);
  grid