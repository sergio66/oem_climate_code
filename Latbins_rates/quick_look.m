plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
plevsA = plevs(1:end-1) - plevs(2:end);
plevsB = log(plevs(1:end-1)./plevs(2:end));
plevs = plevsA./plevsB;
plays = plevs(4:100); plays = flipud(plays);

N = input('Enter how many files to look at : ');
if N == 5
  save_lat = [-65 -45 0 +45 +65];
else
  load latbins.mat;
end

for ix = 1 : N
  loader = ['load ../Output/testx_' num2str(ix) '.mat'];
  eval(loader)
  params(ix,:)     = oem.finalrates;
  params_sigs(ix,:) = oem.finalsigs;
  fitted_rates(ix,:) = oem.fit;
  input_rates(ix,:)  = rateset.rates;
end

thestr = jacobian.qstnames;
chanset = jacobian.chanset;
chanset = jacobian.chanset_used;
lambdax = struct;
  lambdax.lambda      = oem.lambda;
  lambdax.lambda_qst  = oem.lambda_qst;
  lambdax.lambda_Q1   = oem.lambda_Q1;
  lambdax.lambda_temp = oem.lambda_temp;
clear oem rateset jacobian lls

if length(params) == 200
  water      = params(:,007:103);
  water_sigs = params_sigs(:,007:103);
  temp       = params(:,104:200);
  temp_sigs  = params_sigs(:,104:200);
  tracegases      = params(:,1:6);
  tracegases_sigs = params_sigs(:,1:6);
  %thestr = {'CO2','O3','N2O','CH4','CFC11','stemp'};
else
  error('oops')
end

figure(1); clf; plot(water,1:97); set(gca,'ydir','reverse'); title('WV(z)')
figure(2); clf; plot(temp,1:97);  set(gca,'ydir','reverse'); title('T(z)')

figure(1); clf; pcolor(save_lat,1:97,water'); set(gca,'ydir','reverse'); 
  title('WV(z) frac yr-1'); shading flat; colorbar
figure(2); clf; pcolor(save_lat,1:97,temp');  set(gca,'ydir','reverse'); 
  title('T(z) K yr-1'); shading flat; colorbar
figure(3); clf; plot(tracegases,save_lat,'-o','linewidth',2); grid
  set(gca,'ydir','reverse'); title('trace');
  hl = legend(thestr); set(hl,'fontsize',8);
  
if length(input_rates) == 2378
  ff = instr_chans('airs'); g = dogoodchan;
  figure(4); clf; plot(ff(g),input_rates(:,g)-fitted_rates(:,g)); title('AIRS : fitted biases'); grid
  figure(5); clf; plot(ff(g),input_rates(:,g)); title('AIRS : inputs'); grid
else
  ff = instr_chans('iasi'); g = 1:length(ff);
  figure(4); clf; plot(ff(g),input_rates(:,g)-fitted_rates(:,g)); title('IASI : fitted biases'); grid
  figure(5); clf; plot(ff(g),input_rates(:,g)); title('IASI : inputs'); grid
end

if N == 5
  [tracegases(3,:);tracegases_sigs(3,:)]'
  figure(3); axis([-5 +5 -90 +90]); 
  figure(6); clf; plot(water(3,:),plays,temp(3,:),plays,'r','linewidth',2); 
  figure(6); hl = errorbar_x(water(3,:),plays,water_sigs(3,:),'b'); set(hl,'linewidth',2); hold on
  figure(6); hl = errorbar_x(temp(3,:),plays,temp_sigs(3,:),'r'); set(hl,'linewidth',2); hold off
    axis([-0.15 +0.15 0 1000])
    set(gca,'ydir','reverse')
  title('TROPICS (b)water (frac)/yr (r)temp K/yr'); grid
end

iSave = input('save the quick results (-1/+1) : ');
if iSave > 0
  caName = input('Enter name of file : ');
  saver = ['save ' caName ' ff g input_rates fitted_rates water* temp* tracegases* thestr save_lat chanset lambdax plays'];
  eval(saver)
end