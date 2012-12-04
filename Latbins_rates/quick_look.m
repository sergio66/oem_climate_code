plevs = load('airslevels.dat');
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

iWhich = input('(1) OEM or (-1) LLS : ');

if iWhich == 1
  for ix = 1 : N
    loader = ['load ../Output/testx_' num2str(ix) '.mat'];
    eval(loader)
    params(ix,:)      = oem.finalrates;
    params_sigs(ix,:) = oem.finalsigs;
    fitted_rates(ix,:) = oem.fit;
    input_rates(ix,:)  = rateset.rates;
  end

  lambdax = struct;
    lambdax.lambda      = oem.lambda;
    lambdax.lambda_qst  = oem.lambda_qst;
    lambdax.lambda_Q1   = oem.lambda_Q1;
    lambdax.lambda_temp = oem.lambda_temp;

elseif iWhich == -1

  for ix = 1 : N
    loader = ['load ../Output/testx_' num2str(ix) '.mat'];
    eval(loader)
    params(ix,:)      = lls.finalrates;
    params_sigs(ix,:) = lls.finalsigs;
    fitted_rates(ix,:) = lls.fit;
    input_rates(ix,:)  = rateset.rates;
  end

%lambdax = struct;
%  lambdax.lambda      = lls.lambda;
%  lambdax.lambda_qst  = lls.lambda_qst;
%  lambdax.lambda_Q1   = lls.lambda_Q1;
%  lambdax.lambda_temp = lls.lambda_temp;
%clear oem rateset jacobian lls

end

thestr = jacobian.qstnames;
chanset = jacobian.chanset;
chanset = jacobian.chanset_used;

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

figure(1); clf; plot(water,1:97); set(gca,'ydir','reverse'); title('AIRS WV(z)')
figure(2); clf; plot(temp,1:97);  set(gca,'ydir','reverse'); title('AIRS T(z)')

figure(1); clf; pcolor(save_lat,1:97,water'); set(gca,'ydir','reverse'); 
  title('AIRS WV(z) frac yr-1'); shading flat; colorbar
figure(2); clf; pcolor(save_lat,1:97,temp');  set(gca,'ydir','reverse'); 
  title('AIRS T(z) K yr-1'); shading flat; colorbar

figure(1); clf; pcolor(save_lat,plays,water'); set(gca,'ydir','reverse'); 
  title('AIRS WV(z) frac yr-1'); shading flat; colorbar
figure(2); clf; pcolor(save_lat,plays,temp');  set(gca,'ydir','reverse'); 
  title('AIRS T(z) K yr-1'); shading flat; colorbar

figure(3); clf; plot(tracegases,save_lat,'-o','linewidth',2);
  hold on
  shadedErrorBarY(tracegases(:,1),save_lat,tracegases_sigs(:,1),'bo-',1); hold on
  shadedErrorBarY(tracegases(:,2),save_lat,tracegases_sigs(:,2),'go-',1); hold on
  shadedErrorBarY(tracegases(:,3),save_lat,tracegases_sigs(:,3),'ro-',1); hold on
  shadedErrorBarY(tracegases(:,4),save_lat,tracegases_sigs(:,4),'co-',1); hold on
  shadedErrorBarY(tracegases(:,5),save_lat,tracegases_sigs(:,5)*0,'mo-',1); hold on
  shadedErrorBarY(tracegases(:,6),save_lat,tracegases_sigs(:,6),'ko-',1); hold on
  hold off; grid; title('AIRS trace');
  hl = legend(thestr); 
  %legendlinestyles(hl,{'o' 'o'  'o' 'o' 'o' 'o'},{},{'b' 'g'  'r' 'c' 'm' 'k'})
  set(hl,'fontsize',8);
  axis([-5 +5 -90 +90]);

if length(input_rates) == 2378
  ff = instr_chans('airs'); g = dogoodchan;
  figure(4); clf; plot(ff(g),input_rates(:,g)-fitted_rates(:,g)); title('AIRS : fitted biases'); grid
    hold on
      plot(ff(chanset),nanmean(input_rates(:,chanset) - fitted_rates(:,chanset)),'ko-','linewidth',2);
    hold off

  figure(5); clf; plot(ff(g),input_rates(:,g)); title('AIRS : inputs'); grid
    hold on
      plot(ff(chanset),nanmean(input_rates(:,chanset)),'ko-','linewidth',2);
    hold off
else
  ff = instr_chans('iasi'); g = 1:length(ff);
  figure(4); clf; plot(ff(g),input_rates(:,g)-fitted_rates(:,g)); title('IASI : fitted biases'); grid
    hold on
      plot(ff(chanset),nanmean(input_rates(:,chanset) - fitted_rates(:,chanset)),'ko-','linewidth',2);
    hold off

  figure(5); clf; plot(ff(g),input_rates(:,g)); title('IASI : inputs'); grid
    hold on
      plot(ff(chanset),nanmean(input_rates(:,chanset)),'ko-','linewidth',2);
    hold off
end

figure(7); clf;
figure(7); hold on
if length(input_rates) == 2378
  ff = instr_chans('airs'); g = dogoodchan;
else
  ff = instr_chans('iasi'); g = 1:length(ff);
end
xx = find(abs(save_lat) <= 30);
if length(xx) > 0
  plot(ff(g),input_rates(xx,g)-fitted_rates(xx,g),'b');   
end
xx = find(save_lat > 30);
if length(xx) > 0
  plot(ff(g),input_rates(xx,g)-fitted_rates(xx,g),'g');   
end
xx = find(save_lat < -30);
if length(xx) > 0
  plot(ff(g),input_rates(xx,g)-fitted_rates(xx,g),'r');   
end
title('Biases (B) tropics (G) NML (R) SML'); grid
figure(7); hold off

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