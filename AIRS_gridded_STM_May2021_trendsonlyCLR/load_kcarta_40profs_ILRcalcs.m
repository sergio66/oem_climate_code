addpath /home/sergio/MATLABCODE/PLOTMISC

iCKD = 1;
iCKD = 25;
iCKD = 32;

fprintf(1,'iCKD = %2i \n',iCKD)

RRTMbands = [10 250 500 630 700 820 980 1080 1180 1390 1480 1800 2080 2250 2380 2600 3000];

dataDir = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/40profiles_uplook_flux_jacs/';

%for ii = 1 : 40
%  flux = load([dataDir 'CKD32/CALC0_JACOBIANS/allkc5bands_prof' num2str(ii) '.mat']);
%  flux0(ii,:)        = flux.dall;
%end

%% read in WV,CO2,O3,CH4,T jacs for 40 profiles
for ii = 1 : 40
  flux = load([dataDir 'CKD32/CALC0_JACOBIANS/allkc5bands_prof' num2str(ii) '.mat']);
  figure(1); plot(flux.wall,flux.jall(:,1)-flux.dall(:,1))
  figure(2); plot(flux.wall,flux.jall(:,2)-flux.dall(:,1))
  jacILR_WV(ii,:)  = flux.jall(:,1)-flux.dall;
  jacILR_CO2(ii,:) = flux.jall(:,2)-flux.dall;
  jacILR_O3(ii,:)  = flux.jall(:,3)-flux.dall;
  jacILR_CH4(ii,:) = flux.jall(:,4)-flux.dall;
  jacILR_T(ii,:)   = flux.jall(:,5)-flux.dall;
  rad0(ii,:)         = flux.dall;
end

%% read in WV,CO2,O3,CH4,T change in rads for 40 profiles x 4 perts = 140
for ii = 1 : 160
  flux = load([dataDir 'CKD32/PERTURBATIONS/allkc5bands_prof' num2str(ii) '.mat']);
  deltaradILR(ii,:)  = flux.dall;
end

[fc,qcjacWV]  = quickconvolve(flux.wall,jacILR_WV,0.25,0.25);
[fc,qcjacCO2] = quickconvolve(flux.wall,jacILR_CO2,0.25,0.25);
[fc,qcjacO3]  = quickconvolve(flux.wall,jacILR_O3,0.25,0.25);
[fc,qcjacCH4] = quickconvolve(flux.wall,jacILR_CH4,0.25,0.25);
[fc,qcjacT]   = quickconvolve(flux.wall,jacILR_T,0.25,0.25);
[fc,qcRad0]   = quickconvolve(flux.wall,rad0,0.25,0.25);
[fc,qcRadPertWV]   = quickconvolve(flux.wall,deltaradILR((1:40)+0*40,:)-rad0,0.25,0.25);
[fc,qcRadPertCO2]  = quickconvolve(flux.wall,deltaradILR((1:40)+1*40,:)-rad0,0.25,0.25);
[fc,qcRadPertT]    = quickconvolve(flux.wall,deltaradILR((1:40)+2*40,:)-rad0,0.25,0.25);
[fc,qcRadPertAll3] = quickconvolve(flux.wall,deltaradILR((1:40)+3*40,:)-rad0,0.25,0.25);

figure(1); clf; plot(fc,nanmean(qcjacWV,2),'b',fc,nanmean(qcjacT,2),'r',fc,nanmean(qcjacCO2,2),'k','linewidth',2); hl = legend('WV','T','CO2','location','best');
xlabel('Wavenumber cm-1'); ylabel('<Flux Jacobian>');

figure(2); clf; plot(fc,nanmean(qcRad0,2),'b',fc,ttorad(fc,289),'g','linewidth',2)
xlabel('Wavenumber cm-1'); ylabel('<Flux>');
shadeA(0550,0,250,140,'g',0.10);

%{
line([0550 0550],[0,140],'color','k');
line([0800 0800],[0,140],'color','k');
line([0980 0980],[0,140],'color','k');
line([1080 1080],[0,140],'color','k');
line([1080 1080],[0,140],'color','k');
line([1250 1250],[0,140],'color','k');
line([1350 1350],[0,140],'color','k');
line([1950 1950],[0,140],'color','k');
shadeA(0010,0,540,140,'b',0.10);
shadeA(0550,0,250,140,'g',0.10);
shadeA(0800,0,180,140,'b',0.10);
shadeA(0980,0,100,140,'y',0.10);
shadeA(1080,0,170,140,'b',0.10);
shadeA(1250,0,170,140,'y',0.10);
shadeA(1350,0,600,140,'b',0.10);
%}

figure(3); clf; plot(fc,nanmean(qcRadPertWV,2),'b',fc,nanmean(qcRadPertT,2),'r',fc,nanmean(qcRadPertCO2,2),'k','linewidth',2); hl = legend('WV','T','CO2','location','best');
xlabel('Wavenumber cm-1'); ylabel('<Avg Rad Change/yr>');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error('do not need the below')
if ~exist('ilrX')
  disp('need to load 40 x 4 perturbation files , "+" indicates 10, "o" indicates 100')
  dir0 = [dataDir 'CKD' num2str(iCKD) '/PERTURBATIONS/'];
  for ii = 1 : 160
    if mod(ii,100) == 0
      fprintf(1,'o');
    elseif mod(ii,10) == 0
      fprintf(1,'+');
    else
      fprintf(1,'.');
    end
    fin = [dir0 '/allkc5bands_prof' num2str(ii) '.mat'];
    junk = load(fin);
    ilrX(ii) = trapz(junk.wall,junk.fluxall(:,1))/1000;  %% to to W/m2
    for bbb = 1 : length(RRTMbands)-1
      boo = find(junk.wall >= RRTMbands(bbb) & junk.wall < RRTMbands(bbb+1));
      ilrXbands(ii,bbb) = trapz(junk.wall(boo),junk.fluxall(boo,1))/1000;  %% to to W/m2
    end
  end
  fprintf(1,'\n');
end

if ~exist('surfprof')
  disp('need to load 40 files , "+" indicates 10')
  dir0 = [dataDir 'CKD' num2str(iCKD) '/CALC0_JACOBIANS'];
  for ii = 1 : 40
    if mod(ii,10) == 0
      fprintf(1,'+');
    else
      fprintf(1,'.');
    end
    fin = [dataDir '/allkc5bands_prof' num2str(ii) '.mat'];
    junk = load(fin,'surf_properties');
    surfprof(ii,:) = junk.surf_properties;
  
    junk = load(fin);
    rad0(ii) = trapz(junk.wall,junk.dall);
    ilr0(ii) = trapz(junk.wall,junk.fluxall(:,1))/1000;  %% to to W/m2

    for bbb = 1 : length(RRTMbands)-1
      boo = find(junk.wall >= RRTMbands(bbb) & junk.wall < RRTMbands(bbb+1));
      ilr0bands(ii,bbb) = trapz(junk.wall(boo),junk.fluxall(boo,1))/1000;  %% to to W/m2
    end

    for jjj = 1 : 6
      radpert(ii,jjj) = trapz(junk.wall,junk.jall(:,jjj));    
    end
  end
  fprintf(1,'\n');
end
