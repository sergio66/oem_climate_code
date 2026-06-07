figure(075); close
figure(076); close
figure(174); close
figure(175); close
figure(176); close

fig13_400mb.iSmooth = iSmooth;
fig13_400mb.rlat65  = rlat65;
fig13_400mb.rlon73  = rlon73;
fig13_400mb.umbc_wv_400 = smoothn((reshape(fracWV(i400,:)',72,64)') ,1);
fig13_400mb.era5_wv_400 = smoothn((reshape(era5_wvrate(i400,:)',72,64)') ,1);

cxmax = 0.0125;
cxmax = 0.0150;
cxmax = 0.0149;

if iSmooth >= 1
  aslmap(75,rlat65,rlon73,smoothn((reshape(fracWV(i400,:)',72,64)') ,1),     [-90 +90],[-180 +180]); caxis([-1 +1]*cxmax); colormap(llsmap5); % title('dfracWV/dt AIRS\_RT 400 mb'); 
  aslmap(76,rlat65,rlon73,smoothn((reshape(era5_wvrate(i400,:)',72,64)') ,1),[-90 +90],[-180 +180]); caxis([-1 +1]*cxmax); colormap(llsmap5); % title('dfracWV/dt ERA5 400 mb');    
elseif iSmooth == 0.5
  aslmap(75,rlat65,rlon73,smoothdata((reshape(fracWV(i400,:)',72,64)'),"movmean",5),     [-90 +90],[-180 +180]); caxis([-1 +1]*cxmax); colormap(llsmap5); % title('dfracWV/dt AIRS\_RT 400 mb'); 
  aslmap(76,rlat65,rlon73,smoothdata((reshape(era5_wvrate(i400,:)',72,64)'),"movmean",5),[-90 +90],[-180 +180]); caxis([-1 +1]*cxmax); colormap(llsmap5); % title('dfracWV/dt ERA5 400 mb');    
elseif iSmooth <= 0
  aslmap(75,rlat65,rlon73,reshape(fracWV(i400,:)',72,64)',     [-90 +90],[-180 +180]); caxis([-1 +1]*cxmax); colormap(llsmap5); % title('dfracWV/dt AIRS\_RT 400 mb'); 
  aslmap(76,rlat65,rlon73,reshape(era5_wvrate(i400,:)',72,64)',[-90 +90],[-180 +180]); caxis([-1 +1]*cxmax); colormap(llsmap5); % title('dfracWV/dt ERA5 400 mb');    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aslmap(175,fig13_400mb.rlat65,fig13_400mb.rlon73,fig13_400mb.umbc_wv_400,[-90 +90],[-180 +180]); caxis([-1 +1]*cxmax); colormap(llsmap5);
    set(gca,'fontsize',18)
    text(-0.25,-1.625,'yr^{-1}','fontsize',14)
aslmap(176,fig13_400mb.rlat65,fig13_400mb.rlon73,fig13_400mb.era5_wv_400,[-90 +90],[-180 +180]); caxis([-1 +1]*cxmax); colormap(llsmap5);
    set(gca,'fontsize',18)
    text(-0.25,-1.625,'yr^{-1}','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jgr_locate = '/home/sergio/MATLABCODE/oem_pkg_run/MATFILES_for_JGR_trends_paper/fig13.mat';
jgr_locate = '/home/sergio/git/JGR_July2025/PLOTTER/fig13.mat';
if exist(jgr_locate)
  oldfig13 = load(jgr_locate);
  aslmap(275,oldfig13.fig13_400mb.rlat65,fig13_400mb.rlon73,fig13_400mb.umbc_wv_400,[-90 +90],[-180 +180]); caxis([-1 +1]*cxmax); colormap(llsmap5);
    set(gca,'fontsize',18)
    text(-0.25,-1.625,'yr^{-1}','fontsize',14)
  aslmap(276,oldfig13.fig13_400mb.rlat65,fig13_400mb.rlon73,fig13_400mb.era5_wv_400,[-90 +90],[-180 +180]); caxis([-1 +1]*cxmax); colormap(llsmap5);
    set(gca,'fontsize',18)
    text(-0.25,-1.625,'yr^{-1}','fontsize',14)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(174); clf

  set(gcf,'resize','off')
  set(gcf, 'Position',  [100, 100,  560*2, 420]);
  
  ta = tiledlayout(1,2,'TileSpacing','None', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  
  tafov(1) = nexttile;
    if iSmooth <= 0
      aslmapSergio(rlat65,rlon73,smoothn((reshape(fracWV(i400,:)',72,64)') ,1),           [-90 +90],[-180 +180]); caxis([-1 +1]*cxmax); colormap(llsmap5); % title('dfracWV/dt AIRS\_RT 400 mb');
    elseif iSmooth == 0.5
      aslmapSergio(rlat65,rlon73,smoothdata((reshape(fracWV(i400,:)',72,64)'),"movmean",5),[-90 +90],[-180 +180]); caxis([-1 +1]*cxmax); colormap(llsmap5); % title('dfracWV/dt AIRS\_RT 400 mb');
    else
      error('oops')
    end
    cb = colorbar('horizontal');
    currentPosition = get(cb, 'Position');
    currentPosition(1) = currentPosition(1) + 0.05;
    currentPosition(2) = currentPosition(2) - 0.1;    
    currentPosition(3) = 0.4; % Set new width
    set(cb, 'Position', currentPosition);
    text(-0.1,-1.5,'yr^{-1}','fontsize',14);
    set(gca,'fontsize',14);
    
  tafov(2) = nexttile;
    if iSmooth <= 0  
      aslmapSergio(rlat65,rlon73,smoothn((reshape(era5_wvrate(i400,:)',72,64)') ,1),           [-90 +90],[-180 +180]); caxis([-1 +1]*cxmax); colormap(llsmap5); % title('dfracWV/dt ERA5 400 mb');
    elseif iSmooth == 0.5      
      aslmapSergio(rlat65,rlon73,smoothdata((reshape(era5_wvrate(i400,:)',72,64)'),"movmean",5),[-90 +90],[-180 +180]); caxis([-1 +1]*cxmax); colormap(llsmap5); % title('dfracWV/dt AIRS\_RT 400 mb');
    else
      error('oops')
    end      
    cb = colorbar('horizontal');
    currentPosition = get(cb, 'Position');
    currentPosition(1) = currentPosition(1) + 0.05;
    currentPosition(2) = currentPosition(2) - 0.1;        
    currentPosition(3) = 0.4; % Set new width
    set(cb, 'Position', currentPosition);
    text(-0.1,-1.5,'yr^{-1}','fontsize',14);    
    set(gca,'fontsize',14);
    
  strspace = 'none';
  %strspace = 'tight';
  %strspace = 'compact';

  % Get rid of all extra space I can
  ta.Padding = strspace;
  ta.TileSpacing = strspace;
