%%%%%%%%%%%%%%%%%%%%%%%%%
figure(70); clf; 
subplot(221); imagesc(Rallchans(2:6,:)'); colorbar; colormap jet; title('All Chans')
cx = caxis; cx = [0 1]; caxis(cx);
%set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

subplot(222); imagesc(R15umchans(2:6,:)'); colorbar; colormap jet; title('15um Chans')
cx = caxis; cx = [0 1]; caxis(cx);
%set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
%set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

subplot(223); imagesc(Rwinchans(2:6,:)'); colorbar; colormap jet; title('Window Chans')
cx = caxis; cx = [0 1]; caxis(cx);
set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

subplot(224); imagesc(RWVchans(2:6,:)'); colorbar; colormap jet; title('WV Chans')
cx = caxis; cx = [0 1]; caxis(cx);
set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
%set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(71); clf; 
subplot(221); imagesc(Mallchans(2:6,:)'); colorbar; colormap(usa2); 
cx = caxis; cx = [-max(abs(cx)) +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [-0.03 +0.03]; caxis(cx);
title('All Chans')
%set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

subplot(222); imagesc(M15umchans(2:6,:)'); colorbar; colormap(usa2);
cx = caxis; cx = [-max(abs(cx)) +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [-0.03 +0.03]; caxis(cx);
title('15um Chans')
%set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
%set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

subplot(223); imagesc(Mwinchans(2:6,:)'); colorbar; colormap(usa2);
cx = caxis; cx = [-max(abs(cx)) +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [-0.03 +0.03]; caxis(cx);
title('Window Chans')
set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

subplot(224); imagesc(MWVchans(2:6,:)'); colorbar; colormap(usa2);
cx = caxis; cx = [-max(abs(cx)) +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [-0.03 +0.03]; caxis(cx);
title('WV Chans')
set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
%set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(72); clf; 
subplot(221); imagesc(Sallchans(2:6,:)'); colorbar; colormap jet; title('All Chans')
cx = caxis; cx = [0 +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [0 0.02]; caxis(cx);
%set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

subplot(222); imagesc(S15umchans(2:6,:)'); colorbar; colormap jet; title('15um Chans')
cx = caxis; cx = [0 +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [0 0.02]; caxis(cx);
%set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
%set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

subplot(223); imagesc(Swinchans(2:6,:)'); colorbar; colormap jet; title('Window Chans')
cx = caxis; cx = [0 +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [0 0.02]; caxis(cx);
set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

subplot(224); imagesc(SWVchans(2:6,:)'); colorbar; colormap jet; title('WV Chans')
cx = caxis; cx = [0 +max(abs(cx))]; caxis(cx);
cx = caxis; cx = [0 0.02]; caxis(cx);
set(gca,'xtick',[1:5],'xticklabel',varnames(2:6))
%set(gca,'ytick',[1:4],'yticklabel',{'Global','Tropical','Midlat','Polar'})

%%%%%%%%%%%%%%%%%%%%%%%%%
