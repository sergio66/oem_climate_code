function H=shadedErrorBarYLog10(xdata,yplays,xdatasigs,plotopt)

%% optimum for yplays = 1013 mb to 1 mb (hence the log10)

if size(xdata) ~= size(yplays)
  xdata = xdata';
end
if size(xdatasigs) ~= size(yplays)
  xdatasigs = xdatasigs';
end

[mm,nn] = size(xdata);
if mm > 1 & nn > 1
  error('xdata must be 1xN or Nx1');
end

faceAlpha = 0.5;
if nargin == 3
  patchColor = [0 0 1];
  plotopt = ['bo-'];
else
  patchColor = plotopt(1);
end

a = patch([xdata'+xdatasigs' fliplr(xdata'-xdatasigs')],...
          [log10(yplays')  fliplr(log10(yplays'))],1,...
              'facecolor',patchColor,...
              'edgecolor',patchColor,...
              'facealpha',faceAlpha);
  alpha(a,0.5)

hold on
  plot(xdata,log10(yplays),plotopt);
hold off

set(gca,'ydir','reverse'); 
set(gcf, 'Renderer', 'OpenGL')

