function H=shadedErrorBarYLog10(xdata,yplays,xdatasigs,patch_color)

% a first pass sanity check for vector dimensions
if size(yplays,1) ~= 1
  yplays = yplays';
end
if size(xdata,1) ~= 1
  xdata = xdata';
end
if size(xdatasigs,1) ~= 1
  xdatasigs = xdatasigs';
end


% Now let's check to make sure everything is a vector in the same dimensions
if size(yplays,1) ~= 1
  error('Y values must be a vector with 1 in the second dimension, transpose the input')
end
if size(xdata) ~= size(yplays)
  error('X values are a different dimension than Y values')
end
if size(xdatasigs) ~= size(yplays)
  error('X Sigmas are a different dimensions than Y values')
end

faceAlpha = 0.25;
if nargin == 3
  patchColor = [0 0 1];
else
  patchColor = patch_color(1);
end

a = patch([xdata+xdatasigs fliplr(xdata-xdatasigs) xdata(1)+xdatasigs(1)],...
          [log10(yplays)  fliplr(log10(yplays)) log10(yplays(1))],1,...
              'facecolor',patchColor,...
              'edgecolor',patchColor,...
              'facealpha',faceAlpha);

  alpha(a,0.15)
hold on
  plot(xdata,log10(yplays),patchColor,'linewidth',2);
hold off

ticks = [0 1 2 3];
set(gca,'YTick',ticks);
set(gca,'YTickLabel',10.^ticks)
hold on
for i = 2:9
  for j = ticks
    line([-1 1],j+log10(10*[i i])-1,'color','k','linestyle',':')
  end
end
%set(gca,'YTickLabel',strcat('10^',num2str(ticks,'%g')'),'Interpreter','latex');

grid on
set(gca,'ydir','reverse'); 
%set(gcf, 'Renderer', 'OpenGL')

ax = axis;
axis([ax(1) ax(2) min(log10(yplays)) max(log10(yplays))])
