function H=shadedErrorBarYLog10(xdata,yplays,xdatasigs,gah)

if size(xdata) ~= size(yplays)
  xdata = xdata';
end
if size(xdatasigs) ~= size(yplays)
  xdatasigs = xdatasigs';
end

faceAlpha = 0.25;
if nargin == 3
  patchColor = [0 0 1];
else
  patchColor = gah(1);
end

a = patch([xdata'+xdatasigs' fliplr(xdata'-xdatasigs')],...
          [log10(yplays')  fliplr(log10(yplays'))],1,...
              'facecolor',patchColor,...
              'edgecolor',patchColor,...
              'facealpha',faceAlpha);

  alpha(a,0.15)
hold on
  plot(xdata,log10(yplays),gah);
hold off

ticks = [0 1 2 3];
set(gca,'YTick',ticks);
set(gca,'YTickLabel',10.^ticks)
hold on
for i = 2:9
  for j = ticks
    plot([-1 1],j+log10(10*[i i]),'k:')
  end
end
%set(gca,'YTickLabel',strcat('10^',num2str(ticks,'%g')'),'Interpreter','latex');


set(gca,'ydir','reverse'); 
set(gcf, 'Renderer', 'OpenGL')

