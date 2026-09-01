function [x,y,c] = get_fig_image_data(figname,iNewFigOrNot);

if nargin == 1
  iNewFigOrNot = -1;  %% use existing figure, do not create new one
end 

if iNewFigOrNot > 0
  hgload(figname);
else
  clf;
  h = gca;

  fig = openfig(figname,'invisible'); % a figure you have already saved
  copyobj(allchild(get(fig,'CurrentAxes')),h);

  %newfigH = findall(gca); %get the handles of the figure to copy to
  %delete(newfigH(2:end)); %delete all but the figure itself
  %children = allchild(fig); %move children form one fig to the other
  delete(fig); %delete the now empty figure

end
  
h = findobj(gca,'type','surface');
x = get(h,'xdata');
y = get(h,'ydata');
c = get(h,'cdata');

if length(x) == 0   
  h = findobj(gca,'type','image');
  x = get(h,'xdata');
  y = get(h,'ydata');
  c = get(h,'cdata');
end

if iNewFigOrNot < 0
  %whos x y c
  pcolor(x,y,c); shading interp; colorbar
end  


