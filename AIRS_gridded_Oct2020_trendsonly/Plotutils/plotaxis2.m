function abcd = plotaxis(x0,y0);

if nargin == 0
  x0 = 0;
  y0 = 0;
elseif nargin == 1
  y0 = 0;
end

%% this cheesy function plots axis ....
%%% xabcd = plotaxis;

abcd = axis(gca);
x = [abcd(1) abcd(2)];
y = [abcd(3) abcd(4)];

minx = min(x); maxx = max(x);
miny = min(y); maxy = max(y);

x1 = min(minx); x2 = max(maxx); y1 = y0; y2 = y0;   
%line([x1 x2],[y0 y0]);
hold on; plot([x1 x2],[y0 y0],'k','LineWidth',2);

x1 = x0; x2 = x0; y1 = min(miny); y2 = max(maxy);   
%line([x0 x0],[y1 y2]);
hold on; plot([x0 x0],[y1 y2],'k','LineWidth',2);

axis(abcd);
hold off

grid on
