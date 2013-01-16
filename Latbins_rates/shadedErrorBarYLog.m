function H=shadedErrorBarYLog(x,y,errBar,lineProps,transparent)

% function H=shadedErrorBarYLog(x,y,errBarx,lineProps,transparent)
% same as shadedErrorBar  but assumes errorBar is for x variable
% same as shadedErrorBarY but much simpler, does a log plot in Y
%
% Purpose 
% Makes a 2-d line plot with a pretty shaded error bar made
% using patch. Error bar color is chosen automatically.
%
% Inputs
% x - vector of x values 
% y - vector of y values 
% errBarx - if a vector we draw symmetric errorbars. If it has a
%          size of [2,length(y)] then we draw asymmetric error bars
%          with row 1 being the upper bar and row 2 being the lower
%          bar. ** alternatively ** errBar can be a cellArray of
%          two function handles. The first defines which statistic
%          the line should be and the second defines the error
%          bar. 
% lineProps - [optional,'-k' by default] defines the properties of
%             the data line. e.g.:    
%             'or-', or {'-or','markerfacecolor',[1,0.2,0.2]}
% transparent - [optional, 0 by default] if ==1 the shaded error
%               bar is made transparent, which forces the renderer
%               to be openGl. However, if this is saved as .eps the
%               resulting file will contain a raster not a vector
%               image. 
%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Error checking 
error(nargchk(3,6,nargin))

%Process y using function handles if needed to make the error bar
%dynamically
if iscell(errBar) && ~isvector(x)
    fun1=errBar{1};
    fun2=errBar{2};
    errBar=fun2(x);
    x=fun1(x);
elseif ~iscell(errBar) && isvector(x)
    x=x(:)';
else
    error('2nd and 3rd input arguments are not compatible')
end

if isempty(y)
    y=1:length(x);
else
    y=y(:)';
end

if length(y) ~= length(x)
    error('inputs x and y are not of equal lengths')
end

%If only one error bar is specified then we will mirror it, turning it into
%both upper and lower bars. 
if length(errBar)==length(errBar(:))
    errBar=repmat(errBar(:)',2,1);
else
    f=find(size(errBar)==2);
    if isempty(f), error('errBar has the wrong size'), end
    if f==2, errBar=errBar'; end
end

if length(y) ~= length(errBar)
    error('inputs x and y must have the same length as errBar')
end

%Set default options
defaultProps={'-k'};
if nargin<4 || isempty(lineProps)
    lineProps=defaultProps; 
end
if ~iscell(lineProps)
    lineProps={lineProps}; 
end

if nargin<5 || ~isnumeric(transparent)
    transparent=0; 
end

H.mainLine=semilogy(x,y,lineProps{:});

col=get(H.mainLine,'color');
edgeColor=col+(1-col)*0.55;
patchSaturation=0.15; %How de-saturated or transparent to make the patch
if transparent
    faceAlpha=patchSaturation;
    patchColor=col;
    set(gcf,'renderer','openGL')
else
    faceAlpha=1;
    patchColor=col+(1-col)*(1-patchSaturation);
    set(gcf,'renderer','painters')
end

% size([x x])
% size([y y])
% whos x errBar
% plot(x); pause
% plot(errBar'); pause

semilogy(x,y); pause

a = patch([x+errBar(1,:) fliplr(x-errBar(2,:))],...
          [y fliplr(y)],1,...
              'facecolor',patchColor,...
              'edgecolor',patchColor,...
              'facealpha',faceAlpha);
  alpha(a,0.5)
  set(gca,'ydir','reverse');
  set(gcf, 'Renderer', 'OpenGL')
hold on
  semilogy(x,y,lineProps{:})
hold off