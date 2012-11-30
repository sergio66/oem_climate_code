function H=shadedErrorBarY(x,y,errBar,lineProps,transparent)

% function H=shadedErrorBarY(x,y,errBarx,lineProps,transparent)
% same as shadedErrorBar but assumes errorBar is for x variable

%
% Purpose 
% Makes a 2-d line plot with a pretty shaded error bar made
% using patch. Error bar color is chosen automatically.
%
% Inputs
% x - vector of x values [optional, can be left empty]
% y - vector of y values or a matrix of n observations by m cases
%     where m has length(x);
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
% Outputs
% H - a structure of handles to the generated plot objects.     
%
%
% Examples
% x=randn(30,80); y=1:size(x,2);
% shadedErrorBarY(mean(x,1),y,std(x),'g');
% shadedErrorBarY(x,y,{@median,@std},{'r-o','markerfacecolor','r'});    
% shadedErrorBarY(x,[],{@median,@std},{'r-o','markerfacecolor','r'});    
%
% Overlay two transparent lines
% x=randn(30,80)*10; y=(1:size(x,2))-40;
% shadedErrorBarY(x,y,{@mean,@std},'-r',1); 
% hold on
% x=ones(30,1)*y; x=x+0.06*x.^2+randn(size(x))*10;
% shadedErrorBarY(x,y,{@mean,@std},'-b',1); 
% hold off
%
%
% Rob Campbell - November 2009


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Error checking    
error(nargchk(3,5,nargin))

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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Plot the main line. We plot this first in order to extract the RGB values
% for the line colour. I am not aware of a function that does this.
H.mainLine=plot(x,y,lineProps{:});


% Work out the color of the shaded region and associated lines
% Using alpha requires the render to be openGL and so you can't
% save a vector image. On the other hand, you need alpha if you're
% overlaying lines. We therefore provide the option of choosing alpha 
% or a de-saturated solid colour for the patch surface.

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

    
%Calculate the y values at which we will place the error bars
uE=x+errBar(1,:);
lE=x-errBar(2,:);



%Add the error-bar plot elements
holdStatus=ishold;
if ~holdStatus, hold on,  end


%Make the cordinats for the patch
xP=[lE,fliplr(uE)];
yP=[y,fliplr(y)];

%remove any nans otherwise patch won't work
yP(isnan(xP))=[];
xP(isnan(xP))=[];


H.patch=patch(xP,yP,1,'facecolor',patchColor,...
              'edgecolor','none',...
              'facealpha',faceAlpha);


%Make nice edges around the patch. 
H.edge(1)=plot(lE,y,'-','color',edgeColor);
H.edge(2)=plot(uE,y,'-','color',edgeColor);

%The main line is now covered by the patch object and was plotted first to
%extract the RGB value of the main plot line. I am not aware of an easy way
%to change the order of plot elements on the graph so we'll just remove it
%and put it back (yuk!)
delete(H.mainLine)
H.mainLine=plot(x,y,lineProps{:});


if ~holdStatus, hold off, end

