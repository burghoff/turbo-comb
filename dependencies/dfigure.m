function varargout = dfigure(varargin)
% David's version of figure. Additional optional arguments:
% 'DPosition' : multiplies the figure's [height,width] by a value. Good for subplots.
% 'DName'     : finds figures with this name and closes them, then gives new figure that name
% 'NumColors' : when set to n, only uses the first n of the colors

dpl = find(cellfun(@(x)isequal(x,'DPosition'),varargin));
dnl = find(cellfun(@(x)isequal(x,'DName'),varargin));
ncl = find(cellfun(@(x)isequal(x,'NumColors'),varargin));
wsl = find(cellfun(@(x)isequal(x,'WindowStyle'),varargin));
xcld=unique([dpl,dnl,ncl,dpl+1,dnl+1,ncl+1]); fargin=varargin(~ismember([1:length(varargin)],xcld));

% Optional argument DPosition
sf=[1,1];
if ~isempty(dpl)
    sf=varargin{dpl(1)+1};
end
    
% Optional argument DName
ws='normal'; if ~isempty(wsl), ws=varargin{wsl(1)+1}; end
already_exists = 0;
if ~isempty(dnl)
    myname = varargin{dnl(1)+1};
    myh = findobj('Type','figure','-and','Name',myname);
    if ~isempty(myh)
        ws = myh.WindowStyle;
        % close(myh);
        clf(myh);
        already_exists = 1;
        fn = myh;
    end
    fargin{end+1}='Name'; fargin{end+1}=myname;
end

% Make figure
if ~already_exists
    fn = figure(fargin{:},'WindowStyle',ws);
end

try
    theme(fn,"light")
end

my_FontName = 'Arial';
my_LineWidth = 0.75;
my_Color = 'w';
my_ColorOrder = dColor([]);

% Optional argument NumColors
nc = size(my_ColorOrder,1);
if ~isempty(ncl)
    nc = min(nc,max(1,varargin{ncl+1}));
end

set(gcf,'Color',my_Color);
set(gcf,'DefaultAxesColorOrder',my_ColorOrder(1:nc,:));
set(gcf,'DefaultLineLineWidth',my_LineWidth);
set(gcf,'DefaultAxesLineWidth',my_LineWidth);
set(gcf,'DefaultAxesFontName',my_FontName);
set(gcf,'DefaultTextFontName',my_FontName);
set(gcf,'DefaultAxesColor','none');
% set(gcf,'defaultAxesXGrid','on','defaultAxesYGrid','on');

if ~already_exists
    cpos = get(gcf,'Position');
    wd=cpos(3); ht=cpos(4);
    set(gcf,'Position',[cpos(1)+wd/2-wd*sf(2)/2,cpos(2)+ht/2-ht*sf(1)/2,wd*sf(2),ht*sf(1)]);
    movegui(gcf,'center');
end

varargout = {fn};
end

