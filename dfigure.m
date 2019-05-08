function varargout = dfigure(varargin)
% My version of figure. Additional optional arguments:
% 'DPosition' : multiplies the figure's [height,width] by a value. Good for subplots.
% 'DName'     : finds figures with this name and closes them, then gives new figure that name

dpl = find(cellfun(@(x)isequal(x,'DPosition'),varargin));
dnl = find(cellfun(@(x)isequal(x,'DName'),varargin));
wsl = find(cellfun(@(x)isequal(x,'WindowStyle'),varargin));
xcld=unique([dpl,dnl,dpl+1,dnl+1]); fargin=varargin(~ismember([1:length(varargin)],xcld));

% Optional argument DPosition
sf=[1,1];
if ~isempty(dpl)
    sf=varargin{dpl(1)+1};
end
    
% Optional argument DName
ws='normal'; if ~isempty(wsl), ws=varargin{wsl(1)+1}; end
if ~isempty(dnl)
    myname = varargin{dnl(1)+1};
    myh = findobj('Type','figure','-and','Name',myname);
    if ~isempty(myh)
        ws = myh.WindowStyle;
        close(myh);
    end
    fargin{end+1}='Name'; fargin{end+1}=myname;
end
    

% Make figure
fn = figure(fargin{:},'WindowStyle',ws);

my_FontName = 'Tahoma';
my_LineWidth = 1.0;
my_Color = 'w';
my_ColorOrder = dColor([]);

set(gcf,'Color',my_Color);
set(gcf,'DefaultAxesColorOrder',my_ColorOrder);
set(gcf,'DefaultLineLineWidth',my_LineWidth);
set(gcf,'DefaultAxesLineWidth',my_LineWidth);
set(gcf,'DefaultAxesFontName',my_FontName);
set(gcf,'DefaultTextFontName',my_FontName);

cpos = get(gcf,'Position');
wd=cpos(3); ht=cpos(4);
set(gcf,'Position',[cpos(1)+wd/2-wd*sf(2)/2,cpos(2)+ht/2-ht*sf(1)/2,wd*sf(2),ht*sf(1)]);
movegui(gcf,'center');

varargout = {fn};

end

