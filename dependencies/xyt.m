function xyt(varargin)
    if nargin>3
        [ax,xlbl,ylbl,ttl]=deal(varargin{:});
        xlabel(ax,xlbl);
        ylabel(ax,ylbl);
        title(ax,ttl);
    else
        [xlbl,ylbl,ttl]=deal(varargin{:});
        xlabel(xlbl);
        ylabel(ylbl);
        title(ttl);
    end
end