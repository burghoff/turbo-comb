function rgb = dColor(cnum,num_colors)

my_ColorOrder = [55,126,184;...
                 228,26,28;...
                 77,175,74;...
                 152,78,163;...
                 247,129,191;...
                 120,120,120;...
                 255,127,0;...
                 166,86,40;...
                 255,215,0]/255;

if isempty(cnum)
    rgb = my_ColorOrder;
elseif isequal(cnum,'demo')
    dfigure; axes; set(gca,'ColorOrder',my_ColorOrder);
    for i1=1:size(my_ColorOrder,1)
        plot([0,1],i1*[1,1],'Color',my_ColorOrder(i1,:)); hold all;
    end
else
    if nargin<2, num_colors = size(my_ColorOrder,1);
    else,        num_colors = min(num_colors, size(my_ColorOrder,1)); end
    rgb = my_ColorOrder(mod(cnum-1,num_colors)+1,:);
end

             
end

