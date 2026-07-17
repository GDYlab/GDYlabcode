function cb = addcolorbar(axnow,gap,barwidth)

%ADDCOLORBAR add colorbar to the current gca, add to the right (eastoutside)
% restore the position of current gca
if nargin < 1
   axnow = gca; 
end
if nargin < 2
   gap = 0.005; 
end
if nargin < 3
   barwidth = 0.01; 
end
postemp = get(axnow,'position');
cbpos =  postemp;
cbpos(1) = postemp(1) + postemp(3) + gap;
cbpos(3) = barwidth;
cb = colorbar(axnow,'Position',cbpos);
set(axnow,'Position',postemp)
end

