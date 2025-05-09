% generate a figure whose size is half of an US letter 
% h: ratio of height of full US letter paper
% w:ratio of width of full US letter paper
% kefei
function f=figure_letter(h,w,varargin)
if nargin<1
    h=1/2;
end
if nargin<2
    w=1;
end
f=figure(varargin{:});
set(f,'paperunits','inches');
set(f,'papersize',[8.5,11]);
u=get(f,'units');
set(f,'units','inches');
tmp=get(f,'papersize');
set(f,'position',[1,1,tmp(1)*w,tmp(2)*h]);
set(f,'units',u);
