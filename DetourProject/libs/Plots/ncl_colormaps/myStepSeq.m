function cmap=myStepSeq
% function cmap=getNCLColmap(fname,n)
% fname = name of the file
% n = number of colors
% cmap = colormap, interpolated column by column.
fname = 'MPL_StepSeq.rgb';
n = 128;
a = mfilename('fullpath');
b = filesep;
ind = strfind(a,b);
cpath = a(1:ind(end));

[cname,ctype]=strtok(fname,'.');
ctype(1)=[];
fprintf('loading colormap %s%s/%s.%s\n',cpath,ctype,cname,ctype);
z=load([cpath,ctype,b,fname],'-ascii');


if max(z(:))>1, z=z/255; end
nz=length(z(:,1));
cmap=z;
if ~isempty(n)
    kc=linspace(1,nz,n);
    cmap=zeros(n,3);
    for k=1:3; cmap(:,k)=interp1(1:nz,z(:,k),kc); end
end

zn = round(size(cmap,1)/5);
for i =1:zn
    cmap2(i,:) = [1 1 1];
end
cmap = [cmap2;cmap];

end
%	ftag='cmap'; figWin(ftag);
%	hax=sameAxSubplt(1,2,ftag);
%	colormap(hax(1),z); colorbar(hax(1));
%	colormap(hax(2),cmap); colorbar(hax(2));
%   c = getNCLColmap('NCV_manga.rgb',256);
%   c = getNCLColmap('MPL_StepSeq.rgb',128);
		

