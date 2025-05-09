function cmap=flamecmap


a = mfilename('fullpath');
b = filesep;
ind = strfind(a,b);
cpath = a(1:ind(end));

z=load([cpath,'flamecmap.mat']);


cmap = z.cmap;

end


