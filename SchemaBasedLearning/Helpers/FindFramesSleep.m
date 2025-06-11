
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,5,5,5,5,5,4,4,4,4,4]);
FolderNumber=[1:25];
Frames=[];
for ii=FolderNumber
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load alltimes.mat
clusters = importdata('clusters');
elecreg = importdata('elecreg');
% reading and outputting celltype
[c1 c2] = textread('celltype','%s %s');
a1 = strrep(c1,'p','1'); % pyrs is 1
a1 = str2num(cell2mat(strrep(a1,'i','0'))); % int is 0, output is number matrix
a2 = strrep(c2,'b','0'); % bad is 0
a2 = strrep(a2,'g','1'); % good is 1
a2 = str2num(cell2mat(strrep(a2,'f','2'))); % food is 2, output is number matrix
a12 = [a1 a2]; % concatenating the 2 vectors
ind_PFC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==2;
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==1;

sess=1:size(alltimes,1);
t=0.005;   %binning of activity for calculating population firing rate
maxl=0.8; %max duration of frame
minl=0.03; %min duration of frame
ssd=0.015; %smoothing SD
thresh=2; %for frame detections in SDs
vlim=2;
threshmin=0;
Frames=[];
for nn=1:2
if nn ==1    
kkk=1:length(clusters);
kk=kkk(ind_HPC);
cl_interest=clusters(ind_HPC);
if sum(ind_HPC) < 10; continue; end
else
kkk=1:length(clusters);
kk=kkk(ind_PFC);
cl_interest=clusters(ind_PFC);
if sum(ind_PFC) < 10; continue; end
end


for j=1:length(sess)
frametime=[];frames=[];ncels=[];nspikes=[];
cltime=[];
for i=1:length(cl_interest)
clG=cl2mat(cl_interest{i});
% get time, velocity and ID
cltime{i}= [clG.featuredata((clG.featuredata(:,8)>=alltimes(sess(j),1) & clG.featuredata(:,8)<=alltimes(sess(j),2)),[8,11]) ...
           repmat(kk(i),size(clG.featuredata((clG.featuredata(:,8)>=alltimes(sess(j),1) & clG.featuredata(:,8)<=alltimes(sess(j),2)),8),1),1)];
clear clG
end


   cltimevec=[];
   
   for c=1:size(cltime,2)
   cltimevec=[cltimevec; cltime{c}];
   end

   cltimemat=[];
   [~,ind]=sort(cltimevec(:,1));
   cltimemat=cltimevec(ind,:);   

vel=cltimemat(:,2);
cltimemat=cltimemat(:,[1,3]);

% % remove spikes during movement    % kefei
% cltimemat(vel>vlim,:)=[];
% vel(vel>vlim,:)=[];
%-------------%

[~,tmpia]=unique(cltimemat(:,1),'stable');
[activity,timeloc]=hist(cltimemat(:,1),[alltimes(sess(j),1):t:alltimes(sess(j),2)]);

vel=interp1(cltimemat(tmpia,1),vel(tmpia),timeloc);

Fs = 1/(t);
npoints = round(4.*ssd.*Fs);
kernel = normpdf(linspace(-4*ssd, 4*ssd, 2*npoints+1), 0, ssd);
kernel = kernel./sum(kernel);
smoothedactivity= conv(activity, kernel,'same');

avgac=mean(smoothedactivity,'omitnan');
sdac=std(smoothedactivity,'omitnan');

indxvec=(smoothedactivity>(avgac+thresh*sdac));
z=(smoothedactivity-avgac)/sdac;
smactivity=[]; x=(1:size(smoothedactivity,2))';
smactivity=[smoothedactivity' indxvec' timeloc' x vel' z'];

up=[];
b=0;c=0;
for a=2:(size(smoothedactivity,2))
    if smoothedactivity(a)>(avgac+threshmin*sdac) & smoothedactivity(a-1)<=(avgac+threshmin*sdac)
       b=b+1;
       up(b,1)=a;
    end

end

down=[];c=0;
for a=up(1):(size(smoothedactivity,2))
    if smoothedactivity(a)<=(avgac+threshmin*sdac) & smoothedactivity(a-1)>(avgac+threshmin*sdac) %& a>up(a)
       c=c+1;
       down(c,1)=a-1;
    end

end

for a=1:size(down,1)
    precursors{a}=[smactivity(up(a):down(a),:)];
end


b=0;
for a=1:size(precursors,2)
    if precursors{a}(end,3)-precursors{a}(1,3)>=minl & precursors{a}(end,3)-precursors{a}(1,3)<=maxl & sum(precursors{a}(:,2)==1)>=1 %& sum(precursors{a}(:,5)>vlim)<1
        b=b+1;
        frametime{b}=precursors{a};
    end
end
% 
% for a=1:size(frametime,2)
%     mean_vel=(frametime{a}(1,5));
%     time_from_last_mov(a)=frametime{a}(1,3)-timeloc(find(timeloc<=frametime{a}(1,3) & vel>=mean_vel,1,'last'));
% end



for a=1:size(frametime,2)
    frames{a}=cltimemat((cltimemat(:,1)>=frametime{a}(1,3) & cltimemat(:,1)<=frametime{a}(end,3)),:);
end



for a=1:size(frames,2)
ncels(a,1)=size(unique(frames{a}(:,2)),1);
nspikes(a,1)=size(frames{a},1);

end

% histogram(ncels,'BinWidth',1)
% histogram(nspikes,'BinWidth',1)
if nn ==1
Frames(j).HPCframetime=frametime;
Frames(j).HPCframes=frames;
Frames(j).HPCncels=ncels;
Frames(j).HPCnspikes=nspikes;

else
Frames(j).PFCframetime=frametime;
Frames(j).PFCframes=frames;
Frames(j).PFCncels=ncels;
Frames(j).PFCnspikes=nspikes;   
end   
end 
end
save('FramesNolim.mat','Frames','-v7.3')
display(strcat('FrameFindingDoneForFolder',num2str(ii),'...'))
clear Frames nspikes ncels frames smoothedactivity vel cltimevec smactivity cltimemat
end


