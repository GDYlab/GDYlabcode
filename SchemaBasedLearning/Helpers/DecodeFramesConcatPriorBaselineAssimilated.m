clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load PlFieldTempNAll
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,5,5,5,5,5,4,4,4,4,4]);
for i=1:size(PlFieldTemp,2)
     plf{i,1}=PlFieldTemp(i).PlaceFields;
     animal(i,1)=PlFieldTemp(i).Animal;
     day(i,1)=PlFieldTemp(i).Day;
     session(i,1)=PlFieldTemp(i).session;
     track(i,1)=PlFieldTemp(i).Track;
     direction(i,1)=PlFieldTemp(i).Direction;
end

for ii=[1:25]
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load LapInfo
load mazel
load FramesNolim
load LFPDataV2

cluster=importdata('clusters');
elecreg=importdata('elecreg');
tbin=0.02;
GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})
dn=pwd;
fstruct = dir('*FrameHPCShuff*.mat');
if ii==6;
fstruct([1:10], :) = fstruct([1,[3:10],2], :);
end
fstruct1=array2table(fstruct);

tic;
if ~isempty(fstruct)
for i=1:length(fstruct)
   tic;    
    
     filename=fstruct(i).name;
     load(filename);
     dotLocations = find(filename == '.');
     filename1 = filename(1:dotLocations(1)-1)
     sesslocs=find(filename == 'f');
     sess = str2double(filename(sesslocs(end)+1:dotLocations(1)-1));
     
     assignin ( 'caller', 'frameDecode1',eval(filename1));
     clear((filename1));
    
FrameDecode=[]; plf1=[];   
 
if any(ii==[1:5,21:25]) | ii==6 | ii==15; dd=0; else dd=1;end
plfConcat=[]; 
for ll=[1,length(mazel)-dd]    
plf1=[];  
for kk=1:2
for jj=1:8
plf2=[];    
plf1=plf{animal==Animals(ii) & day==Days(ii) & session==ll & track==jj & direction==kk};
plf2(:,1:70)=(plf1(:,1:end-1)+plf1(:,2:end))/2;
plfConcat=[plfConcat, plf2];
end
end 
end
           if isfield(Frames,'HPCframetime'); 
         
           for nn=1:length(frameDecode1)  
           tvec=squeeze(frameDecode1(nn).t(1,1,1,:));
           sp=[]; tvec3=[];
           for a=nn:length(Frames(sess).HPCframes)
           tvec2=Frames(sess).HPCframetime{1,a}(:,3); 
           binns=nanmin(tvec2):tbin:nanmax(tvec2);
           binedges=(binns(1:end-1)+binns(2:end))/2;    
           if tvec(1)==binedges(1); sp= Frames(sess).HPCframes{1, a};tvec3=tvec2; continue; end
           end
           [dprob binedges]=BayesDecode_Frames(tvec3,plfConcat,sp,tbin);
           FrameDecode.HPCProb{nn}=dprob;
           FrameDecode.HPCt{nn}=tvec';
           end           
           end
         
      
temp_var2 = strcat( 'FrameHPCDecodeConcatComb',num2str(sess) );
assignin ( 'caller', temp_var2, FrameDecode); 
display(strcat('DoneForFolder',num2str(ii),'ForFrameSess',num2str(sess)));
toc
GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})
save(strcat(temp_var2,'.mat'),temp_var2,'-v7.3')
clear(temp_var2)
cd(dn)
  end   
 
end


GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})
fstruct = dir('*FramePFCShuff*.mat');
if ii==6;
fstruct([1:10], :) = fstruct([1,[3:10],2], :);
end
fstruct1=array2table(fstruct);
% if exist('FrameDecode.mat')
% delete('FrameDecode.mat');
% end
tic;
if ~isempty(fstruct)
for i=1:length(fstruct)
   tic;    
     ss=i
     filename=fstruct(i).name;
     load(filename);
     dotLocations = find(filename == '.');
     filename1 = filename(1:dotLocations(1)-1)
     sesslocs=find(filename == 'f');
     sess = str2double(filename(sesslocs(end)+1:dotLocations(1)-1));
     assignin ( 'caller', 'frameDecode1',eval(filename1));
     clear((filename1));
    
FrameDecode=[]; plf1=[];   
 
if any(ii==[1:5,21:25]) | ii==8 | ii==13 | ii==18 | ii==6 | ii==15; dd=0; else dd=1;end
plfConcat=[]; 
 for ll=[1,length(mazel)-dd]    
plf1=[];  
for kk=1:2
for jj=1:8
plf2=[];    
plf1=plf{animal==Animals(ii) & day==Days(ii) & session==ll & track==jj & direction==kk};
plf2(:,1:70)=(plf1(:,1:end-1)+plf1(:,2:end))/2;
plfConcat=[plfConcat, plf2];
end
end 
 end    
           if isfield(Frames,'PFCframetime'); 
         
           for nn=1:length(frameDecode1)  
           tvec=squeeze(frameDecode1(nn).t(1,1,1,:));
           sp=[]; tvec3=[];
           for a=nn:length(Frames(sess).PFCframes)
           tvec2=Frames(sess).PFCframetime{1,a}(:,3); 
           binns=nanmin(tvec2):tbin:nanmax(tvec2);
           binedges=(binns(1:end-1)+binns(2:end))/2;    
           if tvec(1)==binedges(1); sp= Frames(sess).PFCframes{1, a};tvec3=tvec2; continue; end
           end
           [dprob binedges]=BayesDecode_Frames(tvec3,plfConcat,sp,tbin);
           FrameDecode.PFCProb{nn}=dprob;
           FrameDecode.PFCt{nn}=tvec';
           end           
           end
            
 
      
temp_var2 = strcat( 'FramePFCDecodeConcatComb',num2str(sess) );
assignin ( 'caller', temp_var2, FrameDecode); 
display(strcat('DoneForFolder',num2str(ii),'ForFrameSess',num2str(sess)));
toc
GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})
save(strcat(temp_var2,'.mat'),temp_var2,'-v7.3')
clear(temp_var2)
cd(dn);
 end
end

end


% for i=9
%      filename=strcat('FrameHPCShuff',num2str(i),'.mat');
%      load(filename);
%      dotLocations = find(filename == '.');
%      filename1 = filename(1:dotLocations(1)-1)
%      sesslocs=find(filename == 'f');
%      sess = str2double(filename(sesslocs(end)+1:dotLocations(1)-1));
%      
%      assignin ( 'caller', 'frameDecode1',eval(strcat('FrameHPCShuff',num2str(2))));
%      clear (filename1)
%      
% 
% temp_var2 = strcat( 'FrameHPCShuff',num2str(i) );
% assignin ( 'caller', temp_var2, frameDecode1); 
% save(strcat(temp_var2,'.mat'),temp_var2,'-v7.3')
% display(strcat('DoneForFolder',num2str(6),'ForFrameSess',num2str(i)));
% clear all;
% end
% 
% for i=9
%      filename=strcat('FramePFCShuff',num2str(i),'.mat');
%      load(filename);
%      dotLocations = find(filename == '.');
%      filename1 = filename(1:dotLocations(1)-1)
%      sesslocs=find(filename == 'f');
%      sess = str2double(filename(sesslocs(end)+1:dotLocations(1)-1));
%      
%      assignin ( 'caller', 'frameDecode1',eval(strcat('FramePFCShuff',num2str(2))));
%      clear (filename1)
%      
% 
% temp_var2 = strcat( 'FramePFCShuff',num2str(i) );
% assignin ( 'caller', temp_var2, frameDecode1); 
% save(strcat(temp_var2,'.mat'),temp_var2,'-v7.3')
% display(strcat('DoneForFolder',num2str(6),'ForFrameSess',num2str(i)));
% clear all;
% end
