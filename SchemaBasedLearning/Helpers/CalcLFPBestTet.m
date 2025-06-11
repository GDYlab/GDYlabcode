function CalcLFPBestTet(FolderNumber)
% read LFP from CSC files of the best channel tetrode in sleep session
disp(strcat('CalcLFPBestTet start ...'));
%% determine the rat number
% for rat 1  = CSC54
% for rat 2 = CSC17
% for rat 3 = TT1 and TT6 for HPC 
% for rat 4= CSC 25,26,27
% for rat 5= CSC4,5,6,7
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,5,5,5,5,5,4,4,4,4,4]);


for ii=FolderNumber
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
extimes= importdata('extimes_all');
NumSess=size(extimes,1);
alltimes=[];
for k=1:(NumSess*2)
    if k ==1
    alltimes(k,:)=[0 extimes(k,1)];
    elseif any(k==[2,4,6,8,10])
    alltimes(k,:)=[extimes(k/2,1); extimes(k/2,2)]; 
    elseif any(k==[3,5,7,9])
    alltimes(k,:)=[extimes((k-1)/2,2); extimes((k+1)/2,1)];     
    end
end

save('alltimes.mat','alltimes');
Ratn=Animals(ii);
CSCn={{49},{17},{23},{115},{17},{49},{17},{23},{115},{17},{49},{17},{23},{115},{17},{49},{17},{23},{115},{17},{49},{17},{23},{115},{17}};
% CSCn={{54},{17},{3},{25},{4},{54},{17},{2},{25},{4},{54},{17},{2},{25},{4}};
colors = {[0.83 0.42 0.4],[0.52 0.45 0.83],[0.68 0.88 0.43],[0.5 0.73 0.93],...
    [0.75 0.42 0.92],[0.94 0.64 0.94],[0.94 0.94 0.52],[0.52 0.94 0.94]}; % colors used for plotting


nplot = 40; % plot LFP examples with length of ddt second  
ddt = 2; % in the unit of second
justplot = 0;


sesall = 1:1:size(alltimes,1); % sleep session in alltimes
nsess = length(sesall);
        P = files{ii};
        S = dir(fullfile(P,'*.ncs'));
        N = {S.name};
        LFPData=[];
if ~justplot
    for is = 1:nsess
       if is==1;
        tst = alltimes(sesall(is),1)+5;
       else
        tst = alltimes(sesall(is),1);  
       end
     
        tend = alltimes(sesall(is),2);
      
        
     
        disp(['Processing Session ',num2str(is),'...'])
%% determine the CSC file to read for the best tetrode
       for ic = 1:length(CSCn{1,ii})
            cscname = ['CSC',num2str(CSCn{1,ii}{1,ic}),'.ncs'];
            X = strcmp(N,cscname);
            % %         X = ~cellfun('isempty',strfind(N,num2str(cscn(ic))));
            % %         I = imread(fullfile(P,N{X}));
            if sum(X)>0
                ncsntemp = N{X};
                % %                 [eeg,ts,fre] = ReadCSC_yuchen(strcat(InPath,ncsntemp),'time', tst, tend);
                % %                 Ptet(it,is).LFP = eeg;
                % %                 Ptet(it,is).t = ts;
                % %                 Ptet(it,is).fs = fre;
                load starttime.mat
                tss=ts_1;
                tst1=(tst*10000+tss)*100; tend1=(tend*10000+tss)*100;
%                 if is==nsess;
%                    [event_time2, event_flag2] = GetEventTime;
%                    tend1=(event_time2(end)-(5*10^6));
%                 end
                            
               
                fromInd=tst1;
                toInd=tend1;
                ddrate=1/30000;
               
                if ii==9 && is > 5
                filename= strcat('CSC115','_1','.ncs');  
                [timestamps,freq,eeg] = getRawCSCData(filename, fromInd, toInd);
                ts1=timestamps(1);
                ts1=ts1/100-tss;
                ts=ts1/10000;
                fre0 = 3e4;  % original sampling frequency;
                fs = 2000;  % target sampling frequency;
                drate = round(fre0/fs);
                for n=2:length(eeg)
                ts(n,1)=ts(n-1,1)+ddrate;
                end
                eeg=eeg(1:drate:end);
                ts=ts(1:drate:end);
                elseif  ii==9 && is == 5
                filename= ncsntemp;               
                [timestamps1,freq,eeg1] = getRawCSCData(filename, fromInd, toInd) ;   
                filename= strcat('CSC115','_1','.ncs'); 
                [timestamps2,freq,eeg2] = getRawCSCData(filename, fromInd, toInd);   
                    
                ts11=timestamps1(1);
                ts11=ts11/100-tss;
                ts111=ts11/10000;
                fre0 = 3e4;  % original sampling frequency;
                fs = 2000;  % target sampling frequency;
                drate = round(fre0/fs);
                for n=2:length(eeg1)
                ts111(n,1)=ts111(n-1,1)+ddrate;
                end
                eeg111=eeg1(1:drate:end);
                ts111=ts111(1:drate:end);
                
                ts12=timestamps2(1);
                ts12=ts12/100-tss;
                ts112=ts12/10000;
                fre0 = 3e4;  % original sampling frequency;
                fs = 2000;  % target sampling frequency;
                drate = round(fre0/fs);
                for n=2:length(eeg2)
                ts112(n,1)=ts112(n-1,1)+ddrate;
                end
                eeg112=eeg2(1:drate:end);
                ts112=ts112(1:drate:end);
                
                ts=[ts111;ts112];eeg=[eeg111;eeg112];
                
                elseif ii==6 &  is == 2
                filename= ncsntemp;               
                [timestamps,freq,eeg] = getRawCSCData_Mod(filename, fromInd, toInd) ;    
                    
               
                ts1=timestamps/100-tss;
                ts=ts1/10000;                          
                fre0 = 3e4;  % original sampling frequency;
                fs = 2000;  % target sampling frequency;
                drate = round(fre0/fs);
%                 for n=2:length(eeg)
%                 ts(n,1)=ts(n-1,1)+ddrate;
%                 end
                eeg=eeg(1:drate:end);
                ts=ts(1:drate:end);                
%               
                else
                filename= ncsntemp;               
                [timestamps,freq,eeg] = getRawCSCData(filename, fromInd, toInd) ;    
                    
                ts1=timestamps(1);
                ts1=ts1/100-tss;
                ts=ts1/10000;
                fre0 = 3e4;  % original sampling frequency;
                fs = 2000;  % target sampling frequency;
                drate = round(fre0/fs);
                for n=2:length(eeg)
                ts(n,1)=ts(n-1,1)+ddrate;
                end
                eeg=eeg(1:drate:end);
                ts=ts(1:drate:end);
                end
%                 ts = ts/10^4; % now it has the unit of seconds
                LFPData(is).LFP = eeg;
                LFPData(is).t = ts;
                LFPData(is).fs = fs;
             
                % get the Delta band
                deltaband = [1 4];
                MyFilt=fir1(8000,deltaband/(fs/2));
                Filtered = Filter0(MyFilt,eeg);
                LFPData(is).firLFP{1} = Filtered;
                LFPData(is).firbands{1} = deltaband;
                
                % get the theta band
                thetaband = [6 10];
                MyFilt=fir1(8000,thetaband/(fs/2));
                Filtered = Filter0(MyFilt,eeg);
                LFPData(is).firLFP{2} = Filtered;
                LFPData(is).firbands{2} = thetaband;
                
                % get the theta band
                Betaband = [15 35];
                MyFilt=fir1(8000,Betaband/(fs/2));
                Filtered = Filter0(MyFilt,eeg);
                LFPData(is).firLFP{5} = Filtered;
                LFPData(is).firbands{5} = Betaband;
                
                % get the ripple band power
                ripband = [140 250];
                MyFilt=fir1(8000,ripband/(fs/2));
                Filtered = Filter0(MyFilt,eeg);
                LFPData(is).firLFP{3} = Filtered;
                LFPData(is).firbands{3} = ripband;
                
                % get the high speed theta band
                htband = [6 12];
                MyFilt=fir1(8000,htband/(fs/2));
                Filtered = Filter0(MyFilt,eeg);
                LFPData(is).firLFP{4} = Filtered;
                LFPData(is).firbands{4} = htband;
                LFPData(is).name=cscname;
                
               
                break
            end
        end
    end
    save('LFPDataV2.mat','LFPData','-v7.3');
else
    load ('LFPDataV2.mat','LFPData')
    nsess = length(LFPData);
end

% for is = 1:nsess
%     tst = alltimes(sesall(is),1);
%     tend = alltimes(sesall(is),2);
%     %% plot the LFP samples
%     tstori = tst;
%     tendori = tend;
%     dt = ceil((tendori-tstori)/nplot);
%     tplot = tstori:dt:tendori;
%     sespath = [files{ii},'/ReplayData',filesep];
%     if ~exist(sespath, 'dir')
%        mkdir(sespath)
%     end
%    
%     for ip = 1:length(tplot)
%         figure()
%         hold on
%         set(gcf,'outerposition',get(0,'screensize'));
%         trange = [tplot(ip),tplot(ip)+ddt];
%         if trange(2)>tend
%             break
%         end
%     
%         ttname = ['CSC',num2str(CSCn{1,ii}{1,1})];
%         try
%             ts = LFPData(is).t;
%             eeg = LFPData(is).LFP;
%             tempind = ts >= trange(1) & ts <= trange(2);
%             if sum(tempind) == 0
%                 continue
%             end
%             tempt = ts(tempind);
%             tempeeg = eeg(tempind);
%             plot(tempt,tempeeg,'color',colors{1+1},'linewidth',0.5)
%             plot([tempt(1) tempt(end)],[0 0],'k--','linewidth',0.5)
%         end
% 
%         ytt = 0;
%         yticks(ytt)
%         yticklabels(ttname)
%         clear ttname
%         tt = ['Segment',num2str(ip)];
%         title(tt)
%         axis tight
% %         File_Path = strcat(sespath,'sleep',num2str(is),tt,'.fig');
% %         saveas(gcf, File_Path);
%         File_Path = strcat(sespath,'sleep',num2str(is),tt,'.jpg');
%         saveas(gcf, File_Path);
%         close all
%     end
%     
% end

disp(strcat('CA_Sleep_loadPyrLFP_Bestchanel: Done','forFolder',num2str(ii)));
end