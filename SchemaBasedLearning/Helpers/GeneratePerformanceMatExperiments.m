cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('Paths')
GenPath=strcat(paths{2},'Results/PlotFigures/Database');
cd(GenPath)
load Rat1Choices
ChoicesRat1=ChoicesRat1(:,3:end);

performance1=(cellfun(@length,ChoicesRat1)-1)
performance=[];
performance_new=[];
performance_old=[];
for i=1:9;
    if i < 5
      for j=1:18
          rr=ChoicesRat1{j,i+3}
          if rr(end)==4 | rr(end)==2 | rr(end)==7
          performance_new{j,i}=rr;
          end
              
      end
       for j=1:18
          rr=ChoicesRat1{j,i+3}
          if rr(end)==1 | rr(end)==5 | rr(end)==8
          performance_old{j,i}=rr;
          end
      end
    else
      for j=1:18
          rr=ChoicesRat1{j,i+3}
          if rr(end)==4 | rr(end)==2 | rr(end)==7
          performance_old{j,i}=rr;
          end
              
      end
       for j=1:18
          rr=ChoicesRat1{j,i+3}
          if rr(end)==1 | rr(end)==5 | rr(end)==8
          performance_new{j,i}=rr;
          end
      end
    end
end
        
    
performance_1{1}=performance1(:,1:3);
performance_1{2}=performance1(:,4:12);
performance_1{3}=performance1(:,13:15);
performance_1{4}=performance1(:,16:19);
% performance{5}=performance1(:,25:30);
% performance{6}=performance1(:,15:17);

per_new1{1}=(cellfun(@length,performance_new)-1);
per_old1{1}=(cellfun(@length,performance_old)-1);

% cd ('C:\Users\babur\Documents\MATLAB\Codes\Behavior')
load Rat2Choices


performance1=(cellfun(@length,ChoicesRat2)-1)
performance1(performance1==-1)=NaN;
performance=[];
performance_new=[];
performance_old=[];
for i=1:8;
    if i < 5
      for j=1:18
          if ~isempty(ChoicesRat2{j,i+3})
          rr=ChoicesRat2{j,i+5}
          if rr(end)==4 | rr(end)==2 | rr(end)==7
          performance_new{j,i}=rr;
          end
          else
          performance_new{j,i}=NaN;  
          end
              
      end
       for j=1:18
          if ~isempty(ChoicesRat2{j,i+3}) 
          rr=ChoicesRat2{j,i+5}
          if rr(end)==1 | rr(end)==5 | rr(end)==8
          performance_old{j,i}=rr;
          end
          else
          performance_old{j,i}=NaN;  
          end    
      end
    else
      for j=1:18
          rr=ChoicesRat2{j,i+3}
          if rr(end)==4 | rr(end)==2 | rr(end)==7
          performance_old{j,i}=rr;
          end
              
      end
       for j=1:18
          rr=ChoicesRat2{j,i+3}
          if rr(end)==1 | rr(end)==5 | rr(end)==8
          performance_new{j,i}=rr;
          end
      end
    end
end
        
    
performance_2{1}=performance1(:,1:3);
performance_2{2}=performance1(:,4:11);
performance_2{3}=performance1(:,12:14);
performance_2{4}=performance1(:,15:18);
% performance_2{5}=performance1(:,25:30);
% performance_2{6}=performance1(:,15:17);

per_new1{2}=(cellfun(@length,performance_new)-1);
per_old1{2}=(cellfun(@length,performance_old)-1);

load Rat3Choices


performance1=(cellfun(@length,ChoicesRat3)-1)
performance1(performance1==-1)=NaN;
performance=[];
performance_new=[];
performance_old=[];
for i=1:6;
    if i < 4
      for j=1:18
          rr=ChoicesRat3{j,i+3}
          if rr(end)==4 | rr(end)==2 | rr(end)==7
          performance_new{j,i}=rr;
          end
              
      end
       for j=1:18
          rr=ChoicesRat3{j,i+3}
          if rr(end)==1 | rr(end)==5 | rr(end)==8
          performance_old{j,i}=rr;
          end
      end
    else
      for j=1:18
          rr=ChoicesRat3{j,i+3}
          if rr(end)==4 | rr(end)==2 | rr(end)==7
          performance_old{j,i}=rr;
          end
              
      end
       for j=1:18
          rr=ChoicesRat3{j,i+3}
          if rr(end)==1 | rr(end)==5 | rr(end)==8
          performance_new{j,i}=rr;
          end
      end
    end
end
        
    
performance_3{1}=performance1(:,1:3);
performance_3{2}=performance1(:,4:9);
performance_3{3}=performance1(:,10:12);
performance_3{4}=performance1(:,13:15);
% performance{5}=performance1(:,25:30);
% performance{6}=performance1(:,15:17);

per_new1{3}=(cellfun(@length,performance_new)-1);
per_old1{3}=(cellfun(@length,performance_old)-1);


load Rat4Choices


performance5=(cellfun(@length,ChoicesRat4)-1)
performance=[];
performance_new=[];
performance_old=[];
for i=1:8;
    if i < 5
      for j=1:18
          rr=ChoicesRat4{j,i+3}
          if rr(end)==4 | rr(end)==2 | rr(end)==7
          performance_new{j,i}=rr;
          end
              
      end
       for j=1:18
          rr=ChoicesRat4{j,i+3}
          if rr(end)==1 | rr(end)==5 | rr(end)==8
          performance_old{j,i}=rr;
          end
      end
    else
      for j=1:18
          rr=ChoicesRat4{j,i+3}
          if rr(end)==4 | rr(end)==2 | rr(end)==7
          performance_old{j,i}=rr;
          end
              
      end
       for j=1:18
          rr=ChoicesRat4{j,i+3}
          if rr(end)==1 | rr(end)==5 | rr(end)==8
          performance_new{j,i}=rr;
          end
      end
    end
end
        
    
performance_5{1}=performance5(:,1:3);
performance_5{2}=performance5(:,4:11);
performance_5{3}=performance5(:,12:14);
performance_5{4}=performance5(:,15:18);
% performance{5}=performance1(:,25:30);
% performance{6}=performance1(:,15:17);

per_new1{4}=(cellfun(@length,performance_new)-1);
per_old1{4}=(cellfun(@length,performance_old)-1);
% % 
% % performance{2}(:,32)=ones(18,1).*NaN; 
% % % performance{3}(:,31)=ones(18,1).*NaN; 
% % performance{4}(:,32)=ones(18,1).*NaN; 
% performance12=[];
% for i=1:5
%     performance12=[performance12; (performance{i})]
% end
% % performance13=[];
% % for i=1:5
% %     performance13=[performance13; (performance{i})]
% % end


%% calculating the day by day error
% cd ('C:\Users\babur\Documents\MATLAB\Codes\Behavior')
load Rat5Choices


performance7=(cellfun(@length,ChoicesRat5)-1)
performance=[];
performance_new=[];
performance_old=[];
for i=1:7;
    if i < 5
      for j=1:18
          rr=ChoicesRat5{j,i+3}
          if rr(end)==4 | rr(end)==2 | rr(end)==7
          performance_new{j,i}=rr;
          end
              
      end
       for j=1:18
          rr=ChoicesRat5{j,i+3}
          if rr(end)==1 | rr(end)==5 | rr(end)==8
          performance_old{j,i}=rr;
          end
      end
    else
      for j=1:18
          rr=ChoicesRat5{j,i+3}
          if rr(end)==4 | rr(end)==2 | rr(end)==7
          performance_old{j,i}=rr;
          end
              
      end
       for j=1:18
          rr=ChoicesRat5{j,i+3}
          if rr(end)==1 | rr(end)==5 | rr(end)==8
          performance_new{j,i}=rr;
          end
      end
    end
end
        
    
performance_7{1}=performance7(:,1:3);
performance_7{2}=performance7(:,4:10);
performance_7{3}=performance7(:,11:13);
performance_7{4}=performance7(:,14:17);
% performance{5}=performance1(:,25:30);
% performance{6}=performance1(:,15:17);

per_new1{5}=(cellfun(@length,performance_new)-1);
per_old1{5}=(cellfun(@length,performance_old)-1);

Err{1}=performance_1; Err{2}=performance_2; Err{3}=performance_3; Err{4}=performance_5; Err{5}=performance_7;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('Paths')
GenPath=strcat(paths{2},'Results/PlotFigures/Database');
cd(GenPath)
save('BehavError.mat','Err','per_new1','per_old1')
