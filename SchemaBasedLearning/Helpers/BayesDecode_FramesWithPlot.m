function [dprobp dprob binedges]=BayedDecode_FramesWithPlot(tvec,plf,sp,tbin)
% tvec is binned timestamps of frames
%plf is the place fields, mcell*nfiringratesin each spatial bin
%sp is the vector with rows of timestamp,cellID 

binns=nanmin(tvec):tbin:nanmax(tvec);
tempmat=zeros(size(plf,1),length(binns)-1);
for i = 1:length(binns)-1
    for icell = 1:size(tempmat,1)
        ind = sp(sp(:,2)==icell,1)>=binns(i) &  sp(sp(:,2)==icell,1)<binns(i+1); % find number of spiking events in the timegrid
        tempmat(icell,i) = sum(ind);  % this is the sp_i in the equation. In sp rows are cells, columns are time stamps
    end
end

dprob=parameter_estimation_simple(tbin, plf', tempmat);
dprob(:,all(dprob))=NaN;

binedges=(binns(1:end-1)+binns(2:end))/2;
% imagesc(dprob)
% colormap hot

tempmat1=movmean(tempmat,[1 1],2);
dprobp=parameter_estimation_simple(tbin, plf', tempmat1);
dprobp(:,all(dprob))=NaN;

% ts=nanmin(tvec);te=nanmax(tvec);
% tbins = get_tbins(ts, te,tbin , percent_overlap);
% tempmat=zeros(size(plf,1),length(tbins)-1);
% for i = 1:length(tbins)-1
%     for icell = 1:size(tempmat,1)
%         ind = sp(sp(:,2)==icell,1)>=tbins(i,1) &  sp(sp(:,2)==icell,1)<tbins(i,2); % find number of spiking events in the timegrid
%         tempmat(icell,i) = sum(ind);  % this is the sp_i in the equation. In sp rows are cells, columns are time stamps
%     end
% end
% 
% dprobp=parameter_estimation_simple(tbin, plf', tempmat);
% dprobp(:,all(dprobp))=NaN;
% 
% binedgesp=(tbins(1:end-1,1)'+tbins(2:end,2)')/2;
% 
% 
% function tbins = get_tbins(ts, te, tbin, percent_overlap)
%     tbins = ts:tbin*(1-percent_overlap):te;
%     if numel(tbins)==1
%         tbins = [tbins, tbins+tbin];
%     end
%     tbins = [tbins', tbins'+tbin];
% 
