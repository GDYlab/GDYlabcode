function  sframes = BayePosDecode_LabCode_Lite(traindirplf,trainmesh,ClF,framet,...
    lposbtrain,tbin,ovlap)
%  BayePosDecode_LabCode_Lite
%  This function run bayesian decoder given several time windows to decode
%  traindirplf  is the training place fields firing rate with row being cells and column
%               being spatial meshgrid
%  trainmesh    is the spatial meshgrid of the training firing rate
%  ClF          is a cell of firing rate matrix with length equal to the
%               number of cells. With in cell it has standard format
%  lposbtrain   is the beginning and end of the linear position for the region of interets in training session (in plfields)
%  framet       is the time range for decoding
%  tbin         temporal bin
%  ovlap        overlap (0-1) of time window
%  OutPut       sframes  decoding results 
     
%% check input
if nargin < 7
    ovlap = 0;
end

%% refine the placemap to the linear position range of interest
posind = trainmesh>=lposbtrain(1) & trainmesh<lposbtrain(2);
traindirplf = traindirplf(:,posind);
tramesh = trainmesh(posind);

%% add some noise to reduce undecodable moments
for ic = 1:size(traindirplf,1)
    plftemp = traindirplf(ic,:);
    if max(plftemp) == 0
        traindirplf(ic,:) = ones(size(plftemp))*1e-8;
    else
        traindirplf(ic,:) = traindirplf(ic,:)+nanmean(traindirplf(ic,:))/1e5;
    end
end
%% start bayesian decoding
% loop over time windows to decode
sframes = struct;
for ii = 1:size(framet,1)
    cells = struct;
    for ic = 1:length(ClF)
        % refine spking events to the frame time and track position
        spkt = ClF{ic}(:,1);
        laptraind = spkt>= framet(ii,1) & spkt<= framet(ii,2);
        % %         dc_spks{ic} = spkt(laptraind);
        cells(ic).time = spkt(laptraind);
    end
    % run decoding
    [dprob, tmesh,nact] = reconstruct_yuchen(framet(ii,1), framet(ii,2), traindirplf', cells, 'tau',tbin,'percent_overlap',ovlap);
    stdp = std(dprob,0,1); % in dprob D1 is space, D2 is time
    pind = stdp<=1e-12; % uniform probability distribution
    dprob(:,pind) = nan;
    sframes(ii).pdf= dprob;
    sframes(ii).tbin= tmesh;
    sframes(ii).spacebin= tramesh;
    sframes(ii).nact = nact;
end


end

