function  sframes = BayePosDecode_LabCode_Lite_CosineExp(trainplf,trainlambda,trainmesh,ClF,framet,...
    lposbtrain,tbin,ovlap,meanact,alllambda)
%  traindirplf
%  This function run bayesian decoder for frames during awake rest within a given session, given track
%  it can also be used to decode theta sequence or sleep frames.
%  decoding based on angle difference and exponential distribution

%  trainplf     is the place fields firing rate with row being cells and column
%               being spatial meshgrid
%  trainlambda  is the parameter of exponential distribution for each
%               position bins
%  trainmesh    is the spatial meshgrid of the training firing rate
%  ClF          is a cell of firing rate matrix with length equal to the
%               number of cells. With in cell it has standard format
%  lposbtrain   is the beginning and end of the linear position for the region of interets in training session (in plfields)
%  framet       is the time range for decoding
%  tbin         temporal bin
%  ovlap        overlap (0-1) of time window
%  meanact      is the mean population vector of decoding state(average over frames)
%  alllambda    is the parameter of exponential distribution of angle
%               difference with meanact
%  OutPut       sframes  decoding results 

addnoise = 1;
%% refine the traning set to the track
posind = trainmesh>=lposbtrain(1) & trainmesh<lposbtrain(2);
trainplf = trainplf(:,posind);
trainlambda = trainlambda(posind);
tramesh = trainmesh(posind);

%% add some noise to reduce undecodable moments
if addnoise
    trainplf = addnoisetoplf(trainplf);
end
%% start bayesian decoding
% frame loop

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
    
    [dprob, normpdf, tmesh, nact, nspk] = reconstruct_CosineExp(framet(ii,1), framet(ii,2), cells,...
        trainplf, meanact, trainlambda, alllambda,'tau',tbin,'percent_overlap',ovlap);

%     stdp = std(dprob,0,1); % in dprob D1 is space, D2 is time
%     pind = stdp<=1e-12; % uniform probability distribution
%     dprob(:,pind) = nan;
    sframes(ii).pdf= dprob;
    sframes(ii).normpdf= normpdf;
    sframes(ii).tbin= tmesh;
    sframes(ii).spacebin= tramesh;
    sframes(ii).nact = nact;
    sframes(ii).nspk = nspk;
end


end

