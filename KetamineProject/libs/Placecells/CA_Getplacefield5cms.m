function [plf,PlFmesh,OutLall,ClFnew] = CA_Getplacefield5cms(clf,mazel,ns,grid,mintidx,varargin)
% this function get the place field given the clusters and behaviour data in a given session
% mazel is a matrix with columns being: xpos, ypos, track num, linear pos, run dir, time, vel, lap dir
% clf is cell of clusters, in each cell, it is a matrix with columns being spktime, linear pos, vel, lap dir, run dir, posx, posy
% ns is the session number
% plf D1 are clusters, D2 are two directions, D3 are firing rate in sptial mesh
% defined in vector PlFmesh

if nargin < 5
mintidx = 0; % whether or not use minimum occupancy time 
end
if nargin < 4
grid = 2; % grid size in space
end
mintime = 0.1; % min occupancy in seconds for each bin/pixel; originally 0.5

args.vlim = 5; % run vel threshold
args.ssd = 2; % smooth parameter
args = parseArgs(varargin, args);

ncells = size(clf,1);

%% process behaviour data

postime = mazel(:,6);
% remove redundant timestamps due to track overlap
[postime,Iuni,~] = unique(postime);
mazel = mazel(Iuni,:);

% posx, posy, poslinear, velocity, direction/lap, direction/run
posvel = mazel(:,[1 2 4 7 8 5]);
% velocity can be negative take magnitude
posvel(:,4) = abs(posvel(:,4));

ClFnew = cell(ncells,1);
%% cluster loop
for ic = 1:ncells
    tempc = clf{ic,ns};
    cltime = tempc(:,1);
    %% interplate firing events to behaviour recording
    clposvel = interp1 (postime, posvel, cltime);

    ClFnew{ic} = [cltime clposvel(:,[3:6,1,2])];

    cl1 = clposvel(find(clposvel(:,4)>=args.vlim & clposvel(:,5)==-1 & clposvel(:,6)~=0), 1:3); % moves in the opposite (dir 1) dir
    pos1 = posvel(find(posvel(:,4)>=args.vlim & posvel(:,5)==-1 & posvel(:,6)~=0), 1:3); % moves in the opposite (dir 1) dir
    cl2 = clposvel(find(clposvel(:,4)>=args.vlim & clposvel(:,5)==1 & clposvel(:,6)~=0), 1:3); % moves in the correct (dir 2) dir
    pos2 = posvel(find(posvel(:,4)>=args.vlim & posvel(:,5)==1 & posvel(:,6)~=0), 1:3); % moves in the correct (dir 2) dir
    
    %% get firing rate in spatial meshgrid
    % define the meshgrid of place field
    PlFmesh=[0:grid:1200-grid];
    
    [cloutL1] = hist(cl1(:,3),PlFmesh);
    [cloutL2] = hist(cl2(:,3),PlFmesh);
    [outL1] = hist(pos1(:,3),PlFmesh);
    [outL2] = hist(pos2(:,3),PlFmesh);
    if mintidx
    outL1(find(outL1>=0 & outL1<30*mintime)) = 0; % added =
    outL2(find(outL2>=0 & outL2<30*mintime)) = 0; % added =
    end
    ratioL1 = 30*cloutL1./outL1; % firing rate in direction 1 running in linearized trajectory
    ratioL2 = 30*cloutL2./outL2; % firing rate in direction 2 running in linearized trajectory
    
    ratioL1(find(outL1==0)) = 0;
    ratioL2(find(outL2==0)) = 0;
    
    %% smooth
    sd = args.ssd; % smooth parameter
    Fs = 1/(grid);
    npoints = round(4.*sd.*Fs);
    kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
    kernel = kernel./sum(kernel);
    
    % ratio1 = conv2(ratio1, kernel, 'same');
    % ratio2 = conv2(ratio2, kernel, 'same');
    ratioL1 = conv2(ratioL1, kernel, 'same');   % smooth
    ratioL2 = conv2(ratioL2, kernel, 'same');
    
    %% output
    plf(ic,1,:) = ratioL1;  % firing rate in linearized trajectory in direction 1 (-1 in mazel)
    plf(ic,2,:) = ratioL2;  % firing rate in linearized trajectory in direction 2 (1 in mazel)
    
end
OutLall = [outL1;outL2]; 
% OutLall gives the occupancy in case it is needed to remove small
% occupancy spatial bins, occupany is the same across all the clusters

end
