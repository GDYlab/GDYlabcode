function mi = function_RunSpikeMutualInfo(ses,mazel,M2,ClF,tbin,spkcri,tlag)
% function_RunSpikeMutualInfo
% Measures pairwise functional connectivity between cells during fast run as
% the mutual information (MI) between one cell's recent past and another
% cell's present spiking (manuscript Fig. 4). For each session it:
%   (1) time-bins spikes (tbin) within high-velocity run segments,
%   (2) for every ordered cell pair with enough spikes, computes the
%       normalized point mutual information nMI(Xpast; Ynow) at a lag tlag,
%   (3) repeats on a within-segment temporal shuffle as a chance control.
% Ketamine (Run2Bi) is expected to weaken this connectivity (lower MI).
%
% inputs:
%   ses,mazel : track geometry / behavioral trajectory (see ../DataStructure.md)
%   M2        : per-cell parameters; column 21 is cell type (0 = interneuron)
%   ClF       : nc x nsess spikes per cell per session (pyr + int; see DataStructure)
%   tbin      : spike time bin, s (default 0.4)
%   spkcri    : minimum spikes (across run bins) for a cell to be included (default 5)
%   tlag      : memory length in bins for the MI (default 1)
%
% output: mi struct with fields
%   mi.sesstr  : 1 x nsess session labels (data order)
%   mi.lognMI  : 1 x nsess cell; log normalized MI over all valid ordered pairs (data)
%   mi.lognMIshf : 1 x nsess cell; same for the temporal shuffle
%   mi.MImat   : 1 x nsess cell; nc x nc normalized-MI matrix (rows past, cols present)
%   mi.nint    : number of interneurons (cells are ordered interneurons-then-pyramidal)
%
% Yuchen Zhou, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com

if nargin < 5 || isempty(tbin);   tbin = 0.4; end
if nargin < 6 || isempty(spkcri); spkcri = 5; end
if nargin < 7 || isempty(tlag);   tlag = 1; end

%% preprocess
ses   = Detour_Ordertracks(ses);
mazel = uniquemazeltime(mazel);
ClF   = cellfun(@double,ClF,'uni',0);
nsess = size(ClF,2);
dirs  = [-1,1];

% order cells interneurons-then-pyramidal (cell type in M2 col 21, 0 = int)
celltype = M2(:,21,1);
intid = find(celltype==0); pyrid = find(celltype~=0);
cellorder = [intid(:); pyrid(:)];
nint = numel(intid); nc = numel(cellorder);

% run-epoch criteria
vcri = 5; tlencri = 1; discri = 10;

% adapt to the number of effective (non-empty) sessions this animal has
[validses,mi.sesstr] = KTM_ValidSessions(mazel);
slotnames = {'Run1','Run2Bi','Run3','Run2woBi'}; % full slot labels (for progress text)
mi.lognMI = cell(1,nsess);
mi.lognMIshf = cell(1,nsess);
mi.MImat = cell(1,nsess);
mi.nint = nint;

fprintf('function_RunSpikeMutualInfo: %d sessions...\n',nsess);
for is = 1:nsess
    if isempty(mazel{is}); continue; end
    fprintf('  session %d/%d (%s)\n',is,nsess,slotnames{min(is,4)});
    smazel = mazel{is};
    mt = smazel(:,6); mlp = smazel(:,4); mvel = smazel(:,7);

    % ---- (1) time-binned spikes within high-velocity run segments ----
    segspk = cell(1,nc);        % segspk{ic}{iseg} = column of spike counts
    spkcount = zeros(1,nc);
    for idir = 1:2
        runidx = abs(mvel) >= vcri;
        ep = GetEpochwithCrioverLen(mt,runidx,tlencri);
        for il = 1:size(ep,1)
            sel = mt>=ep(il,1) & mt<ep(il,2);
            lp = mlp(sel);
            if (lp(end)-lp(1))*dirs(idir) < discri; continue; end % require real traversal
            edges = ep(il,1):tbin:ep(il,2);
            if numel(edges) < 2; continue; end
            v = abs(interp1(mt,mvel,edges)); v(isnan(v)) = 0;
            keep = v >= vcri;                 % drop bins that dip below threshold
            for jc = 1:nc
                ic = cellorder(jc);
                st = ClF{ic,is}(:,1);
                c = histcounts_row_kf(st, [edges(:), edges(:)+tbin]);
                c = c(:); c = c(keep);
                segspk{jc} = cat(1,segspk{jc},{c});
                spkcount(jc) = spkcount(jc) + sum(c);
            end
        end
    end

    % temporal shuffle: permute counts within each cell across all segments
    shfspk = cell(1,nc);
    for jc = 1:nc; shfspk{jc} = local_shfxseg(segspk{jc}); end

    % ---- (2-3) pairwise normalized MI (data and shuffle) ----
    MImat = nan(nc,nc); lognMI = []; lognMIshf = [];
    for ii = 1:nc
        if spkcount(ii) < spkcri; continue; end
        for jj = 1:nc
            if jj==ii || spkcount(jj) < spkcri; continue; end
            % nMI between cell ii's past and cell jj's present
            [~,~,~,~,~,~,~,~,nXY] = PointMutualInfo_multiseg(segspk{ii},segspk{jj},tlag);
            [~,~,~,~,~,~,~,~,snXY] = PointMutualInfo_multiseg(shfspk{ii},shfspk{jj},tlag);
            MImat(ii,jj) = nXY;
            if nXY>0;  lognMI(end+1,1)    = log(nXY);  end %#ok<AGROW>
            if snXY>0; lognMIshf(end+1,1) = log(snXY); end %#ok<AGROW>
        end
    end
    mi.MImat{is} = MImat;
    mi.lognMI{is} = lognMI;
    mi.lognMIshf{is} = lognMIshf;
end
% keep only the effective (non-empty) sessions, aligned with mi.sesstr
mi.lognMI = mi.lognMI(validses);
mi.lognMIshf = mi.lognMIshf(validses);
mi.MImat = mi.MImat(validses);
end

%% ----------------------------------------------------------------------- %%
function sspk = local_shfxseg(spk)
% within-cell temporal shuffle: pool counts across segments, permute, re-split
nseg = numel(spk);
seglen = cellfun(@numel,spk);
allspk = cat(1,spk{:});
allspk = allspk(randperm(numel(allspk)));
cl = [0; cumsum(seglen(:))];
sspk = cell(nseg,1);
for iseg = 1:nseg; sspk{iseg} = allspk(cl(iseg)+1:cl(iseg+1)); end
end
