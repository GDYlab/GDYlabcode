function lfp = function_LFPSpec_RunBispectrum(seglfp,fs,thetarange,frange)
% function_LFPSpec_RunBispectrum
% Quantifies the non-sinusoidal (asymmetric) shape of the hippocampal theta
% rhythm during fast run via the LFP bispectrum (manuscript Fig. 5). For each
% session it:
%   (1) estimates the normalized bispectrum (bicoherence) of the best-tetrode
%       LFP, averaged over high-velocity run segments (segment-weighted by
%       their duration; see basicN2BSpt_Yu),
%   (2) extracts the theta wave asymmetry as the mean imaginary part of the
%       bicoherence over the theta-theta interaction region.
% A non-zero imaginary bicoherence in the theta band reflects a front/back
% asymmetric (saw-tooth) theta wave; ketamine (Run2Bi) changes this asymmetry.
%
% inputs:
%   seglfp     : 1 x nsess cell; seglfp{is} is a 1 x nseg cell of LFP traces
%                (row vectors), one per high-velocity run segment of session is
%                (empty for absent sessions). See ../DataStructure.md.
%   fs         : LFP sampling rate, Hz (default 2000)
%   thetarange : [lo hi] theta band in Hz for the asymmetry measure (default [4 10])
%   frange     : [lo hi] frequency range kept in the returned bicoherence map,
%                Hz (default [0 120])
%
% output: lfp struct with fields
%   lfp.sesstr     : 1 x nvalid session labels (effective sessions only)
%   lfp.ff         : frequency axis (Hz) of the returned bicoherence maps
%   lfp.bic        : 1 x nvalid cell; complex bicoherence map (principal domain,
%                    masked to frange). abs() = bicoherence, imag() = asymmetry.
%   lfp.asym       : 1 x nvalid theta wave asymmetry (mean imag over theta region)
%   lfp.dof        : 1 x nvalid degrees of freedom (total run seconds) per session
%   lfp.thetarange : the theta band used
%
% Yuchen Zhou, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com

if nargin < 2 || isempty(fs);         fs = 2000; end
if nargin < 3 || isempty(thetarange); thetarange = [4 10]; end
if nargin < 4 || isempty(frange);     frange = [0 120]; end

nsess  = numel(seglfp);
nseg   = fix(fs*0.8);      % 0.8 s analysis window (segment), as in the manuscript
dofcap = 400;              % stop accumulating once a session has this many seconds

% adapt to the number of effective (non-empty) sessions this animal has
[validses,lfp.sesstr] = KTM_ValidSessions(seglfp);
slotnames = {'Run1','Run2Bi','Run3','Run2woBi'};
bic_all  = cell(1,nsess);
asym_all = nan(1,nsess);
dof_all  = nan(1,nsess);
ff = [];

fprintf('function_LFPSpec_RunBispectrum: %d sessions...\n',nsess);
for is = 1:nsess
    segs = seglfp{is};
    if isempty(segs); continue; end
    fprintf('  session %d/%d (%s): %d run segments\n',is,nsess,slotnames{min(is,4)},numel(segs));

    % ---- (1) duration-weighted average bispectrum over run segments ----
    allb = []; alla1 = []; alla2 = []; dofall = 0;
    for j = 1:numel(segs)
        seg = double(segs{j}(:));
        if numel(seg) < nseg; continue; end       % too short for one window
        [~,b,a1,a2,frq] = basicN2BSpt_Yu(seg,nseg,hanning(nseg),nseg,50);
        dof = round(numel(seg)/fs);
        if isempty(allb)
            allb = b./dof; alla1 = a1./dof; alla2 = a2./dof;
        else
            allb  = (allb.*dofall  + b.*dof) ./(dofall+dof);
            alla1 = (alla1.*dofall + a1.*dof)./(dofall+dof);
            alla2 = (alla2.*dofall + a2.*dof)./(dofall+dof);
        end
        dofall = dofall + dof;
        if dofall > dofcap; break; end
    end
    if isempty(allb); continue; end
    ff = frq * fs;                                  % cycles/sample -> Hz
    bic = allb ./ sqrt(alla1.*alla2);               % normalized bispectrum

    % restrict to frange and keep the principal (non-redundant) domain
    ind = ff >= frange(1) & ff <= frange(2);
    ss = bic(ind,ind);
    ss = triu(ss); ss = fliplr(ss); ss = triu(ss); ss = fliplr(ss);
    bic_all{is} = ss;

    % ---- (2) theta wave asymmetry = mean imag(bicoherence) over theta band ----
    asym_all(is) = local_thetaasym(bic,ff,thetarange);
    dof_all(is)  = dofall;
end

lfp.ff  = ff(ff>=frange(1) & ff<=frange(2));
lfp.ff  = lfp.ff(:)';
lfp.bic  = bic_all(validses);
lfp.asym = asym_all(validses);
lfp.dof  = dof_all(validses);
lfp.thetarange = thetarange;
end

%% ----------------------------------------------------------------------- %%
function asym = local_thetaasym(bic,ff,thetarange)
% mean imaginary part of the bicoherence over the theta-theta triangular region
ind = find(ff <= thetarange(2) & ff >= thetarange(1));
stmp = bic(ind,ind); fftmp = ff(ind);
stmp = triu(stmp); stmp = fliplr(stmp); stmp = triu(stmp); stmp = fliplr(stmp);
px = [thetarange(1) thetarange(2) thetarange(2) thetarange(1)];
py = [thetarange(1) thetarange(1) thetarange(2) thetarange(1)];
asymnow = 0; count = 0;
for i = 1:numel(fftmp)
    for j = 1:numel(fftmp)
        if inpolygon(fftmp(j),fftmp(i),px,py)
            asymnow = asymnow + imag(stmp(i,j));
            count = count + 1;
        end
    end
end
asym = asymnow / count;
end
