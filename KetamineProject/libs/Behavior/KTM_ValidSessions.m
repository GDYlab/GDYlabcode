function [validses,sesstr] = KTM_ValidSessions(sessioncell)
% KTM_ValidSessions
% Returns the non-empty session slots in a sample dataset and their labels,
% so the tutorial adapts to however many effective run sessions an animal has.
%
% The Run2-split data always stores up to 4 session slots whose meaning is
% fixed by index: 1=Run1, 2=Run2Bi, 3=Run3, 4=Run2woBi. Animals with only 3
% effective sessions (e.g. URat1, URat3) leave slot 4 (Run2woBi) empty.
%
% input:  sessioncell - a per-session cell array (e.g. mazel); a slot is
%                       "valid" if its element is non-empty
% outputs: validses   - row vector of valid (non-empty) slot indices
%          sesstr     - 1 x numel(validses) labels for those slots
allnames = {'Run1','Run2Bi','Run3','Run2woBi'};
nslot = numel(sessioncell);
validses = find(~cellfun(@isempty,sessioncell(:)'));
names = allnames(1:min(nslot,numel(allnames)));
sesstr = names(validses);
end
