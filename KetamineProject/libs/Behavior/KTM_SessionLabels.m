function labels = KTM_SessionLabels(nsess)
% KTM_SessionLabels
% Correct data-order session labels for the ketamine Run2-split data.
% 4 sessions: {Run1, Run2Bi, Run3, Run2woBi}; 3 sessions: {Run1, Run2Bi, Run3}.
% (Run2woBi follows Run3 in data order; the manuscript figures re-order it
%  before Run3 for display only.)
switch nsess
    case 4
        labels = {'Run1','Run2Bi','Run3','Run2woBi'};
    case 3
        labels = {'Run1','Run2Bi','Run3'};
    otherwise
        labels = arrayfun(@(k) sprintf('Run%d',k),1:nsess,'uni',0);
end
end
