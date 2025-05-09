function detses = Det_FindDetTSes(tnum,ses)
% this function give the detour session of the detour track
detses = nan;
for is = 1:length(ses)
    if ~ismember(tnum,ses(is).tra_p)
        detses = is;
        return
    end
end
end