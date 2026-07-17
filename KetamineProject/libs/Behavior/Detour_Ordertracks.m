function ses = Detour_Ordertracks(ses)
% order tracks by linear position for each session
nsess = length(ses);
for is = 1:nsess
    trac = mean(ses(is).tralim,2);
    [~,traodr] = sort(trac,'ascend');
    ses(is).tra_p = ses(is).tra_p(traodr);
    ses(is).tralim = ses(is).tralim(traodr,:);
    if isfield(ses,'traxlim')
        ses(is).traxlim = ses(is).traxlim(traodr,:);
    end
    if isfield(ses,'traylim')
        ses(is).traylim = ses(is).traylim(traodr,:);
    end
    if isfield(ses,'tra_plfd')
        ses(is).tra_plfd = ses(is).tra_plfd(traodr);
    end
    if isfield(ses,'tramarks')
        ses(is).tramarks = ses(is).tramarks(traodr);
    end
end
end

