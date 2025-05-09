function len = Detour_GetDetourSegLen(ses,dettrack)
% get detour track first and last segment lengths
len = cell(1,4); % for all the tracks, only detour track will have value
nsess = length(ses);
tralinear = [1,2,3,4]; % track number of original linear tracks;
for it = dettrack
    xs = NaN;
    for is = 1:nsess
        if ~ismember(it,ses(is).tra_p)
            xs = is; break      % xs is the session with detour track
        end
    end
    newtrack = setdiff(ses(xs).tra_p,tralinear);
    xtraind = find(ses(xs).tra_p == newtrack);
    % get segments length:
    markpos = sort(ses(xs).tramarks(xtraind).posl,'ascend');
    len{it} = [markpos(1)-ses(xs).tralim(xtraind,1),ses(xs).tralim(xtraind,2)-markpos(4)];  % first segment and last segment
end


end

