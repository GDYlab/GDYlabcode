function prepost = Det_SleepBehPrePost(sleeps,behs)
% this function determine whether the sleep session is pre or post the
% behaviour session in detour experiment
% sleeps is the sleep session
% behs is the behavior session

%% check input
if ~ismember(sleeps,[1,2,3,4]) || ~ismember(behs,[1,2,3,4])
error('Invalid sleep or behaviour session number')
end

if sleeps == 1
    prepost = 'pre';
    return
end

if sleeps == 4
    prepost = 'post';
    return
end

if sleeps == 2
    if behs == 1
        prepost = 'post';
    else
        prepost = 'pre';
    end
    return
end

if sleeps == 3
    if behs == 4
        prepost = 'pre';
    else
        prepost = 'post';
    end
    return
end

end