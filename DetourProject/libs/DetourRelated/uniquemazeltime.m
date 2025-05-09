function mazel = uniquemazeltime(mazel)

% this function remove redundant time stamps in mazel
if iscell(mazel)
    for is = 1:length(mazel)
        smazel = mazel{is};
        mazelt = smazel(:,6);
        % remove redundant timestamps due to track overlap
        [~,Iuni,~] = unique(mazelt);
        smazel = smazel(Iuni,:);
        mazel{is} = smazel;
    end
else
    mazelt = mazel(:,6);
    % remove redundant timestamps due to track overlap
    [~,Iuni,~] = unique(mazelt);
    mazel = mazel(Iuni,:);
end

end