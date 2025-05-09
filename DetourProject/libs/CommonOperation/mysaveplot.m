function mysaveplot(OutPath,stt,fhandle,onlyjpg,overwrite)
if nargin < 4
    onlyjpg = 0;
end
if nargin < 5
    overwrite = 1;
end


if ~onlyjpg
    File_Path = strcat(OutPath,stt,'.fig');
    saveas(fhandle, File_Path);
end

% Check if file exists, and append "_v2", "_v3", etc., if necessary
File_Path = strcat(OutPath,stt,'.jpg');
if overwrite
    saveas(fhandle, File_Path);
else
    version = 1;
    while exist(File_Path, 'file')
        version = version + 1;
        File_Path = fullfile(OutPath, sprintf('%s_v%d.jpg', stt, version));
    end
    saveas(fhandle, File_Path);
end

end