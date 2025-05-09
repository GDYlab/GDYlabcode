function a = becolumn(a)

% if isempty(a)
%     a = a;
%     return
% end
if ~isvector(a)
%    warning('input is not a vector, taking all the segments') 
   a = a(:);
end


if isrow(a)
    a = a';
end
if iscolumn(a)
    a = a;
end

end

