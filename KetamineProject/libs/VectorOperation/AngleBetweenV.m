function ThetaInDegrees = AngleBetweenV(u,v,dim)
% compute the angle in degree between population vectors based on inner
% product along dimension dim

if nargin < 3
    CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    ThetaInDegrees = real(acosd(CosTheta));
else
    dotprod = dot(u,v,dim);
    unorm = vecnorm(u,2,dim);
    vnorm = vecnorm(v,2,dim);
    CosTheta = max(min(dotprod./(unorm.*vnorm),1),-1);
    ThetaInDegrees = real(acosd(CosTheta));
end

end