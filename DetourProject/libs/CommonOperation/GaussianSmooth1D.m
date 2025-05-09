function smoothedactivity = GaussianSmooth1D(dt,ssd,activity)

Fs = 1/dt;
npoints = round(4.*ssd.*Fs);
kernel = normpdf(linspace(-4*ssd, 4*ssd, 2*npoints+1), 0, ssd);
kernel = kernel./sum(kernel);
smoothedactivity= conv(activity, kernel,'same');

end