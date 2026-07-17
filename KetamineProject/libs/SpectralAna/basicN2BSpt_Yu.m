function [b, b_ori,a12,aa12, frq] = basicN2BSpt_Yu (y,  nfft, wind, nsamp, overlap)
% Coding: Alex Sheremet 2014
%BICOHER - Direct (FD) method for estimating bicoherence
%	[bsp,frq] = bicoher (y,  nfft, wind, segsamp, overlap)
%	y     - data vector or time-series
%	nfft - fft length [default = power of two > nsamp]
%	       actual size used is power of two greater than 'nsamp'
%	wind - specifies the time-domain window to be applied to each
%	       data segment; should be of length 'segsamp' (see below);
%		otherwise, the default Hanning window is used.
%	nsamp - samples per segment [default: such that we have 8 segments]
%	        - if x is a matrix, segsamp is set to the number of rows
%	overlap - percentage overlap, allowed range [0,99]. [default = 50];
%	        - if x is a matrix, overlap is set to 0.
%	b       - estimated normalized bispectrum: an nfft x nfft array, with origin
%	          at the center, and axes pointing down and to the right.
%	bs      - estimated std of normalized bispectrum: an nfft x nfft array, with origin
%	          at the center, and axes pointing down and to the right.
%	frq   - vector of frequencies associated with the rows and columns
%	          of bsp;  sampling frequency is assumed to be 1.

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.7 $
%  A. Swami   January 20, 1995

%     RESTRICTED RIGHTS LEGEND
% Use, duplication, or disclosure by the Government is subject to
% restrictions as set forth in subparagraph (c) (1) (ii) of the
% Rights in Technical Data and Computer Software clause of DFARS
% 252.227-7013.
% Manufacturer: United Signals & Systems, Inc., P.O. Box 2374,
% Culver City, California 90231.
%
%  This material may be reproduced by or for the U.S. Government pursuant
%  to the copyright license under the clause at DFARS 252.227-7013.

% based on the bicoherence code in HOSA, modified by Alex to return std

% --------------------- parameter checks -----------------------------

	[ly, nrecs] = size(y);
	if (ly == 1) y = y(:);  ly = nrecs; nrecs = 1; end

	if (exist('nfft') ~= 1)            nfft = 128; end
	if (exist('overlap') ~= 1)      overlap = 50;  end
	overlap = max(0,min(overlap,99));
	if (nrecs > 1)                  overlap = 0;   end
	if (exist('nsamp') ~= 1)          nsamp = 0;  end
	if (nrecs > 1)                    nsamp = ly;  end

	if (nrecs == 1 & nsamp <= 0)
		nsamp = fix(ly/ (8 - 7 * overlap/100));
	end
	if (nfft  < nsamp)   nfft = 2^nextpow2(nsamp); end

	overlap  = fix( nsamp * overlap/100);
	dKSeg = nsamp - overlap;
	nrecs    = fix ( (ly*nrecs - overlap) / dKSeg);

% ----------------------------------------------------------------------
	if (exist('wind') ~= 1) wind = hanning(nsamp); end
	[rw,cw] = size(wind);
	if (min(rw,cw) ~= 1 | max(rw,cw) ~= nsamp)
		disp(['Segment size  is ',int2str(nsamp)])
		disp(['"wind" array  is ',int2str(rw),' by ',int2str(cw)])
		disp(['Using default Hanning window'])
		wind = hanning(nsamp);
	end
	wind = wind(:);
% ---------------- accumulate triple products ----------------------


	mask = hankel([1:nfft],[nfft,1:nfft-1] );   % the hankel mask (faster)
	ind  = [1:nsamp];
	b=zeros(nfft,nfft); aa12=zeros(nfft,nfft); a12=zeros(nfft,nfft);
	
	
% modified by Alex for normalization type 2
% B=<A(f1)A(f2)conj(A(f1+f2))> / sqrt(<|A(f1)A(f1)|^2><|A(f1)|^2> )
% for 0<= B <=1
	for k = 1:nrecs
		ys=y(ind);
		ys=(ys(:)-mean(ys)) .* wind;
		Yf=fft(ys,nfft)/nsamp;
        
		aa12=aa12+abs(Yf*Yf.').^2;
        tmp = Yf(mask);
		a12=a12+abs(conj(tmp)).^2;
		b=b+(Yf*Yf.').*conj(tmp);
		ind=ind+dKSeg;
	end
	aa12=aa12/nrecs; a12=a12/nrecs;
    b_ori = b/nrecs;   b_ori = fftshift(b_ori);  % Yu
	b=b/nrecs./sqrt(a12.*aa12);
	b=fftshift(b);
    a12 = fftshift(a12);    aa12 = fftshift(aa12);   % Yu
    

	if (rem(nfft,2) == 0)
		frq = [-nfft/2:(nfft/2-1)]'/nfft;
	else
		frq = [-(nfft-1)/2:(nfft-1)/2]'/nfft;
	end
