% A_fPsfF.m
% 
% Function to encode sparse spatial data into non-sparse Fourier,
% coefficients, with PSF taken into account in Fourier domain.
%
% Input:
% x             sparse function (1D)
% dim           Dimensions of original data ([x, y] vector for standard image)
% picks         indices of known spatial values within full size data
% OTF           optical transfer function, as obtained from makeOTF() function
% 
% Output:
% b             measurement vector
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

function b = A_fPsfF(x, dim, picks, OTF)

    K = length(picks);
    if (size(picks, 1) == 1)
        picks = picks';
    end
    
    %Applying PSF in Fourier domain
    if (length(dim) == 2)
        %Standard 2D image
    	xs = reshape(x, dim(1), dim(2));
        %Apply PSF in direct transform
	    xs = fft2(xs);
	    xs = xs.*OTF;
	    xs = ifft2(xs);    
	else
	    xs = reshape(x, dim(1), dim(2), dim(3));
        xs = fft(xs, [], 2);
        xs = fft(xs, [], 3);
        OTF = permute(OTF, [3 1 2]); %Adding singleton
        xs = bsxfun(@times, xs, OTF);

        xs = ifft(xs, [], 2);
        xs = ifft(xs, [], 3);
    end
   	%Reducing signal to relevant values
   	xs = real(xs(picks));	
    
    %Perform FFT (1/sqrt(K) for symmetric transform)
    fx = 1/sqrt(K) * fft(xs);
    %Encode values by taking non-redundant ones from real signal
    %0: constant value
    %1:K/2: symmetric part
    %K/2+1: pivot value
    %1 and K/2+1 are real (imag = 0)
    b = sqrt(2) * [real(fx(1:K/2+1, :)); imag(fx(2:K/2, :))];
    
end
