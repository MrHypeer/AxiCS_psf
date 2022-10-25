% A_fPsfD.m
% 
%Adjoint of function to encode sparse spatial data into non-sparse Fourier
%coefficients, with PSF taken into account in Fourier domain
%through deconvolution.
%
% Input:
% b             measurement vector
% dim           Dimensions of original data ([x, y] vector for standard image)
% picks         indices of known spatial values within full size data
% OTF           optical transfer function, as obtained from makeOTF() function
% 
% Output:
% x             sparse function (1D)
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

function x = At_fPsfD(b, dim, picks, OTF, mask)
    K = length(picks);

    if (length(dim) == 2)
    	N = dim(1)*dim(2);
    else
        N = dim(1)*dim(2)*dim(3);
    end
    
	if ( floor(K/2) ~= ceil(K/2) )
        error('The amount of samples must be even');
	end
    
    %Reconstructing Fourier signal from non-redundant observations
    fx = zeros(K, 1);
    %Real part, with pivot values
    fx(1:K/2+1) = b(1:K/2+1)/sqrt(2); %1/sqrt(2) for proper scaling
	%Imaginary part, taking into account 0 values
	fx(2:K/2) = fx(2:K/2) + 1i*b(K/2+2:K)/sqrt(2); %1/sqrt(2) for proper scaling
	%Copying symmetric real part
	fx(K/2+2:end) = real(fx(K/2:-1:2));
	%Copying anti-symmetric imaginary part
	fx(K/2+2:end) = fx(K/2+2:end) - 1i*imag(fx(K/2:-1:2));
	
	%Passing in spatial again
	x = zeros(N,1);
	%sqrt(K) factor for symmetrical transform
	x(picks) = sqrt(K)*real(ifft(fx));
	
	%Applying PSF in 2D Fourier domain
    if (length(dim) == 2)
        xs = reshape(x, dim(1), dim(2));
        
        %Regularized division
        xs = fft2(xs)./OTF;
        %Truncating the spectrum to remove amplified noise
        xs = xs.*mask;   
        xs = ifft2(xs);
    else
        xs = reshape(x, dim(1), dim(2), dim(3));
        xs = fft(xs, [], 2);
        xs = fft(xs, [], 3);
        OTF = permute(OTF, [3 1 2]); %Adding singleton
        %Regularized division
        xs = bsxfun(@rdivide, xs, OTF);
        %Truncating spatial spectrum to remove noise
        mask = permute(mask, [3 1 2]); %Adding singleton
        xs = bsxfun(@times, xs, mask);
        
        xs = ifft(xs, [], 2);
        xs = ifft(xs, [], 3);   
    end 
    
    %Putting in line again, and removing imag part 
    %(should be negligible as our data is real)
    x = real(xs(:));
end
