% At_sf.m
% 
% Adjoint of function to encode sparse spatial data into non-sparse 
% Fourier coefficients.
%
% Input:
% b             measurement vector
% N             original signal size
% picks         indices of known spatial values within full size data
% 
% Output:
% x             sparse function (1D)
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

function x = At_sf(b, N, picks)
    %Length of measurement vector
    K = length(picks);

    if (K > N)
        error('The amount of samples must be smaller than full data size');
    end
    if (K ~= length(b))
        error('The length of indices and measurement vector must be identical');
    end
    
    if ( floor(K/2) ~= ceil(K/2) )
        error('The amount of samples must be even');
    end
    
    %Recreating Fourier signal from redundant values
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
end
