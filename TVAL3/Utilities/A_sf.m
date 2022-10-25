% A_sf.m
% 
% Function to encode sparse spatial data into non-sparse Fourier coefficients.
%
% Input:
% x             sparse function (1D)
% picks         indices of known spatial values within full size data
%
% Output:
% b             measurement vector
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

function b = A_sf(x, picks)
    %Size of measurement vector
    K = length(picks);

    if (K > length(x))
        error('The amount of samples must be smaller than full data size');
    end
    
    %Reducing signal to relevant values
    xs = x(picks);
    %Perform FFT (1/sqrt(K) for symmetric transform)
    fx = 1/sqrt(K) * fft(xs);
    
    %Encode values by taking non-redundant ones from real signal
    %0: constant value
    %1:K/2: symmetric part
    %K/2+1: pivot value
    %1 and K/2+1 are real (imag = 0)
    b = sqrt(2) * [real(fx(1:K/2+1, :)); imag(fx(2:K/2, :))];  
end 