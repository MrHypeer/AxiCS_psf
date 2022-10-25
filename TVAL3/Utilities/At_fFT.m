% At_fFT.m
%
% Adjoint of function to extract a subset of Fourier coefficients from shuffled data.
%
% Input:
% b             measurement vector
% N             size of output data
% freq          frequencies to retain
% idxShuffle    indices (shuffled) from data (as obtained from randperm, for example)
%
% Output:
% x             output data (1D)
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

function x = At_fFT(b, N, freq, idxShuffle)

M = length(b);
ft = zeros(N,1);
%Put known frequencies in vector
ft(freq) = sqrt(2)*b(1:M/2) + 1i*sqrt(2)*b(M/2+1:M);
x = zeros(N,1);
%Remove imag part (should be negligible as x is real)
xs = sqrt(N)*real(ifft(ft));

%Shuffle back
x(idxShuffle) = xs;

end

