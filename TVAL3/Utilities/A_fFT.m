% A_fFT.m
%
% Function to extract a subset of Fourier coefficients from shuffled data.
%
% Input:
% x             input data (1D)
% freq          frequencies to retain
% idxShuffle    indices (shuffled) from data (as obtained from randperm, for example)
%
% Output:
% b             measurement vector
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

function b = A_fFT(x, freq, idxShuffle)

N = length(x);

xs = x(idxShuffle);
ft = 1/sqrt(N)*fft(xs);
b = sqrt(2) * [real(ft(freq)); imag(ft(freq))];

end