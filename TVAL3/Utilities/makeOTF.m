%makeOTF.m
%
%Generates an OTF function from its corresponding PSF. The function is
%essentially a convenience function, taking care of the circular shifts and
%resizing involved.
%The function also takes care of Fourier truncation, and mask apodization
%to reduce spatial ripples.

%Input
% PSF       Spatial matrix containing the PSF function, centered. Can be
%           small and contain only significant values
% dim       Dimensions of image on which to use OTF
% lim       Limit for truncation, defined as the smallest allowed value in
%           OTF intensity. Setting to 0 disables truncation
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

function [OTF, mask] = makeOTF(PSF, dim, lim)

if (sum(mod(dim, 2)) > 0)
    error('Dimensions must be even');
end

%Normalizing PSF if required
if (sum(PSF(:)) ~= 1)
    PSF = PSF./sum(PSF(:));
end

%Preparing sizes
OTF = zeros(dim);
N = size(PSF, 1);
N2 = floor(N/2)+1;
M = size(PSF, 2);
M2 = floor(M/2)+1;

%Putting PSF at origin
OTF(1:N2, 1:M2) = PSF(N2:end, M2:end);
OTF(end-N2+2:end, 1:M2) = PSF(1:N2-1, M2:end);
OTF(1:N2, end-M2+2:end) = PSF(N2:end, 1:M2-1);
OTF(end-N2+2:end, end-M2+2:end) = PSF(1:N2-1, 1:M2-1);

%Tranferring to Fourier domain
OTF = fft2(OTF);

%Finding limit positions
limX = find(abs(OTF(:, 1)) < lim, 1)-1;
limY = find(abs(OTF(1, :)) < lim, 1)-1;
if (isempty(limX))
    limX = size(OTF, 1);
end
if (isempty(limY))
    limY = size(OTF, 2);
end

%Elliptic binary mask
N2 = floor(size(OTF, 1)/2);
M2 = floor(size(OTF, 2)/2);
[x, y] = meshgrid(-N2:N2-1, -M2:M2-1);
dist = sqrt( (x/limX).^2 + (y/limY).^2)';
if (~isempty(limX) && ~isempty(limY))
    mask = double(dist<1);
else
    mask = ones(size(dist));
end

%Smoothing border to minimize ringing
sk = 11;
ker = ones(sk, sk)./sk.^2;
mask = conv2(mask, ker, 'same');
%Putting mask at origin
mask = fftshift(mask);

end
