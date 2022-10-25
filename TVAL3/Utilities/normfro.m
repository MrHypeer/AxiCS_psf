% normfro.m
% 
% Function to compute Frobenius norm for ND matrices (Matlab built-in works
% only for 2D).
%
% %The Matlab implementation is based on the SVD definition, which is not
% usable for 3D matrix, we use use the simple sum definition.
%
% Input:
% m             input matrix
%
% Output:
% normVal       norm
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

function normVal = normfro(m)
    normVal = sqrt(sum(m(:).^2));
end


 