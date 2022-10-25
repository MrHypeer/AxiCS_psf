% setVersion.m
% 
% Simple function to set and display a version number for the package.
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

function version = setVersion(show)

version = 'Laser-scanning CS, based on TVAL3. Version 1.1';

if (exist('show', 'var') && show)
    disp(version);
end

end