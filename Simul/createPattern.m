% createPattern.m
% 
% Function to prepare a simulated image for estimating algorithm performance
%
% Input:
% show          variable determining if intermediate images should be shown
% 
% Output:
% pat           created pattern image
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

function pat = createPattern(show)

if (~exist('show', 'var'))
    show = 0;
end

%% Pattern 1
x = linspace(0, 12*pi, 512);
xx = repmat(x, 512, 1);

xx = imrotate(xx, 10);

xx = angle(exp(1i*xx));
pat1 = (xx(129:384, 129:384)+pi)./(2*pi);
if (show)
    subplot(221), imagesc(pat1);
end

%% Pattern 2
loc = [225, 205, 40;
       65, 95, 30;
       165, 120, 20;
       65, 225, 25];
       
pat2 = zeros(306, 306);
[x, y] = meshgrid(1:306, 1:306);
for k = 1:size(loc, 1)
    dist = sqrt( (x-loc(k, 1)).^2 + (y-loc(k, 2)).^2 ) < loc(k, 3);
    pat2 = pat2 + dist;
end
x = linspace(-1, 1, 25);
[x, y] = meshgrid(x, x);
kernel = exp( -( x.^2 + y.^2 )./(2*0.5^2) );

pat2 = conv2(pat2, kernel, 'same');
pat2 = pat2(26:281, 26:281);
pat2 = pat2./max(pat2(:));
if (show)
    subplot(222), imagesc(pat2);
end

%% Pattern 3
loc = [160, 320;
       130, 260;
       210, 340;
       320, 200;
       250, 220;
       210, 170];

pat3 = zeros(512, 512);
for k = 1:size(loc, 1)
    pat3( loc(k, 1)-5:loc(k, 1)+4 , loc(k, 2)-15:loc(k, 2)+14) = 1;
end

pat3 = imrotate(pat3, 15);

pat3 = pat3(158:413, 158:413);
if (show)
    subplot(223), imagesc(pat3);
end

pat = pat1 + pat2 + pat3;
%Normalization (pat2 and pat3 do not superimpose)
pat = pat./2;

if (show)
    subplot(224), imagesc(pat)
end

end