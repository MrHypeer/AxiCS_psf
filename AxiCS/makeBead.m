function image = makeBead(sparsity, show)
%%%%%%%%
% this code generates bead sample in 512*512 image
% density can be controlled using
%%%%%%%%

%% define params
x_pxl = 512; y_pxl = 512;

%% generate image
image = zeros(x_pxl, y_pxl);
point2gen = floor(sparsity/10000*x_pxl*y_pxl); % sparsity: x over 100*100 pxls

%Comment out for different random realizations
if (~exist('randSeed_axics.mat', 'file'))
    s = rng;
    save('randSeed_axics.mat',  's');
else
    load randSeed_axics.mat
end
rng(s);

idx = randperm(x_pxl*y_pxl, point2gen);
image(idx) = 1;

if show
    figure
    imagesc(image); title('random point generation')
    xlabel('x/pxl'); ylabel('y/pxl')
end
end

