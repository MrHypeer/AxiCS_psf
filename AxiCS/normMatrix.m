function norm_value = normMatrix(input)
%%% this function morms input into range 0-255(uint8), canbe 3d and 2d input
%%% default input is (x,y,z)
norm_value = uint8(zeros(size(input)));
if length(size(input)) == 3
    % process 3d data
    for ii_stack = 1:size(input,3)
        max_int = max(max(input(:,:,ii_stack)));
        min_int = min(min(input(:,:,ii_stack)));
        norm_value(:,:,ii_stack) = (input(:,:,ii_stack)-min_int)./(max_int ...
            -min_int)*255;
    end
elseif length(size(size(input) == 2))
    % process 2d data
    max_int = max(max(input));
    min_int = min(min(input));
    norm_value = (input-min_int)./(max_int ...
        -min_int)*255;
else 
    % error
    error('Input dim shoule be 2D or 3D')
end
end
