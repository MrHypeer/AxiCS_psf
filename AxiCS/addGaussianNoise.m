function img = addGaussianNoise(data, mean, var, txtname)
    % add gaussian noise on image with designated params
    % input: data in x-y format
    % mean/var, gaussian noise params
    % txtname, to which noise params will be stored
    data_x = size(data,1);
    data_y = size(data,2);
    noise_set = length(mean)*length(var);
    img = zeros(data_x, data_y, noise_set); % to store output with noise
    
    data = data-min(data(:));
    data = data/max(data(:))/1.2; % divide 1.2 for noise addition
    
    % open file to store noise information
    fileID = fopen(txtname,'a');
    
    count = 1;
    for ii_mean = 1:length(mean)
       for ii_var = 1:length(var)
           count = count + 1;
           img(:,:,(ii_mean-1)*length(var)+ii_var) = imnoise(data, 'gaussian', ...
               mean(ii_mean), var(ii_var));
           formatSpec = '%2.0u: gaussian nosie mean/var is %4.2f/%4.2f\n';
           fprintf(fileID, formatSpec,...
                [count, mean(ii_mean), var(ii_var)]);
       end
    end
    fclose(fileID); % close file
end
