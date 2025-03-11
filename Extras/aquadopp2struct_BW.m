function AD = aquadopp2struct_BW(filepath)
    % Read CSV file using readtable
    opts = detectImportOptions(filepath, 'Delimiter', ';');
    T = readtable(filepath, opts);
    
    % Extract column names
    headers = T.Properties.VariableNames;
    
    % Initialize struct
    AD = struct();
    
    % Convert DateTime to MATLAB datenum
    AD.SDN = datenum(T.DateTime);    
    % Assign standard variables
    AD.Batt     = T.Battery;
    AD.Heading  = T.Heading;
    AD.Pitch    = T.Pitch;
    AD.Roll     = T.Roll;
    AD.Pres     = T.Pressure;
    AD.TC       = T.Temperature;
    
    % Identify speed and direction columns and extract depth bins
    speed_cols = contains(headers, 'Speed');
    dir_cols = contains(headers, 'Dir');
    
    bin_depth = nan(sum(speed_cols), 1);
    speed_data = nan(height(T), sum(speed_cols));
    dir_data = nan(height(T), sum(dir_cols));
    
    speed_names = headers(speed_cols);
    dir_names = headers(dir_cols);
    
    for i = 1:length(speed_names)
        % Extract depth bin from column name
        tokens = regexp(speed_names{i}, '\((\d+\.\d+)m\)', 'tokens');
        if ~isempty(tokens)
            bin_depth(i) = str2double(tokens{1}{1});
        end
        
        % Store speed and direction data
        speed_data(:, i) = T.(speed_names{i});
        dir_data(:, i) = T.(dir_names{i});
    end
    
    % Assign processed data to struct
    AD.bin_depth = bin_depth;
    AD.vel = speed_data';
    AD.dir = dir_data';
    
    % Compute u and v components
    AD.uspd = nan(size(AD.vel));
    AD.vspd = nan(size(AD.vel));
    
    for i = 1:length(AD.bin_depth)
        AD.dir(i, :) = pol2compass(AD.dir(i, :));
        [u, v] = compass2cart(AD.dir(i, :), AD.vel(i, :));
        AD.uspd(i, :) = u;
        AD.vspd(i, :) = v;
    end
    
    % Store variable names
    AD.vars = {'Batt', 'SDN', 'Heading', 'Pitch', 'Roll', 'TC', 'Pres', 'vel', 'dir', 'uspd', 'vspd', 'bin_depth'};
    AD.vector_vars = {'Batt', 'SDN', 'Heading', 'Pitch', 'Roll', 'TC', 'Pres'};
    AD.matrix_vars = {'vel', 'dir', 'uspd', 'vspd'};
end

