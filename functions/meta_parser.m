function DeploymentInfo = meta_parser(meta)
% Parse meta data from BEAMS for use in post-processing
    DeploymentInfo = struct;

    idx = find(contains(meta, 'k0ext:'));
    
    if ~isempty(idx)
        % Extract the string containing 'k0ext:'
        str = meta{idx};
    
        % Use regular expression to parse the numeric value
        DeploymentInfo.k0ext = sscanf(str,' k0ext:        	%f');
    else
        disp('k0ext: not found in meta.');
    end
    
    % Parse k2
    idx = find(contains(meta, 'k2ext:'));
    
    if ~isempty(idx)
        % Extract the string containing 'k0ext:'
        str = meta{idx};
    
        % Use regular expression to parse the numeric value
        DeploymentInfo.k2ext = sscanf(str,' k2ext:        	%f');
    else
        disp('k0ext: not found in meta.');
    end

    % Parse number of reps per pump cycle
    idx = find(contains(meta, 'Pump1 cycles:'));
    
    if ~isempty(idx)
        % Extract the string containing 'k0ext:'
        str = meta{idx};
    
        % Use regular expression to parse the numeric value
        DeploymentInfo.nreps = sscanf(str,' Pump1 cycles:     	%f');
    else
        disp('Pump1 cycles not found in meta.');
    end
    
    % Parse sample period
    idx = find(contains(meta, ' Sampling period (s)'));
    
    if ~isempty(idx)
        % Extract the string containing 'k0ext:'
        str = meta{idx};
    
        % Use regular expression to parse the numeric value
        DeploymentInfo.sampleperiod = sscanf(str,' Sampling period (s)	%f');
    else
        disp('Sample period not found in meta.');
    end
  

    % Parse sample interval
    idx = find(contains(meta, ' Pump1 on time (s):'));
    
    if ~isempty(idx)
        % Extract the string containing 'k0ext:'
        str = meta{idx};
    
        % Use regular expression to parse the numeric value
        DeploymentInfo.sampinterval = sscanf(str,' Pump1 on time (s):	%f');
    else
        disp('Sample interval not found in meta.');
    end

    % min2avg = sampinterval * nreps / 60; 