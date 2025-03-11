function AD = aquadopp2mat(filepath, delimit, daterange)

% AD = aquadopp2mat(filepath)
% 
% depth_res is the depth resolution of the aquadopp in cm. Usually should
% be 3cm.
%
% Assumes a CSV file, i.e. delimiter must be commas.
%
% Yui Takeshita
% Scripps Insitution of Oceanography
% Created: 3/12/2014

% I had to change the first column of aquadopp file (date&time) to numeric
% to use dlmread. I still can't figure out how to efficiently use
% textscan...


% filename = 'test.csv';
% trex = dlmread(filename, ',', 1,0);
% SDN = dlmread(filename, ',', [1, 0, size(trex,1), 0]);

% % textscan... I want to be able to write it so that I can just get a mxn
% % matrix of string. 

% ************OPEN FILE, EXTRACT DEPTH BIN AND DATA************************
file_format = '%s'; % declare file_format;
fid = fopen(filepath); % create file id
headers = textscan(fgetl(fid), '%s', 'delimiter', delimit); % Calculate how many columns file has.
headers = headers{:}; % put each header into its own cell.
n_col = length(headers); % number of columns
% extract depths of bins
i_bins = 10:2:n_col;
bin_depth = nan(length(i_bins), 1); % presize vector 
for i = 1:length(i_bins)
    header_str = headers{i_bins(i)};
    istart = strfind(header_str, '(')+1;
    iend = strfind(header_str, ')')-2;
    bin_depth(i) = str2double(header_str(istart:iend));
end
% go back to first line
frewind(fid); 
for i = 1:n_col-1; file_format = [file_format,'%f']; end; % loop for file format
trex = textscan(fid, file_format, 'headerlines', 1, 'delimiter', delimit, 'collectoutput', 1);
% trex = trex{:};
fclose(fid);
% ************************************************************************

% ======================Build data structure===============================

% Separate into numerical and string. in this case, the first column is
% string (SDN), and the rest are numeric.
AD.SDN = datenum(trex{1});
% Start numeric import
trex = trex{2};
AD.Batt     = trex(:,1);
AD.Heading  = trex(:,2);
AD.Pitch    = trex(:,3);
AD.Roll     = trex(:,4); 
AD.Pres     = trex(:,5);
AD.TC       = trex(:,6);
AD.vel      = trex(:,9:2:end)';
AD.dir      = trex(:,10:2:end)';
% AD.veldpth  = [10+depth_res:depth_res:10+size(AD.vel,1)*depth_res]';
AD.bin_depth= bin_depth;

AD.uspd = NaN(size(AD.vel));
AD.vspd = NaN(size(AD.vel));

% Now calculate u and v components. 
for i = 1:length(AD.bin_depth)
    AD.dir(i,:)  = pol2compass(AD.dir(i,:));
    [u, v] = compass2cart(AD.dir(i,:), AD.vel(i,:));
    AD.uspd(i,:) = u;
    AD.vspd(i,:) = v;   
end

AD.uv = sqrt(AD.uspd.^2 + AD.vspd.^2); %BW edit
% 
% 
% % If angle is in Q1 (0-90) or Q3 (180-270), then u is sin(theta). otherwise
% % cos(theta).
% % Get indicies for each quadrant.
% iQ1 = AD.dir <= 90;
% iQ2 = AD.dir > 90  & AD.dir <= 180;
% iQ3 = AD.dir > 180 & AD.dir <= 270;
% iQ4 = AD.dir > 270;
% % indicies to use for trigonometry calculations
% iQ1Q3 = iQ1|iQ3;
% iQ2Q4 = iQ2|iQ4;
% % presize matricies
% AD.uspd = NaN(size(AD.vel));
% AD.vspd = NaN(size(AD.vel));
% % Calculate speed in u and v direction.
% AD.uspd(iQ1Q3) = sin(AD.dir(iQ1Q3)).*AD.vel(iQ1Q3);
% AD.uspd(iQ2Q4) = cos(AD.dir(iQ2Q4)).*AD.vel(iQ2Q4);
% AD.vspd(iQ1Q3) = cos(AD.dir(iQ1Q3)).*AD.vel(iQ1Q3);
% AD.vspd(iQ2Q4) = sin(AD.dir(iQ2Q4)).*AD.vel(iQ2Q4);

iuse = inrange(AD.SDN, daterange);

AD.SDN(~iuse) = [];
AD.Batt(~iuse) = [];
AD.Heading(~iuse) = [];
AD.Pitch(~iuse) = [];
AD.Roll(~iuse) = [];
AD.Pres(~iuse) = [];
AD.TC(~iuse) = [];
AD.vel(:,~iuse) = [];
AD.dir(:,~iuse) = [];
AD.uspd(:,~iuse) = [];
AD.vspd(:,~iuse) = [];

vector_vars = {'Batt', 'SDN', 'Heading', 'Pitch', 'Roll', 'TC', 'Pres'};
matrix_vars = {'vel', 'dir', 'uspd', 'vspd', 'uv'}; % BW edit added uv
AD.vars = {vector_vars{:}, matrix_vars{:}, 'bin_depth'};
AD.vector_vars = vector_vars;
AD.matrix_vars = matrix_vars;

return













