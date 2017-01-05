function [Results] = percentile_bootstrap(cond1_data, cond2_data, varargin)
% 
% Function to run a percentile bootstrap on single-sample or paired-samples
% data. For each of B iterations a bootstrap sample is drawn. For
% paired-samples cases  in this function a bootstrap sample is taken from 
% the difference scores. The X% trimmed mean is calculated for each iteration
% (default 20% trimming) and (1 - alpha) percent confidence intervals are 
% calculated for the resulting bootstrap distribution (default 95% confidence intervals).
%
%
% PUT WILCOX REF HERE
%
% Inputs:
%
%   cond1_data          A vector of data from condition 1, with the length
%                       equal to the number of subjects in the dataset.
%
%   cond2_data          A vector of data from condition 2, with the length
%                       equal to the number of subjects in the dataset
%                       (This is not needed for single sample tests).
%
%   null_value          The null hypothesis of the difference between
%                       conditions (or size of the effect for single-samples tests). 
%                       In most cases this is set to zero (null hypothesis
%                       of zero difference or effect).
%
%   alpha               Nominal alpha level. This determines the width of 
%                       the (two-tailed) confidence interval. Default is
%                       0.05, corresponding to confidence interval limits of 
%                       the 2.5th and 97.5th percentiles of the bootstrap
%                       distribution.
% 
%   n_bootstrap_samples     Number of bootstrap samples to draw. Default is 10000.
%                           At least 1000 is recommended for the p = 0.05 alpha level,
%                           and at least 5000 is recommended for the p = 0.01 alpha level.
%                           This is due to extreme events at the tails being very rare, 
%                           needing many bootstrap samples to accurately
%                           estimate them.
% 
%   trimming_percent    Percentage of highest and lowest value to trim when
%                       using the trimmed mean. Must be between 0 and 50%.
%                       Default trimming is 20%.
%                       
%
% Outputs:
% "Results" structure containing:
%   
%   trimmed_mean_difference     The X% trimmed mean of the difference
%                               scores (or trimmed mean of the single
%                               sample data).
%
%   CI          (1 - alpha / 2) percent confidence interval derived from
%               the bootstrap distribution
%
%   h           Outcome of null-hypothesis statistical testing based on
%               the derived confidence interval. 1 = reject null
%               hypothesis, 0 = do not reject null hypothesis.
%
%   p           The p-value derived from the bootstrap distribution. This
%               is calculated as the percentage of bootstrap sample trimmed 
%               means that are lower/higher than the null hypothesis value if the 
%               20% trimmed mean difference of the original data is larger/smaller 
%               the null hypothesis value.  
%               
%   This structure a record of settings used: n_bootstrap_samples, alpha, trimming_percent
%
%
% Example:      [Results] = percentile_bootstrap(cond1_data, cond2_data, 'null_value', 0, 'alpha', 0, 'n_bootstrap_samples', 10000, 'trimming_percent', 20)
%
% 
% Copyright (c) 2016 Daniel Feuerriegel
%
%
% This file is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%


%% Defaults (if not specified in input arguments)
if nargin < 2 % if condition 2 data not given (single sample case)
    cond2_data = zeros(size(cond1_data, 1), size(cond1_data, 2));
end

%% Handling variadic inputs
% Define defaults at the beginning
options = struct(...
    'null_value', 0,...
    'alpha', 0.05,...
    'n_bootstrap_samples', 10000,...
    'trimming_percent', 20);
    
% Read the acceptable names
option_names = fieldnames(options);

% Count arguments
n_args = length(varargin);
if round(n_args/2) ~= n_args/2
   error([mfilename ' needs property name/property value pairs'])
end

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
   inp_name = lower(pair{1}); % make case insensitive

   % Overwrite default options
   if any(strcmp(inp_name, option_names))
      options.(inp_name) = pair{2};
   else
      error('%s is not a recognized parameter name', inp_name)
   end
end
clear pair
clear inp_name

% Renaming variables for use below:
null_value = options.null_value;
alpha = options.alpha;
n_bootstrap_samples = options.n_bootstrap_samples;
trimming_percent = options.trimming_percent;
clear options;

if trimming_percent < 0 || trimming_percent > 50 % if trimming percentage over 50 or less than zero is chosen
    error('trimming parameter needs to be between 0 and 50 percent')
end

%% Seed the random number generator
rng('shuffle'); % Seed random number generator based on computer clock

%% Bootstrap procedure
% Set percentiles based on alpha level (simple multiplication)
significance_threshold = alpha * 100; % 05 for p < 0.05; 01 for p < 0.01; 0.1 for p < 0.001

% Calculate the difference between condition 1 and 2 datasets (or condition
% 1 and zero for single-sample case)
diff_scores = cond1_data - cond2_data;

% Calculate the trimmed mean of the difference scores
trimmed_mean_difference = trimmean(diff_scores, trimming_percent);

% Draw the required number of bootstrap samples
[~, bootsam] = bootstrp(n_bootstrap_samples, 'mean', diff_scores(:));

% Preallocate matrices
bootstrap_samples_data = zeros(size(bootsam, 1), size(bootsam, 2));
bootstrap_trimmed_mean = zeros(size(bootsam, 1), 1);

for k = 1:size(bootsam, 2)
    
    for l = 1:size(bootsam, 1);
        bootstrap_samples_data(l,k) = diff_scores(bootsam(l, k));
    end % for for l
    
    % Calculate the trimmed mean for each bootstrap sample
    bootstrap_trimmed_mean(k) = trimmean(bootstrap_samples_data(:, k), trimming_percent);
    
end % of for k

% Sort trimmed means from bootstrap samples
sorted_trimmed_means = sort(bootstrap_trimmed_mean);

% Calculate confidence interval limits based on percentiles
CI(1) = prctile(sorted_trimmed_means, significance_threshold / 2);
CI(2) = prctile(sorted_trimmed_means, (100 - significance_threshold / 2));

% Calculate the p-value based on the number of bootstrap samples that do
% not include zero
if trimmed_mean_difference > null_value % if a positive difference direction
    p = mean(sorted_trimmed_means(:) <= null_value);
    p = p .* 2; % Two-tailed correction
    if p > 1 % If p is larger than 1 then set to 1
        p = 1;
    end
    
elseif trimmed_mean_difference < null_value % if a negative difference direction
    p = mean(sorted_trimmed_means(:) >= null_value);
    p = p .* 2; % Two-tailed correction
    if p > 1 % If p is larger than 1 then set to 1
        p = 1;
    end
    
elseif trimmed_mean_difference == 0 % if difference is exactly zero
    p = 1;
end

% Determine statistical significance based on p-values
if p < alpha
    h = 1;
else
    h = 0;
end


% Copy settings and results into a structure to output
Results.trimmed_mean_difference = trimmed_mean_difference;
Results.CI = CI;
Results.h = h;
Results.p = p;
Results.n_bootstrap_samples = n_bootstrap_samples;
Results.alpha = alpha;
Results.trimming_percent = trimming_percent;

end % of function