% evaluate_type1errors_pcboot
%
% This script evaluates the type 1 error rate of the percentile bootstrap
% method, to ensure that it has a type 1 error rate at or close to the
% nominal alpha level. This can also be used to compare Student's
% paired-samples t test and Yuen's paired-samples t, to see whehter type 1
% error rates differ.
% 
% 

% Housekeeping
clear all;
close all;

% Settings
sampleSizesToUse = [10 15 20];
alpha = 0.05;
nIterations = 10000; % Number of tests to run
trimming = 20; % Percent trimming to use for Yuen's t

% Preallocate
Type1Errors.pcBoot = zeros(nIterations, length(sampleSizesToUse));
Type1Errors.students_t = zeros(nIterations, length(sampleSizesToUse));
Type1Errors.yuen = zeros(nIterations, length(sampleSizesToUse));

p_values.pcBoot = zeros(nIterations, length(sampleSizesToUse));
p_values.students_t = zeros(nIterations, length(sampleSizesToUse));
p_values.yuen = zeros(nIterations, length(sampleSizesToUse));

% Seed random number generator based on computer clock
rng('shuffle');

for sampleSize = 1:length(sampleSizesToUse)
    for iteration = 1:nIterations

        fprintf('Running sample size %i iteration %i...\n', sampleSizesToUse(sampleSize), iteration);

        % Generate a sample drawn from a normal distribution, mean 0
        temp = randn(sampleSizesToUse(sampleSize), 1);

        % Create dummy dataset of zeros (as only using one randomly
        % generated dataset)
        temp2 = zeros(sampleSizesToUse(sampleSize), 1);
        
        % Run percentile bootstrap
        Results = percentile_bootstrap(temp, temp2, 'trimming_percent', trimming);

        % Log type 1 errors
        Type1Errors.pcBoot(iteration, sampleSize) = Results.h;
        p_values.pcBoot(iteration, sampleSize) = Results.p;

        % Run student's t
        [Type1Errors.students_t(iteration, sampleSize), p_values.students_t(iteration, sampleSize)] = ttest(temp, 0, 'Alpha', alpha);

        % Run Yuen's t
        [Type1Errors.yuen(iteration, sampleSize), p_values.yuen(iteration, sampleSize)] = yuend_ttest(temp, temp2, trimming, alpha);
        
    end % of for iteration
end % of for sampleSize

for i = 1:length(sampleSizesToUse)
    Type1ErrorRate.pcBoot(i) = sum(Type1Errors.pcBoot(:, i)) / nIterations;
    Type1ErrorRate.students_t(i) = sum(Type1Errors.students_t(:, i)) / nIterations;
    Type1ErrorRate.yuen(i) = sum(Type1Errors.yuen(:, i)) / nIterations;
end


% Plot results
figure('Name', 'Type 1 Error Rates by Method');
plot(Type1ErrorRate.students_t);
hold on;
plot(Type1ErrorRate.yuen);
plot(Type1ErrorRate.pcBoot);


