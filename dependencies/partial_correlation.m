function [residual_data, residual_behav] = partial_correlation(data, behavioral, covariates)
% PARTIAL_CORRELATION Residualize data and behavioral variables for partial correlation
%
% This function removes the linear effects of covariates from both the
% neural data and behavioral variables, enabling partial correlation analysis
% that controls for confounding variables.
%
% Inputs:
%   data - Neural data matrix [n_subjects x n_channels]
%   behavioral - Behavioral variable vector [n_subjects x 1]
%   covariates - Covariate matrix [n_subjects x n_covariates]
%
% Outputs:
%   residual_data - Residualized neural data [n_subjects x n_channels]
%   residual_behav - Residualized behavioral variable [n_subjects x 1]
%
% Method:
%   For each variable (behavioral and each channel), we:
%   1. Fit a linear model: Y = X*beta + residuals
%      where X = [ones(n,1), covariates] (design matrix with intercept)
%   2. Extract residuals: residuals = Y - X*beta
%   3. Use residuals for correlation analysis
%
% This implements partial correlation by removing shared variance with
% covariates before computing correlations.

    n_subjects = size(data, 1);
    n_channels = size(data, 2);
    
    % Create design matrix with intercept
    X = [ones(n_subjects, 1), covariates];
    
    % Residualize behavioral variable
    beta_behav = X \ behavioral;
    residual_behav = behavioral - X * beta_behav;
    
    % Residualize each channel
    residual_data = zeros(size(data));
    for ch = 1:n_channels
        y = data(:, ch);
        
        % Skip if all NaN
        if all(isnan(y))
            residual_data(:, ch) = y;
            continue;
        end
        
        % Handle NaN values
        valid_idx = ~isnan(y);
        if sum(valid_idx) < size(X, 2) + 1
            % Not enough valid data points for regression
            residual_data(:, ch) = y;
            continue;
        end
        
        % Fit model on valid data
        beta_ch = X(valid_idx, :) \ y(valid_idx);
        
        % Compute residuals
        residuals = nan(n_subjects, 1);
        residuals(valid_idx) = y(valid_idx) - X(valid_idx, :) * beta_ch;
        
        residual_data(:, ch) = residuals;
    end
    
    fprintf('Partial correlation: Residualized %d channels and behavioral variable\n', n_channels);
end