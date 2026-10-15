function r = calculate_correlation(mask_X, mask_Y, X, Y, method)
    if nargin < 5
        method = 'spearman'; % Default
    end
    % 1. Safety Checks
    if sum(mask_X) == 0 || sum(mask_Y) == 0
        r = 0; return; 
    end
    
    selected_data_X = X(:, mask_X == 1);
    selected_data_Y = Y(:, mask_Y == 1);
    
    % 2. Calculate Representative Value per Participant
    if strcmp(method, 'pls')
        % --- STRATEGY 7: PARTIAL LEAST SQUARES (PLS) ---
        try
            [~, ~, Xscores, Yscores] = plsregress(selected_data_X, selected_data_Y, 1);
            mean_val_X = Xscores; 
            mean_val_Y = Yscores; 
        catch
            mean_val_X = mean(selected_data_X, 2); 
            mean_val_Y = mean(selected_data_Y, 2);
        end
        
    elseif strcmp(method, 'pca')
        % --- STRATEGY 6: PCA AGGREGATION ---
        try
            [coeff_X, score_X] = pca(selected_data_X, 'NumComponents', 1);
            mean_val_X = score_X; 
            
            [coeff_Y, score_Y] = pca(selected_data_Y, 'NumComponents', 1);
            mean_val_Y = score_Y;
            
            if isempty(mean_val_X) || isempty(mean_val_Y)
                 mean_val_X = mean(selected_data_X, 2); 
                 mean_val_Y = mean(selected_data_Y, 2);
            end
        catch
            mean_val_X = mean(selected_data_X, 2); 
            mean_val_Y = mean(selected_data_Y, 2);
        end
    else
        % Default: Simple Mean
        mean_val_X = mean(selected_data_X, 2);
        mean_val_Y = mean(selected_data_Y, 2);
    end
    
    clean_X = mean_val_X;
    clean_Y = mean_val_Y;
    
    % Safety Check: If we deleted too many people, return 0
    if length(clean_X) < 15
        r = 0; return;
    end

    % --- CRITICAL FIX 1: Check for "Blobbing" ---
    persistent ridge_std_check_shown;
    if strcmp(method, 'ridge') && isempty(ridge_std_check_shown)
        fprintf('[RIDGE STD CHECK] std(clean_X)=%.6f, std(clean_Y)=%.6f\n', std(clean_X), std(clean_Y));
        ridge_std_check_shown = true;
    end

    if std(clean_X) < 1e-5 || std(clean_Y) < 1e-5
        r = 0; 
    else
        if strcmp(method, 'robust')
            % --- STRATEGY 4: ROBUST REGRESSION ---
            try
                warning('off', 'stats:statrobustfit:IterationLimit'); 
                [~, stats] = robustfit(clean_X, clean_Y);
                warning('on', 'stats:statrobustfit:IterationLimit'); 
                
                rmse = stats.s; 
                r_val = 1 / (1 + rmse); 
            catch
                r_val = 0; 
            end
            
        elseif strcmp(method, 'ridge')
            % --- STRATEGY 8: RIDGE REGRESSION ---
            try
                weights = ridge(clean_Y, selected_data_X, 0.1, 0);
                % ridge(y, X, k, 0) returns [intercept; coeff_1; ...; coeff_p]
                if isrow(weights), weights = weights'; end
                intercept = weights(1);
                feature_weights = weights(2:end);
                y_pred = selected_data_X * feature_weights + intercept;
                rmse = sqrt(mean((clean_Y - y_pred).^2));
                r_val = 1 / (1 + rmse);
            catch ME
                r_val = 0;
            end
            
        elseif strcmp(method, 'lasso')
            % --- STRATEGY 9: LASSO REGRESSION ---
            try
                [weights, fitInfo] = lasso(selected_data_X, clean_Y, 'Lambda', 0.05);
                y_pred = selected_data_X * weights + fitInfo.Intercept;
                rmse = sqrt(mean((clean_Y - y_pred).^2));
                r_val = 1 / (1 + rmse);
            catch
                r_val = 0;
            end
            
        elseif strcmp(method, 'svr')
            % --- STRATEGY 10: SUPPORT VECTOR REGRESSION ---
            try
                svm_mdl = fitrsvm(selected_data_X, clean_Y, 'KernelFunction', 'linear', 'Standardize', true);
                y_pred = predict(svm_mdl, selected_data_X);
                rmse = sqrt(mean((clean_Y - y_pred).^2));
                r_val = 1 / (1 + rmse);
            catch
                r_val = 0;
            end
            
        else
            % Default: Spearman Correlation
            r_val = corr(clean_X, clean_Y, 'Type', 'Spearman', 'Rows', 'complete');
        end
        
        r = abs(r_val);
        
        % --- STRATEGY 1: SPARSITY PENALTY ---
        USE_SPARSITY_PENALTY = 0;  % Set to 1 to enable, 0 to disable

        n_genes_x = sum(mask_X);
        n_genes_y = sum(mask_Y);
        total_genes = length(mask_X) + length(mask_Y);

        penalty_factor = 0.1;  % Sparsity penalty: uniform across all methods
        sparsity_cost = USE_SPARSITY_PENALTY * penalty_factor * ((n_genes_x + n_genes_y) / total_genes);

        r = r - sparsity_cost;
        
        if r < 0, r = 0; end
    end
end
