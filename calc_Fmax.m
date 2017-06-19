%Calculate Fmax permutation effect. This function is called by FmaxGND.
%For more inforamtion see:
%>> help FmaxGND
%
%REQUIRED INPUTS
% data          - An electrode x time points x conditions x subjects array of ERP
%                 data. Array will vary in number of dimensions based on how many
%                 factors there are
% dims          - Dimensions of the data array involved in the effect to be
%                 calculated. For example, if data is an electrode x time points
%                 x Factor A x Factor B x subjects array and you want to
%                 calculated the main effect of A, dims = 3. If you want to
%                 calculate the AxB interaciton, dims  = [3, 4].
% n_perm        - Number of permutations to used to calculate null distribution
% alpha         - Alpha level of the test
% int_method    - A string that should be either 'exact' or 'approximate'.
%                 If 'exact', the method of restricted permutations will
%                 be used to conduct a test that controls the Type I error
%                 rate at alpha (assuming enough permutations). 
%                 If 'approximate', the method of permutation of residuals 
%                 will be used to conduct a test with Type I error rate 
%                 asymptotic to alpha as noise decreases and/or number of 
%                 subjects increases.
%
%OUTPUT
% test_results - A struct with results of the Fmax test
%
%
%VERSION DATE: 12 June 2016
%AUTHOR: Eric Fields, Tufts University (Eric.Fields@tufts.edu)
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 
%This is a beta version of this software. It needs additional testing 
%and SHOULD NOT be considered error free.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.
%This function incorporates some code from the Mass Univariate Toolbox, 
%Copyright (c) 2015, David Groppe: https://github.com/dmgroppe/Mass_Univariate_ERP_Toolbox/blob/master/LICENSE

%%%%%%%%%%%%%%%%%%%  REVISION LOG   %%%%%%%%%%%%%%%%%%%
% 11/28/16  - Moved from FmaxGND. See FmaxGND revision log for information
%             on earlier versions
% 6/12/17   - Added estimated alpha; added verblevel reports



function test_results = calc_Fmax(data, dims, n_perm, alpha, int_method)
%Reduce data as a necessary and send to the right function for statistical
%test
% 1. Average across any factors not involved in the effect being calculated
% 2. For interaction effects using the exact method, subtract across
%    2-level factors until the data is reduced to a one-way design
% 3. Send the reduced data to the appropriate ANOVA function
    
    %Average across factors not involved in this effect
    if length(dims) < ndims(data) - 3
        %Put the factors to average across as the initial dimensions
        num_dims = ndims(data);
        all_dims = 1:ndims(data);
        reorder = [all_dims(~ismember(all_dims, [1,2,dims,num_dims])), all_dims(ismember(all_dims, [1,2,dims,num_dims]))];
        reduced_data = permute(data, reorder);
        %Reduce all the factors to average across to a single dimension
        dim_sizes = size(reduced_data);
        num_dims_to_avg = ndims(data) - length(dims) - 3;
        reduced_data = reshape(reduced_data, [prod(dim_sizes(1:num_dims_to_avg)), dim_sizes((num_dims_to_avg+1):end)]);
        %Take the mean across that dimension
        reduced_data = mean(reduced_data, 1);
        %Get rid of the extra singleton dimension
        reduced_data = reshape(reduced_data, dim_sizes((num_dims_to_avg+1):end));
    else
        reduced_data = data;
    end
    
    %For exact interactions, reduce data to one-way design via subtraction
    if length(dims) > 1 && strcmpi(int_method, 'exact')
        %Get/check design structure
        dim_sizes = size(reduced_data);
        factor_levels = dim_sizes(3:(ndims(reduced_data)-1));
        assert(sum(factor_levels > 2) < 2);
        %Put the factors to subtract across as the initial dimensions
        [~, factor_order] = sort(factor_levels);
        reorder = [factor_order(1:end-1)+2, 1, 2, factor_order(end)+2 ndims(reduced_data)];
        reduced_data = permute(reduced_data, reorder);
        %Flatten data and subtract until the data is reduced to the right size
        reduced_data = reshape(reduced_data, 1, []);
        while length(reduced_data) > size(data,1)*size(data,2)*max(factor_levels)*size(data,ndims(data))
            reduced_data = reduced_data(1:2:length(reduced_data)) - reduced_data(2:2:length(reduced_data));
        end
        %Put back in regular format
        reduced_data = reshape(reduced_data, [size(data,1), size(data,2), max(factor_levels), size(data,ndims(data))]);
    end

    %Calculate appropriate ANOVA
    if length(dims) == 1 || strcmpi(int_method, 'exact')
        test_results = oneway(reduced_data, n_perm, alpha);
    elseif length(dims) == 2
        test_results = twoway_approx_int(reduced_data, n_perm, alpha);
    elseif length(dims) == 3
        test_results = threeway_approx_int(reduced_data, n_perm, alpha);
    end

end


function test_results = oneway(data, n_perm, alpha)
%Permutation one-way ANOVA
% 1. Randomly permute conditions within each subject across all time
%    points and electrodes
% 2. For each permutation, perform one-way ANOVA across time points and
%    and electrodes and  save the largest F from each permutation (Fmax)
% 3. Compare Fobs for unpermuted data to distribution of Fmax to reject or 
%    fail to reject null for each time point an electrode

    global VERBLEVEL

    %Make sure there's only one factor
    assert(ndims(data) == 4);
    
    %Some useful numbers
    [n_electrodes, n_time_pts, n_conds, n_subs] = size(data);
    
    %Calculate degrees of freedom
    %(Always the same, so no point calculating in the loop)
    dfA   = n_conds - 1;
    dfBL  = n_subs - 1;
    dfRES = dfA * dfBL;
    
    %Perform n_perm permutations
    F_dist = NaN(n_perm, n_electrodes, n_time_pts);
    for i = 1:n_perm
        
        %Permute the data
        if i ==1
            perm_data = data;
        else
            for n = 1:n_subs
                perm_data(:, :, :, n) = data(:, :, randperm(size(data, 3)), n);
            end
        end
        
        %Calculate sums of squares
        SSyint = (sum(sum(perm_data, 3), 4).^2) / (n_conds * n_subs);
       %SSTO   = sum(sum(perm_data.^2, 3), 4) - SSyint;
        SSA    = (sum(sum(perm_data, 4).^2, 3) / n_subs) - SSyint;
        SSBL   = (sum(sum(perm_data, 3).^2, 4) / n_conds) - SSyint;
        SSRES  = sum(sum(perm_data.^2, 3), 4) - SSA - SSBL - SSyint;
        %assert(all(abs(SSTO(:) - (SSA(:) + SSBL(:) + SSRES(:))) < 1e-9));
        
        %Calculate F
        SSA(SSA < 1e-12) = 0; %Eliminates large F values that result from floating point error 
        F_dist(i, :, :) = (SSA/dfA) ./ (SSRES/dfRES);
        
        if VERBLEVEL
            if i == 1
                fprintf('Permutations completed: ')
            elseif i == n_perm
                fprintf('%d\n', i)
            elseif ~mod(i, 1000)
                fprintf('%d, ', i)
            end
        end
        
    end
    
    %Calculate Fobs, Fmax distribution, and Fmax critical value
    F_obs = reshape(F_dist(1, :, :), [n_electrodes, n_time_pts]);
    Fmax_dist = max(max(F_dist, [], 2), [], 3);
    Fmax_dist = sort(Fmax_dist);
    Fmax_crit = Fmax_dist(ceil((1-alpha) * length(Fmax_dist)));

    %Null hypothesis test
    h = F_obs > Fmax_crit;
    est_alpha = mean(Fmax_dist>Fmax_crit);
    if VERBLEVEL
        fprintf('Estimated alpha level is %f\n', est_alpha);
    end
    
    
    %Calculate p-value
    p = NaN(n_electrodes, n_time_pts);
    for e = 1:n_electrodes
        for t = 1:n_time_pts
            p(e,t) = mean(Fmax_dist >= F_obs(e,t));
        end
    end
    assert(isequal(h, p<=alpha));
    
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs;
    test_results.Fmax_crit = Fmax_crit;
    test_results.df = [dfA, dfRES];
    test_results.estimated_alpha = est_alpha;
    
end


function test_results = twoway_approx_int(data, n_perm, alpha)
%Use permutation of residuals method to calculate an approximate test of
%an interaction effect in a factorial design
% 1. Subtract main effects within each subject from all data points to obtain permutation residuals
% 2. Randomly permute all conditions within each subject across all time
%    points and electrodes
% 3. For each permutation, perform factorial ANOVA and save the largest F for the 
%    interaction effect across all time points and electrodes (Fmax)
% 4. Compare Fobs for unpermuted data to distribution of Fmax 
%    to reject or fail to reject null for each time point

    global VERBLEVEL

    %Make sure we're dealing with a two-way design
    assert(ndims(data) == 5);

    %Some useful numbers
    [n_electrodes, n_time_pts, n_conds_A, n_conds_B, n_subs] = size(data);

    %Subtract main effects within each subject so that the data is 
    %exchangeable under the null hypothesis for the interaction
    int_res = NaN(size(data));
    for p = 1:n_conds_A
        for q = 1:n_conds_B
            int_res(:,:, p, q, :) = data(:,:,p,q,:) ...
                                     - mean(data(:,:,p,:,:), 4) ...
                                     - mean(data(:,:,:,q,:), 3) ...
                                     + mean(mean(data, 3), 4);
        end
    end
    assert(abs(mean(int_res(:))) < 1e-9);
    
    %Calculate degrees of freedom
    %(Always the same, so no point calculating in the loop)
    dfBL     = n_subs - 1;
    dfA      = n_conds_A - 1;
    %dfAerr   = dfA * dfBL;
    dfB      = n_conds_B - 1;
    %dfBerr   = dfB * dfBL;
    dfAxB    = dfA * dfB;
    dfAxBerr = dfAxB * dfBL;
    %dfRES    = (num_subs - 1) * (num_conds_A * num_conds_B - 1);

    %Re-arrange data for permutation
    flat_data = reshape(int_res, [n_electrodes, n_time_pts, n_conds_A*n_conds_B, n_subs]);

    %Perform n_perm permutations
    F_dist = NaN(n_perm, n_electrodes, n_time_pts);
    flat_perm_data = NaN(size(flat_data));
    for i = 1:n_perm
        %Permute the data
        if i == 1
            perm_data = int_res;
        else
            for s = 1:n_subs
                flat_perm_data(:, :, :, s) = flat_data(:, :, randperm(size(flat_data, 3)), s);
            end
            perm_data = reshape(flat_perm_data, n_electrodes, n_time_pts, n_conds_A, n_conds_B, n_subs);
        end
        %Calculate F at each time point and electrode combination
        
        %Calculate sums of squares
        SSyint   = (sum(sum(sum(perm_data, 3), 4), 5).^2)/(n_conds_A*n_conds_B*n_subs);
        %SSTO     = sum(sum(sum(perm_data.^2, 3), 4), 5) - SSyint;
        SSA      = sum(sum(sum(perm_data, 4), 5).^2, 3)/(n_conds_B*n_subs) - SSyint;
        SSB      = sum(sum(sum(perm_data, 3), 5).^2, 4)/(n_conds_A*n_subs) - SSyint;
        SSBL     = sum(sum(sum(perm_data, 3), 4).^2, 5)/(n_conds_A*n_conds_B) - SSyint;
        SSAxB    = sum(sum(sum(perm_data, 5).^2, 3), 4)/n_subs - SSA - SSB -SSyint;
        SSAxBL   = sum(sum(sum(perm_data, 4).^2, 3), 5)/n_conds_B - SSA - SSBL - SSyint;
        SSBxBL   = sum(sum(sum(perm_data, 3).^2, 4), 5)/n_conds_A - SSB - SSBL - SSyint;
        SSAxBxBL = sum(sum(sum(perm_data.^2, 3), 4), 5) - SSA - SSB - SSBL - SSAxB - SSAxBL - SSBxBL - SSyint;
        %SSRES    = sum(sum(sum(perm_data.^2, 3), 4), 5) - SSA - SSB - SSBL - SSAxB - SSyint;

        %Doublechecking that the numbers match up
        %assert(all(SSRES - (SSAxBL + SSBxBL + SSAxBxBL) < 1e-9)); %SSRES is equal to its three subcomponents
        %assert(all(SSTO  - (SSRES + SSBL + SSA + SSB + SSAxB) < 1e-9)); %sums of squares add up

        %Calculate F
        SSAxB(SSAxB < 1e-12) = 0; %Eliminates large F values that result from floating point error 
        F_dist(i, :, :) = (SSAxB/dfAxB) ./ (SSAxBxBL/dfAxBerr);
        
        if VERBLEVEL
            if i == 1
                fprintf('Permutations completed: ')
            elseif i == n_perm
                fprintf('%d\n', i)
            elseif ~mod(i, 1000)
                fprintf('%d, ', i)
            end
        end

    end

    %Calculate Fobs, Fmax distribution, and Fmax critical value
    F_obs = reshape(F_dist(1,:,:), [n_electrodes, n_time_pts]);
    Fmax_dist = max(max(F_dist, [], 2), [], 3);
    Fmax_dist = sort(Fmax_dist);
    Fmax_crit = Fmax_dist(ceil((1-alpha) * length(Fmax_dist)));

    %Null hypothesis test
    h = F_obs > Fmax_crit;
    est_alpha = mean(Fmax_dist>Fmax_crit);
    if VERBLEVEL
        fprintf('Estimated alpha level is %f\n', est_alpha);
    end
    
    
    %Calculate p-value
    p = NaN(n_electrodes, n_time_pts);
    for e = 1:size(data, 1)
        for t = 1:size(data, 2)
            p(e,t) = mean(Fmax_dist >= F_obs(e,t));
        end
    end
    assert(isequal(h, p<=alpha));
    
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs;
    test_results.Fmax_crit = Fmax_crit;
    test_results.df = [dfAxB, dfAxBerr];
    test_results.estimated_alpha = est_alpha;

end


function test_results = threeway_approx_int(data, n_perm, alpha)
%Use permutation of residuals method to calculate an approximate test of
%an interaction effect in a factorial design
% 1. Subtract main effects within each subject from all data points to obtain 
%    interaction residuals
% 2. Subtract two-way interaction effect form residuals calculated above to
%    obtain three-way interaction residuals
% 3. Randomly permute all conditions within each subject across all time
%    points and electrodes
% 4. For each permutation, perform factorial ANOVA and save the largest F for the 
%    three-way interaction effect across all time points and electrodes (Fmax)
% 5. Compare Fobs for the unpermuted data to distribution of Fmax to reject 
%    or fail to reject null for each time point and electrode
    
    global VERBLEVEL

    %Make sure we're dealing with a two-way design
    assert(ndims(data) == 6);

    %Some useful numbers
    [n_electrodes, n_time_pts, n_conds_A, n_conds_B, n_conds_C, n_subs] = size(data);

    %Subtract main effects then two-way effects within each subject so that 
    %the data is exchangeable under the null hypothesis for the three-way
    %interaction
    int_res1 = NaN(size(data));
    int_res = NaN(size(data));
    %Subtract main effects
    for p = 1:n_conds_A
        for q = 1:n_conds_B
            for r = 1:n_conds_C
                int_res1(:,:,p,q,r,:) = data(:,:,p,q,r,:) ... 
                                        - mean(mean(data(:,:,p,:,:,:), 4), 5) ...
                                        - mean(mean(data(:,:,:,q,:,:), 3), 5) ...
                                        - mean(mean(data(:,:,:,:,r,:), 3), 4) ...
                                        + 2*mean(mean(mean(data, 3), 4), 5);
            end
        end
    end
    %Subtract two-way effects
    for p = 1:n_conds_A
        for q = 1:n_conds_B
            for r = 1:n_conds_C
                int_res(:,:,p,q,r,:) = int_res1(:,:,p,q,r,:) ... 
                                        - mean(int_res1(:,:,p,q,:,:), 5) ...
                                        - mean(int_res1(:,:,p,:,r,:), 4) ...
                                        - mean(int_res1(:,:,:,q,r,:), 3) ...
                                        + 2*mean(mean(mean(int_res1, 3), 4), 5);
            end
        end
    end
    assert(abs(mean(int_res(:))) < 1e-9)
    
    %Calculate degrees of freedom
    %(Always the same, so no point calculating in the loop)
    dfBL       = n_subs - 1;
    dfA        = n_conds_A - 1;
   %dfAerr     = dfA * dfBL;
    dfB        = n_conds_B - 1;
   %dfBerr     = dfB * dfBL;
    dfC        = n_conds_C - 1;
   %dfCerr     = dfC * dfBL;
   %dfAxB      = dfA * dfB;
   %dfAxBerr   = dfAxB * dfBL;
   %dfAxC      = dfA * dfC;
   %dfAxCerr   = dfAxC * dfBL;
   %dfBxC      = dfB * dfC;
   %dfBxCerr   = dfBxC * dfBL;
    dfAxBxC    = dfA * dfB * dfC;
    dfAxBxCerr = dfAxBxC * dfBL;
   %dfRES      = (n_subs - 1) * (n_conds_A * n_conds_B * n_conds_C - 1);

    %Re-arrange data for permutation
    flat_data = reshape(int_res, [n_electrodes, n_time_pts, n_conds_A*n_conds_B*n_conds_C, n_subs]);

    %Perform n_perm permutations
    F_dist = NaN(n_perm, n_electrodes, n_time_pts);
    flat_perm_data = NaN(size(flat_data));
    for i = 1:n_perm

        %Permute the data
        if i == 1
            perm_data = int_res;
        else
            for s = 1:n_subs
                flat_perm_data(:, :, :, s) = flat_data(:, :, randperm(size(flat_data, 3)), s);
            end
            perm_data = reshape(flat_perm_data, size(int_res));
        end
        
        %Calculate F at each time point and electrode combination

        %Calculate sums of squares
        SSyint     = (sum(sum(sum(sum(perm_data, 3), 4), 5), 6).^2)/(n_conds_A*n_conds_B*n_conds_C*n_subs);
       %SSTO       = sum(sum(sum(sum(perm_data.^2, 3), 4), 5), 6) - SSyint;
        SSA        = sum(sum(sum(sum(perm_data, 4), 5), 6).^2, 3)/(n_conds_B*n_conds_C*n_subs) - SSyint;
        SSB        = sum(sum(sum(sum(perm_data, 3), 5), 6).^2, 4)/(n_conds_A*n_conds_C*n_subs) - SSyint;
        SSC        = sum(sum(sum(sum(perm_data, 3), 4), 6).^2, 5)/(n_conds_A*n_conds_B*n_subs) - SSyint;
        SSBL       = sum(sum(sum(sum(perm_data, 3), 4), 5).^2, 6)/(n_conds_A*n_conds_B*n_conds_C) - SSyint;
        SSAxB      = sum(sum(sum(sum(perm_data, 5), 6).^2, 3), 4)/(n_conds_C*n_subs) - SSA - SSB -SSyint;
        SSAxC      = sum(sum(sum(sum(perm_data, 4), 6).^2, 3), 5)/(n_conds_B*n_subs) - SSA - SSC -SSyint;
        SSBxC      = sum(sum(sum(sum(perm_data, 3), 6).^2, 4), 5)/(n_conds_A*n_subs) - SSB - SSC -SSyint;
        SSAxBxC    = sum(sum(sum(sum(perm_data, 6).^2, 3), 4), 5)/n_subs - SSA - SSB - SSC - SSAxB - SSAxC - SSBxC - SSyint;
        SSAxBL     = sum(sum(sum(sum(perm_data, 4), 5).^2, 3), 6)/(n_conds_B*n_conds_C) - SSA - SSBL - SSyint;
        SSBxBL     = sum(sum(sum(sum(perm_data, 3), 5).^2, 4), 6)/(n_conds_A*n_conds_C) - SSB - SSBL - SSyint;
        SSCxBL     = sum(sum(sum(sum(perm_data, 3), 4).^2, 5), 6)/(n_conds_A*n_conds_B) - SSC - SSBL - SSyint;
        SSAxBxBL   = sum(sum(sum(sum(perm_data, 5).^2, 3), 4), 6)/n_conds_C - SSA - SSB - SSBL - SSAxB - SSAxBL - SSBxBL - SSyint;
        SSAxCxBL   = sum(sum(sum(sum(perm_data, 4).^2, 3), 5), 6)/n_conds_B - SSA - SSC - SSBL - SSAxC - SSAxBL - SSCxBL - SSyint;
        SSBxCxBL   = sum(sum(sum(sum(perm_data, 3).^2, 4), 5), 6)/n_conds_A - SSB - SSC - SSBL - SSBxC - SSBxBL - SSCxBL - SSyint;
        SSAxBxCxBL = sum(sum(sum(sum(perm_data.^2, 3), 4), 5), 6) - SSA - SSB - SSC - SSBL - SSAxB - SSAxC - SSBxC - SSAxBL - SSBxBL - SSCxBL - SSAxBxC - SSAxBxBL - SSAxCxBL - SSBxCxBL - SSyint;
       %SSRES      = sum(sum(sum(sum(perm_data.^2, 3), 4), 5), 6) - SSA - SSB - SSC - SSBL - SSAxB - SSAxC - SSBxC - SSAxBxC - SSyint;

        %Doublechecking that the numbers match up
        %assert(all(SSRES - (SSAxBL + SSBxBL + SSCxBL + SSAxBxBL + SSAxCxBL + SSBxCxBL + SSAxBxCxBL) < 1e-9)); %SSRES is equal to its three subcomponents
        %assert(all(SSTO  - (SSRES + SSBL + SSA + SSB + SSC + SSAxB + SSAxC + SSBxC + SSAxBxC) < 1e-9)); %sums of squares add up

        SSAxBxC(SSAxBxC < 1e-12) = 0; %Eliminates large F values that result from floating point error 
        F_dist(i, :, :) = (SSAxBxC/dfAxBxC) ./ (SSAxBxCxBL/dfAxBxCerr);
        
        if VERBLEVEL
            if i == 1
                fprintf('Permutations completed: ')
            elseif i == n_perm
                fprintf('%d\n', i)
            elseif ~mod(i, 1000)
                fprintf('%d, ', i)
            end
        end

    end
    
    %Calculate Fobs, Fmax distribution, and Fmax critical value
    F_obs = reshape(F_dist(1,:,:), [n_electrodes, n_time_pts]);
    Fmax_dist = max(max(F_dist, [], 2), [], 3);
    Fmax_dist = sort(Fmax_dist);
    Fmax_crit = Fmax_dist(ceil((1-alpha) * length(Fmax_dist)));

    %Null hypothesis test
    h = F_obs > Fmax_crit;
    est_alpha = mean(Fmax_dist>Fmax_crit);
    if VERBLEVEL
        fprintf('Estimated alpha level is %f\n', est_alpha);
    end
    
    %Calculate p-value
    p = NaN(n_electrodes, n_time_pts);
    for e = 1:size(data, 1)
        for t = 1:size(data, 2)
            p(e,t) = mean(Fmax_dist >= F_obs(e,t));
        end
    end
    assert(isequal(h, p<=alpha));
    
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs;
    test_results.Fmax_crit = Fmax_crit;
    test_results.df = [dfAxBxC, dfAxBxCerr];
    test_results.estimated_alpha = est_alpha;

end