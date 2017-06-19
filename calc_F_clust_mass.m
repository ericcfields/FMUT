%Calculate F cluster mass permutation effect. This function is called by FclustGND.
%For more information, see:
%>> help FclustGND
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
% n_perm        - Number of permutations to use to calculate null distribution
% alpha         - Alpha level of the test
% int_method    - A string that should be either 'exact' or 'approximate'.
%                 If 'exact', the method of restricted permutations will
%                 be used to conduct a test that controls the Type I error
%                 rate at alpha (assuming enough permutations). 
%                 If 'approximate', the method of permutation of residuals 
%                 will be used to conduct a test with Type I error rate 
%                 asymptotic to alpha as noise decreases and/or number of 
%                 subjects increases.
% chan_hood      - A 2D symmetric binary matrix that indicates
%                  which channels are considered neighbors of other 
%                  channels. E.g., if chan_hood(2,10)=1, then Channel 2 
%                  and Channel 10 are nieghbors. You can produce a 
%                  chan_hood matrix using the function spatial_neighbors.m. 
% thresh_p       - The test-wise p-value threshold for cluster inclusion. If
%                  a channel/time-point has a F-value that corresponds to an
%                  uncorrected p-value greater than thresh_p, it is assigned
%                  a p-value of 1 and not considered for clustering.
%
%OUTPUT
% test_results - A struct with results of the cluster mass test
%
%
%VERSION DATE: 13 June 2017
%AUTHOR: Eric Fields, Tufts University (Eric.Fields@tufts.edu)
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 
%This is a beta version of this software. It needs additional testing 
%and SHOULD NOT be considered error free.

%Copyright (c) 2016, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.
%This function incorporates some code from the Mass Univariate Toolbox, 
%Copyright (c) 2015, David Groppe: https://github.com/dmgroppe/Mass_Univariate_ERP_Toolbox/blob/master/LICENSE

%%%%%%%%%%%%%%%%%%%  REVISION LOG   %%%%%%%%%%%%%%%%%%%
% 11/28/16  - Moved from FclustGND. See FclustGND revision log for earlier
%             versions
% 4/17/17   - Added estimated alpha
% 6/12/17   - Added verblevel reports
% 6/13/17   - null test for clusters are now logicals


function test_results = calc_F_clust_mass(data, dims, n_perm, alpha, int_method, chan_hood, thresh_p)
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
        test_results = oneway(reduced_data, n_perm, alpha, chan_hood, thresh_p);
    elseif length(dims) == 2
        test_results = twoway_approx_int(reduced_data, n_perm, alpha, chan_hood, thresh_p);
    elseif length(dims) == 3
        test_results = threeway_approx_int(reduced_data, n_perm, alpha, chan_hood, thresh_p);
    end

end


function test_results = oneway(data, n_perm, alpha, chan_hood, thresh_p)
    
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
    
    %Threshold for clutster inclusion
    thresh_F = finv(1-thresh_p, dfA, dfRES); 
    
    %Perform n_perm permutations
    clust_mass_dist = zeros(n_perm, 1);
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
        F = (SSA/dfA) ./ (SSRES/dfRES);
        
        %Find clusters in F
        [clust_ids, n_clust] = find_clusters(F, thresh_F, chan_hood, 1);
        
        %Save observed data
        if i == 1
            F_obs = F;
            clust_ids_obs = clust_ids;
            n_clust_obs = n_clust;
        end
        
        %Find largest cluster mass for this permutation
        for c = 1:n_clust,
            use_mass = sum(F(clust_ids == c));
            if use_mass > clust_mass_dist(i),
                clust_mass_dist(i) = use_mass;
            end
        end
        
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
    
    %Get cluster mass critical value
    clust_mass_dist = sort(clust_mass_dist);
    clust_mass_crit = clust_mass_dist(ceil((1-alpha) * n_perm));
    
    %Find cluster mass for each observed cluster and compare to critical
    %value
    clust_mass_obs = NaN(1, n_clust_obs);
    clust_pval = NaN(1, n_clust_obs);
    null_test = false(1, n_clust_obs);
    p = ones(n_electrodes, n_time_pts);
    for c = 1:n_clust_obs
        clust_mass_obs(c) = sum(F_obs(clust_ids_obs == c));
        if clust_mass_obs(c) > clust_mass_crit
            null_test(c) = true;
        end
        clust_pval(c) = mean(clust_mass_dist >= clust_mass_obs(c));
        p(clust_ids_obs == c) = clust_pval(c);
    end
    h = (p <= alpha);
    est_alpha = mean(clust_mass_dist>clust_mass_crit);
    if VERBLEVEL
        fprintf('Estimated alpha level is %f\n', est_alpha);
    end
    
    
    clust_info = struct('null_test', null_test, ...
                        'pval', clust_pval, ...
                        'clust_mass', clust_mass_obs, ...
                        'clust_ids', clust_ids_obs);
    
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs; 
    test_results.df = [dfA, dfRES];
    test_results.clust_info = clust_info;
    test_results.estimated_alpha = est_alpha;

end


function test_results = twoway_approx_int(data, n_perm, alpha, chan_hood, thresh_p)
%Use permutation of residuals method to calculate an approximate test of
%an interaction effect in a factorial design

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
    assert(abs(mean(int_res(:))) < 1e-9)
    
    %Calculate degrees of freedom
    %(Always the same, so no point calculating in the loop)
    %dfBL     = n_subs - 1;
    %dfA      = n_conds_A - 1;
    %dfAerr   = dfA * dfBL;
    %dfB      = n_conds_B - 1;
    %dfBerr   = dfB * dfBL;
    dfAxB    = (n_conds_A - 1) * (n_conds_B - 1);
    dfAxBerr = dfAxB * (n_subs - 1);
    %dfRES    = (num_subs - 1) * (num_conds_A * num_conds_B - 1);
    
    %Threshold for clutster inclusion
    thresh_F = finv(1-thresh_p, dfAxB, dfAxBerr); 

    %Re-arrange data for permutation
    flat_data = reshape(int_res, [n_electrodes, n_time_pts, n_conds_A*n_conds_B, n_subs]);

    %Perform n_perm permutations
    clust_mass_dist = zeros(n_perm, 1);
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
        F = (SSAxB/dfAxB) ./ (SSAxBxBL/dfAxBerr);
        
        %Find clusters in F
        [clust_ids, n_clust] = find_clusters(F, thresh_F, chan_hood, 1);
        
        %Save observed data
        if i == 1
            F_obs = F;
            clust_ids_obs = clust_ids;
            n_clust_obs = n_clust;
        end
        
        %Find largest cluster mass for this permutation
        for c = 1:n_clust,
            use_mass = sum(F(clust_ids == c));
            if use_mass > clust_mass_dist(i),
                clust_mass_dist(i) = use_mass;
            end
        end
        
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
    
        %Get cluster mass critical value
    clust_mass_dist = sort(clust_mass_dist);
    clust_mass_crit = clust_mass_dist(ceil((1-alpha) * n_perm));
    
    %Find cluster mass for each observed cluster and compare to critical
    %value
    clust_mass_obs = NaN(1, n_clust_obs);
    clust_pval = NaN(1, n_clust_obs);
    null_test = zeros(1, n_clust_obs);
    p = ones(n_electrodes, n_time_pts);
    for c = 1:n_clust_obs
        clust_mass_obs(c) = sum(F_obs(clust_ids_obs == c));
        if clust_mass_obs(c) > clust_mass_crit
            null_test(c) = 1;
        end
        clust_pval(c) = mean(clust_mass_dist >= clust_mass_obs(c));
        p(clust_ids_obs == c) = clust_pval(c);     
    end
    h = (p <= alpha);
    est_alpha = mean(clust_mass_dist>clust_mass_crit);
    if VERBLEVEL
        fprintf('Estimated alpha level is %f\n', est_alpha);
    end
    
    
    clust_info = struct('null_test', null_test, ...
                        'pval', clust_pval, ...
                        'clust_mass', clust_mass_obs, ...
                        'clust_ids', clust_ids_obs);
    
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs; 
    test_results.df = [dfAxB, dfAxBerr];
    test_results.clust_info = clust_info;
    test_results.estimated_alpha = est_alpha;
    
end


function test_results = threeway_approx_int(data, n_perm, alpha, chan_hood, thresh_p)
%Use permutation of residuals method to calculate an approximate test of
%an interaction effect in a factorial design

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
   %dfBL       = n_subs - 1;
   %dfA        = n_conds_A - 1;
   %dfAerr     = dfA * dfBL;
   %dfB        = n_conds_B - 1;
   %dfBerr     = dfB * dfBL;
   %dfC        = n_conds_C - 1;
   %dfCerr     = dfC * dfBL;
   %dfAxB      = dfA * dfB;
   %dfAxBerr   = dfAxB * dfBL;
   %dfAxC      = dfA * dfC;
   %dfAxCerr   = dfAxC * dfBL;
   %dfBxC      = dfB * dfC;
   %dfBxCerr   = dfBxC * dfBL;
    dfAxBxC    = (n_conds_A - 1) * (n_conds_B - 1) * (n_conds_C - 1);
    dfAxBxCerr = dfAxBxC * (n_subs - 1);
   %dfRES      = (n_subs - 1) * (n_conds_A * n_conds_B * n_conds_C - 1);
   
    %Threshold for clutster inclusion
    thresh_F = finv(1-thresh_p, dfAxBxC, dfAxBxCerr); 

    %Re-arrange data for permutation
    flat_data = reshape(int_res, [n_electrodes, n_time_pts, n_conds_A*n_conds_B*n_conds_C, n_subs]);

    %Perform n_perm permutations
    clust_mass_dist = zeros(n_perm, 1);
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
        F = (SSAxBxC/dfAxBxC) ./ (SSAxBxCxBL/dfAxBxCerr);

        %Find clusters in F
        [clust_ids, n_clust] = find_clusters(F, thresh_F, chan_hood, 1);
        
        %Save observed data
        if i == 1
            F_obs = F;
            clust_ids_obs = clust_ids;
            n_clust_obs = n_clust;
        end
        
        %Find largest cluster mass for this permutation
        for c = 1:n_clust,
            use_mass = sum(F(clust_ids == c));
            if use_mass > clust_mass_dist(i),
                clust_mass_dist(i) = use_mass;
            end
        end
        
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
    
    %Get cluster mass critical value
    clust_mass_dist = sort(clust_mass_dist);
    clust_mass_crit = clust_mass_dist(ceil((1-alpha) * n_perm));
    
    %Find cluster mass for each observed cluster and compare to critical
    %value
    clust_mass_obs = NaN(1, n_clust_obs);
    clust_pval = NaN(1, n_clust_obs);
    null_test = zeros(1, n_clust_obs);
    p = ones(n_electrodes, n_time_pts);
    for c = 1:n_clust_obs
        clust_mass_obs(c) = sum(F_obs(clust_ids_obs == c));
        if clust_mass_obs(c) > clust_mass_crit
            null_test(c) = 1;
        end
        clust_pval(c) = mean(clust_mass_dist >= clust_mass_obs(c));
        p(clust_ids_obs == c) = clust_pval(c);     
    end
    h = (p <= alpha);
    est_alpha = mean(clust_mass_dist>clust_mass_crit);
    if VERBLEVEL
        fprintf('Estimated alpha level is %f\n', est_alpha);
    end
    
    
    clust_info = struct('null_test', null_test, ...
                        'pval', clust_pval, ...
                        'clust_mass', clust_mass_obs, ...
                        'clust_ids', clust_ids_obs);
    
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs;
    test_results.df = [dfAxBxC, dfAxBxCerr];
    test_results.clust_info = clust_info;
    test_results.estimated_alpha = est_alpha;
    
end