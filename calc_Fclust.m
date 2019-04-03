%Calculate F-test cluster mass permutation effect.
%
%EXAMPLE USAGE
% >> test_results = calc_Fclust(data, [], [3, 4], 1e4, 0.05, chan_hood, 0.05)
%
%REQUIRED INPUTS
% data          - An electrode x time points x conditions x subjects array of ERP
%                 data. Array will vary in number of dimensions based on how many
%                 factors there are
% cond_subs     - Array giving the number of subjects in each condition of
%                 the between subjects factor. For example, if cond_subs is
%                 [8, 9], then there should be 17 subjects with the first 8
%                 being in condition A and the next 9 being in condition B.
%                 For fully within-subjects designs cond_subs = []
% dims          - Dimensions of the data array involved in the effect to be
%                 calculated. For example, if data is an electrode x time points
%                 x Factor A x Factor B x subjects array and you want to
%                 calculated the main effect of A, dims = 3. If you want to
%                 calculate the AxB interaciton, dims  = [3, 4].
% n_perm        - Number of permutations to use to calculate the null distribution
% alpha         - Family-wise alpha level of the test
% chan_hood     - A 2D symmetric binary matrix that indicates
%                 which channels are considered neighbors of other 
%                 channels. E.g., if chan_hood(2,10)==1, then Channel 2 
%                 and Channel 10 are nieghbors. You can produce a 
%                 chan_hood matrix using the function spatial_neighbors.m. 
% thresh_p      - The test-wise p-value threshold for cluster inclusion. If
%                 a channel/time-point has a F-value that corresponds to an
%                 uncorrected p-value greater than thresh_p, it is assigned
%                 a p-value of 1 and not considered for clustering.
%
%OUTPUT
% test_results - A struct with results of the cluster mass test
%
%
%VERSION DATE: 3 April 2019
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function test_results = calc_Fclust(data, cond_subs, dims, n_perm, alpha, chan_hood, thresh_p)

    %Some useful numbers
    n_electrodes = size(data, 1);
    n_time_pts   = size(data, 2);
    
    %%% Calculate ANOVA %%%
    
    if ~isempty(cond_subs) && ~isequal(cond_subs, 0) && length(cond_subs) > 1
        [F_dist, df_effect, df_res, exact_test] = perm_spANOVA(data, cond_subs, dims, n_perm);
    else
        [F_dist, df_effect, df_res, exact_test] = perm_rbANOVA(data, dims, n_perm);
    end
    F_obs = reshape(F_dist(1, :, :), [n_electrodes, n_time_pts]);
    
    %%% Calculate cluster correction %%%
    
    thresh_F = finv(1-thresh_p, df_effect, df_res); 
    [h, p, clust_info, est_alpha] = Fclust_corr(F_obs, F_dist, alpha, chan_hood, thresh_F);
    
    %%% Output %%%
    
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs; 
    test_results.df = [df_effect, df_res];
    test_results.clust_info = clust_info;
    if exact_test
        test_results.estimated_alpha = est_alpha;
        test_results.exact_test = true;
    else
        test_results.estimated_alpha = NaN;
        test_results.exact_test = false;
    end
    
end

function [h, p, clust_info, est_alpha] = Fclust_corr(F_obs, F_dist, alpha, chan_hood, thresh_F)
%Calculate cluster correction given the permutation F_distribution

    global VERBLEVEL

    %% Find clusters and cluster ditribution
    
    if VERBLEVEL
        fprintf('Calculating clusters . . . ')
    end
    
    %Some useful numbers
    [n_perm, n_electrodes, n_time_pts] = size(F_dist);
    
    %Find clusters and cluster mass distribution
    clust_mass_dist = zeros(n_perm, 1);
    for i = 1:n_perm
        F = reshape(F_dist(i, :, :), [n_electrodes, n_time_pts]);
        %Find clusters in F
        [clust_ids, n_clust] = find_clusters(F, thresh_F, chan_hood, 1);
        %Save observed data
        if i == 1
            clust_ids_obs = clust_ids;
            n_clust_obs = n_clust;
        end
        %Find largest cluster mass for this permutation
        for c = 1:n_clust
            use_mass = sum(F(clust_ids == c));
            if use_mass > clust_mass_dist(i)
                clust_mass_dist(i) = use_mass;
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
        fprintf('DONE\n');
        fprintf('Estimated alpha level is %f\n', est_alpha);
    end
    
    clust_info = struct('null_test', null_test, ...
                        'pval', clust_pval, ...
                        'clust_mass', clust_mass_obs, ...
                        'clust_ids', clust_ids_obs);

end
