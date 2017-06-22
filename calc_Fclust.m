%Calculate F cluster mass permutation effect.
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
% chan_hood     - A 2D symmetric binary matrix that indicates
%                 which channels are considered neighbors of other 
%                 channels. E.g., if chan_hood(2,10)=1, then Channel 2 
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
%VERSION DATE: 22 June 2017
%AUTHOR: Eric Fields, Tufts University (Eric.Fields@tufts.edu)
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 
%This is a beta version of this software. It needs additional testing 
%and SHOULD NOT be considered error free.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

%%%%%%%%%%%%%%%%%%%  REVISION LOG   %%%%%%%%%%%%%%%%%%%
% 11/28/16  - Moved from FclustGND. See FclustGND revision log for earlier
%             versions
% 4/17/17   - Added estimated alpha
% 6/12/17   - Added verblevel reports
% 6/13/17   - null test for clusters are now logicals
% 6/22/17   - major re-organization of code to reduce repeated code

function test_results = calc_Fclust(data, dims, n_perm, alpha, int_method, chan_hood, thresh_p)

    global VERBLEVEL

    %Some useful numbers
    n_electrodes = size(data, 1);
    n_time_pts   = size(data, 2);
    
    %% Calcualte ANOVA
    
    %Eliminate factors not involved in this effect by averaging or
    %reduce exact interaction to one-way via subtraction
    reduced_data = reduce_data(data, dims, int_method);
    
    %Calculate the ANOVA (F-obs and the permutation distribution)
    [F_dist, df_effect, df_res] = perm_rbANOVA(reduced_data, n_perm);
    
    
    %% Find clusters and cluster ditribution
    
    %Find clusters and cluster mass distribution
    thresh_F = finv(1-thresh_p, df_effect, df_res); 
    clust_mass_dist = zeros(n_perm, 1);
    for i = 1:n_perm
        F = reshape(F_dist(i, :, :), n_electrodes, n_time_pts);
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
    
    %% Output
    
    clust_info = struct('null_test', null_test, ...
                        'pval', clust_pval, ...
                        'clust_mass', clust_mass_obs, ...
                        'clust_ids', clust_ids_obs);
    
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs; 
    test_results.df = [df_effect, df_res];
    test_results.clust_info = clust_info;
    test_results.estimated_alpha = est_alpha;

end
