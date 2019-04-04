%Calculate Fmax corrected F-test.
%
%EXAMPLE USAGE
% >> test_results = calc_Fmax(data, [], [3, 4], 1e4, 0.05)
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
%
%OUTPUT
% test_results - A struct with results of the Fmax test
%
%
%VERSION DATE: 4 April 2019
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 

%Copyright (c) 2017-2019, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function test_results = calc_Fmax(data, cond_subs, dims, n_perm, alpha)
     
    %%% Calculate ANOVA %%%
    
    if ~isempty(cond_subs) && ~isequal(cond_subs, 0) && length(cond_subs) > 1
        [F_obs, F_dist, df_effect, df_res, exact_test] = perm_spANOVA(data, cond_subs, dims, n_perm);
    else
        [F_obs, F_dist, df_effect, df_res, exact_test] = perm_rbANOVA(data, dims, n_perm);
    end
    
    %%% Calculate Fmax correction %%%
    
    [h, p, Fmax_crit, est_alpha] = Fmax_corr(F_obs, F_dist, alpha);
    
    %%% Output %%%
    
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs;
    test_results.Fmax_crit = Fmax_crit;
    test_results.df = [df_effect, df_res];
    if exact_test
        test_results.estimated_alpha = est_alpha;
        test_results.exact_test = true;
    else
        test_results.estimated_alpha = NaN;
        test_results.exact_test = false;
    end

end
