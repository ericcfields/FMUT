%Calculate parametric ANOVA at each time point and electrode with various
%multiple comparisons corrections available.
%
%REQUIRED INPUTS
% data          - An electrode x time points x conditions x subjects array of ERP
%                 data. Array will vary in number of dimensions based on how many
%                 factors there are.
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
% alphaORq      - Test-wise alpha level (no correction), family-wise alpha
%                 level (bonferroni or sidak correction) or False Discovery
%                 Rate (bh, bk, or bky corrections)
%
%OPTIONAL INPUTS
% correction    - A string indicating a multiple comparisons correction.
%                 Options are 'none', 'bonferroni', 'sidak', 'bh' (classic 
%                 Benjamini & Hochberg 1995) procedure, 'by' (Benjamini & 
%                 Yekutieli, 2001), or 'bky' (Benjamini, Krieger, &
%                 Yekutieli, 2006). {default: 'none'}.
%
%OUTPUT
% test_results - A struct with results of the mass univariate ANOVA
%
%
%VERSION DATE: 24 July 2017
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 
%This is a beta version of this software. It needs additional testing 
%and SHOULD NOT be considered error free.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

%%%%%%%%%%%%%%%%%%%  REVISION LOG   %%%%%%%%%%%%%%%%%%%
% 4/9/17  - Added documentation
% 6/22/17 - Updated to eliminate repeated code
% 6/23/17 - Added multiple comparisons corrections
% 7/13/17 - Updated to reflect that ANOVA sub-functions no longer require
%           int_method input
% 7/14/17 - Now handle between subjects factors
% 7/24/17 - Moved reduced_data to ANOVA functions

function test_results = calc_param_ANOVA(data, cond_subs, dims, alphaORq, correction)

    if nargin < 5
        correction = 'none';
    elseif ~any(strcmpi(correction, {'none', 'bonferroni', 'sidak', 'bh', 'by', 'bky'}))
        error('correction must be ''none'', ''bonferroni'', ''sidak'', ''bh'', ''by'', or ''bky''');
    end

    %Some useful numbers
    n_electrodes = size(data, 1);
    n_time_pts   = size(data, 2);
    
    %% Calculate ANOVA

    %Calculate the ANOVA (F-obs and the permutation distribution)
    if ~isempty(cond_subs) && ~isequal(cond_subs, 0) && length(cond_subs) > 1
        [F_dist, df_effect, df_res] = perm_spANOVA(data, cond_subs, dims, 1);
    else
        [F_dist, df_effect, df_res] = perm_rbANOVA(data, dims, 1);
    end
    F_obs = reshape(F_dist(1, :, :), [n_electrodes, n_time_pts]);
    uncorr_p = 1 - fcdf(F_obs, df_effect, df_res);
    
    
    %% Calculate F-crit, p-values, and null test with appropriate correction
    
    switch correction
        case 'none'
            F_crit = finv(1-alphaORq, df_effect, df_res);
            p = uncorr_p;
            h = F_obs > F_crit;
        case 'bonferroni'
            F_crit = finv(1 - (alphaORq/(n_electrodes * n_time_pts)), ... 
                          df_effect, df_res);
            p = uncorr_p * (n_electrodes * n_time_pts);
            p(p>1) = 1;
            h = F_obs > F_crit;
        case 'sidak'
            F_crit = finv((1-alphaORq).^(1/(n_electrodes*n_time_pts)), ... 
                          df_effect, df_res);
            p = 1 - (1 - uncorr_p) .^ (n_electrodes * n_time_pts);
            h = F_obs > F_crit;
        case 'bh'
            [h, crit_p, ~, p] = fdr_bh(uncorr_p, alphaORq, 'pdep', 'no');
        case 'by'
            [h, crit_p, ~, p] = fdr_bh(uncorr_p, alphaORq, 'dep', 'no');
        case 'bky'
            [h, crit_p] = fdr_bky(uncorr_p, alphaORq, 'no');
            p = NaN;   
    end
    if any(strcmpi(correction, {'bh', 'by', 'bky'}))
        if crit_p == 0
            F_crit = NaN;
        else
            F_crit = finv(1-crit_p, df_effect, df_res);
        end
    end
    if ~strcmpi(correction ,'bky')
        assert(isequal(h, p<=alphaORq));
    end

    %Create results struct
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs;
    test_results.F_crit = F_crit;
    test_results.df = [df_effect, df_res];

end
