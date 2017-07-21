%Calculate Fmax corrected F-test.
%For more inforamtion see: 
%>> help FmaxGND
%and
%>> help FmaxGRP
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
%VERSION DATE: 21 July 2017
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
% 11/28/16  - Moved from FmaxGND. See FmaxGND revision log for information
%             on earlier versions
% 6/12/17   - Added estimated alpha; added verblevel reports
% 7/10/17   - Now supports between subjects factors
% 7/13/17   - Updated for elimination of int_method input
% 7/21/17   - Added step-down procedure


function test_results = calc_Fmax(data, cond_subs, dims, n_perm, alpha, step_down)
    
    global VERBLEVEL
    
    if nargin < 6
        step_down = false;
    end

    %Some useful numbers
    n_electrodes = size(data, 1);
    n_time_pts   = size(data, 2);
    
    
    %% Calculate ANOVA
    
    %Eliminate factors not involved in this effect and reduce interactions
    %via subtraction
    [reduced_data, new_dims] = reduce_data(data, dims);
    
    %Calculate the ANOVA (F-obs and the permutation distribution)
    if ~isempty(cond_subs) && ~isequal(cond_subs, 0) && length(cond_subs) > 1
        if ndims(reduced_data) == 3
            [F_dist, df_effect, df_res] = perm_crANOVA(reduced_data, cond_subs, n_perm);
            exact_test = true;
        else
            [F_dist, df_effect, df_res, exact_test] = perm_spANOVA(reduced_data, cond_subs, new_dims, n_perm);
        end
    else
        [F_dist, df_effect, df_res, exact_test] = perm_rbANOVA(reduced_data, n_perm);
    end
    F_obs = reshape(F_dist(1, :, :), [n_electrodes, n_time_pts]);
    
    
    %% Calculate Fmax correction
    
    if step_down
        
        %Null distributions and critial values
        flat_F_dist = reshape(F_dist, [n_perm, n_electrodes*n_time_pts]);
        step_down_dist = sort(sort(flat_F_dist, 2, 'descend'), 1);
        Fmax_crit = step_down_dist(ceil((1-alpha) * size(step_down_dist, 1)), :)';
        
        %Get ordered F-values
        [sorted_F_obs, idx] = sort(F_obs(:), 'descend');
        
        %Null hypothesis test
        h = sorted_F_obs > Fmax_crit;
        stop_point = find(~h, 1);
        h(stop_point:end) = false;
        h = reshape(h(idx), [n_electrodes, n_time_pts]);
        est_alpha = mean(step_down_dist(:,1) > Fmax_crit(1));
        if VERBLEVEL
            fprintf('Estimated alpha level is %f\n', est_alpha);
        end
        
        %Calculate p-values
        p = NaN(size(sorted_F_obs));
        for i = 1:length(sorted_F_obs)
            p(i) = mean(step_down_dist(:, i) > sorted_F_obs(i));
        end
        p(stop_point:end) = NaN;
        p = reshape(p(idx), [n_electrodes, n_time_pts]);

        assert(isequal(h, p<=alpha));
        
    else
    
        %Calculate Fmax distribution and Fmax critical value
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
        
    end
    
    
    %% Output
    
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

