%Calculate Fmax corrected F-test.
%For more inforamtion see: >> help FmaxGND
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
% n_perm        - Number of permutations to use to calculate the null distribution
% alpha         - Family-wise alpha level of the test
%
%OUTPUT
% test_results - A struct with results of the Fmax test
%
%
%VERSION DATE: 13 July 2017
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
% 7/13/17   - Updated for elimination of int_method input


function test_results = calc_Fmax(data, dims, n_perm, alpha)
    
    global VERBLEVEL

    %Some useful numbers
    n_electrodes = size(data, 1);
    n_time_pts   = size(data, 2);
    
    
    %% Calculate ANOVA
    
    %Eliminate factors not involved in this effect by averaging or
    %reduce exact interaction to one-way via subtraction
    reduced_data = reduce_data(data, dims);
    
    %Calculate the ANOVA (F-obs and the permutation distribution)
    [F_dist, df_effect, df_res] = perm_rbANOVA(reduced_data, n_perm);
    F_obs = reshape(F_dist(1, :, :), [n_electrodes, n_time_pts]);
    
    
    %% Calculate Fmax correction
    
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
    
    
    %% Output
    
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs;
    test_results.Fmax_crit = Fmax_crit;
    test_results.df = [df_effect, df_res];
    if ndims(reduced_data) == 4
        test_results.estimated_alpha = est_alpha;
        test_results.exact_test = true;
    else
        test_results.estimated_alpha = NaN;
        test_results.exact_test = false;
    end

end

