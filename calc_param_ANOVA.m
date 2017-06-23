%Calculate parametric ANOVA at each time point and electrode. Results are
%not multiple comparison corrected.
%
%REQUIRED INPUTS
% data          - An electrode x time points x conditions x subjects array of ERP
%                 data. Array will vary in number of dimensions based on how many
%                 factors there are.
% dims          - Dimensions of the data array involved in the effect to be
%                 calculated. For example, if data is an electrode x time points
%                 x Factor A x Factor B x subjects array and you want to
%                 calculated the main effect of A, dims = 3. If you want to
%                 calculate the AxB interaciton, dims  = [3, 4].
% alpha         - Test-wise alpha level
%
%OUTPUT
% test_results - A struct with results of the mass univariate ANOVA
%
%
%VERSION DATE: 23 June 2017
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


function test_results = calc_param_ANOVA(data, dims, alpha)

    %Some useful numbers
    n_electrodes = size(data, 1);
    n_time_pts   = size(data, 2);

    %Average across factors not involved in this effect
    reduced_data = reduce_data(data, dims, 'none');
    
    %Calculate ANOVA with single permutation to get F-obs
    [F_dist, df_effect, df_res] = perm_rbANOVA(reduced_data, 1);
    F_obs = reshape(F_dist, n_electrodes, n_time_pts);
    
    %Calculate F_crit
    F_crit = finv(.95, df_effect, df_res);

    %Null hypothesis test
    h = F_obs > F_crit;
    
    %Calculate p-value
    p = 1 - fcdf(F_obs, df_effect, df_res);
    assert(isequal(h, p<=alpha));
    
    %Create results struct
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs;
    test_results.F_crit = F_crit;
    test_results.df = [df_effect, df_res];

end
