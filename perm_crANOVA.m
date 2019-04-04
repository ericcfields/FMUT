%Calculate F-observed and the empirical F-distribution for a one-way
%between subjects ANOVA
%
%EXAMPLE USAGE
% >> [F_obs, F_dist, df_effect, df_res] = perm_crANOVA(data, [16, 16], 1e4)
%
%REQUIRED INPUTS
% data          - An electrode x time points x conditions x subjects array of ERP
%                 data. Array will vary in number of dimensions based on how many
%                 factors there are
% cond_subs     - Array giving the number of subjects in each condition of
%                 the between subjects factor. For example, if cond_subs is
%                 [8, 9], then there should be 17 subjects with first 8
%                 being in condition A and the next 9 being in condition B
% n_perm        - Number of permutations to conduct
%
%OUTPUT
% F_obs         - electrode x time point matrix of unpermuted F-values
% F_dist        - F-values at each time point and electrode for each
%                 permutation. The first permutation is F-observed.
% df_effect     - numerator degrees of freedom
% df_res        - denominator degrees of freedom
%
%
%VERSION DATE: 4 April 2019
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.


function [F_obs, F_dist, df_effect, df_res] = perm_crANOVA(data, cond_subs, n_perm)

    global VERBLEVEL

    %Make sure there's only one factor
    assert(ndims(data) == 3);

    %Some useful numbers
    [n_electrodes, n_time_pts, n_subs] = size(data);
    n_conds = length(cond_subs);
    if sum(cond_subs) ~= n_subs
        error('The number of subjects in the ''cond_subs'' input doesn''t match the number of subjects in the data');
    end

    %Calculate degrees of freedom
    %(Always the same, so no point calculating in the loop)
    dfA   = n_conds - 1;
    dfRES = n_subs - n_conds;

    %Perform n_perm permutations
    F_dist = NaN(n_perm, n_electrodes, n_time_pts);

    for i = 1:n_perm

        %Permute the data
        if i ==1
            perm_data = data;
        else
            perm_data = data(:, :, randperm(size(data, 3)));
        end

        %Calculate sums of squares
        A = 0;
        for a = 1:n_conds
            first = sum(cond_subs(1:a)) - cond_subs(a) + 1;
            last  = sum(cond_subs(1:a));
            A = A + ((sum(perm_data(:, :, first:last), 3).^2) / cond_subs(a));
        end
        SSyint = (sum(perm_data, 3).^2) / n_subs;
        SSTO   = sum(perm_data.^2, 3) - SSyint;
        SSA    = A - SSyint;
        SSRES  = SSTO - SSA;

        %Calculate F
        SSA(SSA < 1e-12) = 0; %Eliminates large F values that result from floating point error 
        F_dist(i, :, :) = (SSA/dfA) ./ (SSRES/dfRES);

        if VERBLEVEL
            if i == 1 && n_perm > 1
                fprintf('Permutations completed: ')
            elseif i == n_perm && n_perm > 1
                fprintf('%d\n', i)
            elseif ~mod(i, 1000)
                fprintf('%d, ', i)
            end
        end

    end
    
    %Extract unpermuted F-values
    F_obs = reshape(F_dist(1, :, :), [n_electrodes, n_time_pts]);

    %degrees of freedom
    df_effect = dfA;
    df_res    = dfRES;

end
