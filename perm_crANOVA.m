%Permutation one-way completely randomized (between subjects) ANOVA
%
%VERSION DATE: 3 July 2017
%AUTHOR: Eric Fields

function [F_dist, df_effect, df_res] = perm_crANOVA(data, cond_subs, n_perm)

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

    for i = 1:n_perm;

        %Permute the data
        if i ==1
            perm_data = data;
        else
            perm_data = data(:, :, randperm(size(data, 3)));
        end

        %Calculate sums of squares
        A = 0;
        for a = 1:n_conds;
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

    %degrees of freedom
    df_effect = dfA;
    df_res    = dfRES;

end