%Calculate F-observed and the empirical F-distribution for an ANOVA with
%one between subjects factor and up to two within subjects factors
%
%REQUIRED INPUTS
% data          - An electrode x time points x conditions x subjects array of ERP
%                 data. Array will vary in number of dimensions based on how many
%                 factors there are
% cond_subs     - Array giving the number of subjects in each condition of
%                 the between subjects factor. For example, if cond_subs is
%                 [8, 9], then there should be 17 subjects with first 8
%                 being in condition A and the next 9 being in condition B
% dims          - Dimensions involved in the effect. Given the structure of
%                 the data specified above, for a two-way ANOVA 3 indicates 
%                 the within-subjects factors, 4 indicates the between 
%                 subjects factors, and [3, 4] indicates the interaction
% n_perm        - Number of permutations to conduct
%
%OUTPUT
% F_dist        - F-values at each time point and electrode for each
%                 permutation. The first permutation is F-observed.
% df_effect     - numerator degrees of freedom
% df_res        - denominator degrees of freedom
%
%
%VERSION DATE: 5 July 2017
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 
%This is a beta version of this software. It needs additional testing 
%and SHOULD NOT be considered error free.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function [F_dist, df_effect, df_res] = perm_spANOVA(data, cond_subs, dims, n_perm)

    if ndims(data) == 3 && dims == 3
        [F_dist, df_effect, df_res] = perm_crANOVA(data, cond_subs, n_perm);
    elseif ndims(data) == 4
        [F_dist, df_effect, df_res] = twoway(data, cond_subs, dims, n_perm);
    elseif ndims(data) == 5
        watchtit('Three-way split plot ANOVA is not implemented yet.');
    end
    
end

function [F_dist, df_effect, df_res] = twoway(data, cond_subs, dims, n_perm)

    global VERBLEVEL

    %Check array structure
    assert(ndims(data) == 4);

    %Some useful numbers
    [n_electrodes, n_time_pts, n_conds_B, n_subs] = size(data);
    n_conds_A = length(cond_subs);
    if sum(cond_subs) ~= n_subs
        error('The number of subjects in the ''cond_subs'' input doesn''t match the number of subjects in the data');
    end
    
    %Interaction residuals
    if length(dims) == 2
        int_res = get_int_res_sp2(data, cond_subs);
    end

    %Calculate degrees of freedom
    %(Always the same, so no point calculating in the loop)
    dfA    = n_conds_A - 1;
    dfB    = n_conds_B - 1;
    dfBL   = n_subs - n_conds_A;
    dfAxB  = dfA * dfB;
    dfBxBL = dfB * dfBL;

    %Perform n_perm permutations
    F_dist = NaN(n_perm, n_electrodes, n_time_pts);
    for i = 1:n_perm;

        %Permute the data
        if length(dims) == 1
            if i == 1
                perm_data = data;
            elseif dims == 3
                for n = 1:n_subs
                    perm_data(:, :, :, n) = data(:, :, randperm(size(data, 3)), n);
                end
            elseif dims == 4 
                perm_data = data(:, :, :, randperm(size(data, 4)));
            end
        elseif length(dims) ==2
            if i ==1
                perm_data = int_res;
            else
                for n = 1:n_subs
                    perm_data(:, :, :, n) = int_res(:, :, randperm(size(data, 3)), n);
                end
                perm_data = perm_data(:, :, :, randperm(size(data, 4)));
            end
        end

        %Calculate sums of squares
        A = 0; AS = 0; AB = 0; ABS = 0;
        for p = 1:n_conds_A;
            first = sum(cond_subs(1:p)) - cond_subs(p) + 1;
            last  = sum(cond_subs(1:p));
            A   = A   + sum(sum(perm_data(:, :, :, first:last), 3), 4).^2 / (cond_subs(p) * n_conds_B);
            AS  = AS  + sum(sum(perm_data(:, :, :, first:last), 3).^2, 4) / n_conds_B;
            AB  = AB  + sum(sum(perm_data(:, :, :, first:last), 4).^2, 3) / cond_subs(p);
            ABS = ABS + sum(sum(perm_data(:, :, :, first:last).^2, 3), 4);
        end
        SSyint = (sum(sum(perm_data, 3), 4).^2) / (n_subs * n_conds_B);
        %SSTO   = sum(sum(perm_data.^2, 3), 4) - SSyint;
        SSA    = A - SSyint;
        SSB    = sum(sum(perm_data, 4).^2, 3) / n_subs - SSyint;
        SSBL   = AS - A;
        SSAxB  = AB - A - SSB; 
        SSBxBL = ABS - AB - AS + A;
        
        %Calculate F
        if length(dims) == 1
            if dims == 3
                SSAxB(SSAxB < 1e-12) = 0; %Eliminates large F values that result from floating point error 
                F_dist(i, :, :) = (SSAxB/dfAxB) ./ (SSBxBL/dfBxBL);
            elseif dims == 4
                SSA(SSA < 1e-12) = 0; %Eliminates large F values that result from floating point error 
                F_dist(i, :, :) = (SSA/dfA) ./ (SSBL/dfBL);
            end
        elseif length(dims) == 2
            SSAxB(SSAxB < 1e-12) = 0; %Eliminates large F values that result from floating point error 
            F_dist(i, :, :) = (SSAxB/dfAxB) ./ (SSBxBL/dfBxBL);
        end
        
        %Report permutations completed to command window
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
    if length(dims) == 1
        if dims == 3
            df_effect = dfB;
            df_res    = dfBxBL;
        elseif dims == 4
            df_effect = dfA;
            df_res    = dfBL;
        end
    elseif length(dims) == 2
        df_effect = dfAxB;
        df_res    = dfBxBL;
    end
    
end

function int_res = get_int_res_sp2(data, cond_subs)

    n_conds_B = size(data, 3);
    n_conds_A = length(cond_subs);

    int_res = NaN(size(data));
    for p = 1:n_conds_A
        first = sum(cond_subs(1:p)) - cond_subs(p) + 1;
        for q = 1:n_conds_B
            for n = 1:cond_subs(p)
                sub = first+n-1;
                int_res(:,:,q,sub) = data(:,:,q,sub) ...
                                   - mean(data(:,:,q, :), 4) ...
                                   - mean(data(:, :, :, sub), 3) ...
                                   + mean(mean(data, 3), 4);
            end
        end
    end
    
    assert(abs(mean(int_res(:))) < 1e-9);

end
