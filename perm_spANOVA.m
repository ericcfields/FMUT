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
%VERSION DATE: 14 July 2017
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
    elseif ndims(data) > 4
        [F_dist, df_effect, df_res] = threeway(data, cond_subs, dims, n_perm);
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
        int_res = get_int_res(data, cond_subs, dims);
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
                SSB(SSB < 1e-12) = 0; %Eliminates large F values that result from floating point error 
                F_dist(i, :, :) = (SSB/dfB) ./ (SSBxBL/dfBxBL);
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


function [F_dist, df_effect, df_res] = threeway(data, cond_subs, dims, n_perm)

    global VERBLEVEL

    %Check array structure
    assert(ndims(data) == 5);

    %Some useful numbers
    [n_electrodes, n_time_pts, n_conds_B, n_conds_C, n_subs] = size(data);
    n_conds_A = length(cond_subs);
    if sum(cond_subs) ~= n_subs
        error('The number of subjects in the ''cond_subs'' input doesn''t match the number of subjects in the data');
    end
    
    %Interaction residuals
    int_res = get_int_res(data, cond_subs, dims);
    
    %Re-arrange data for within-subjects permutation
    flat_data = reshape(int_res, [n_electrodes, n_time_pts, n_conds_B*n_conds_C, n_subs]);

    %Calculate degrees of freedom
    %(Always the same, so no point calculating in the loop)
    dfA      = n_conds_A - 1;
    dfB      = n_conds_B - 1;
    dfC      = n_conds_C - 1;
    dfBL     = n_subs - n_conds_A;
    %dfAxB    = dfA * dfB;
    %dfAxC    = dfA * dfC;
    dfBxC    = dfB * dfC;
    dfAxBxC  = dfA * dfB * dfC;
    %dfBxBL   = dfB * dfBL;
    %dfCxBL   = dfC * dfBL;
    dfBxCxBL = dfB * dfC * dfBL;
    
    F_dist = NaN(n_perm, n_electrodes, n_time_pts);
    flat_perm_data = NaN(size(flat_data));
    for i = 1:n_perm
        
        %Permute the data
        if i == 1
            perm_data = int_res;
        else
            for s = 1:n_subs
                flat_perm_data(:, :, :, s) = flat_data(:, :, randperm(size(flat_data, 3)), s);
            end
            perm_data = reshape(flat_perm_data, [n_electrodes, n_time_pts, n_conds_B, n_conds_C, n_subs]);
            if any(dims == 5)
                perm_data = perm_data(:, :, :, :, randperm(n_subs));
            end
        end
        
        %Calculate sums of squares
        A = 0; AB = 0; AC = 0; AS = 0; ABS = 0; ACS = 0; ABC = 0; ABCS = 0;
        for p = 1:n_conds_A
            first = sum(cond_subs(1:p)) - cond_subs(p) + 1;
            last  = sum(cond_subs(1:p));
            A    = A    + sum(sum(sum(perm_data(:, :, :, :, first:last), 3), 4), 5).^2 / (cond_subs(p) * n_conds_B * n_conds_C);
            AB   = AB   + sum(sum(sum(perm_data(:, :, :, :, first:last), 4), 5).^2, 3) / (cond_subs(p) * n_conds_C);
            AC   = AC   + sum(sum(sum(perm_data(:, :, :, :, first:last), 3), 5).^2, 4) / (cond_subs(p) * n_conds_B);
            AS   = AS   + sum(sum(sum(perm_data(:, :, :, :, first:last), 3), 4).^2, 5) / (n_conds_B * n_conds_C);
            ABS  = ABS  + sum(sum(sum(perm_data(:, :, :, :, first:last), 4).^2, 3), 5) / n_conds_C;
            ACS  = ACS  + sum(sum(sum(perm_data(:, :, :, :, first:last), 3).^2, 4), 5) / n_conds_B;
            ABC  = ABC  + sum(sum(sum(perm_data(:, :, :, :, first:last), 5).^2, 3), 4) / cond_subs(p);
            ABCS = ABCS + sum(sum(sum(perm_data(:, :, :, :, first:last).^2, 3), 4), 5);
        end
        SSyint   = sum(perm_data(:)).^2 / (n_subs * n_conds_B * n_conds_C);
        SSA      = A - SSyint;
        SSB      = sum(sum(sum(perm_data, 4), 5).^2, 3) / (n_subs * n_conds_C) - SSyint;
        SSC      = sum(sum(sum(perm_data, 3), 5).^2, 4) / (n_subs * n_conds_B) - SSyint;
        %SSBL     = AS - A;
        SSAxB    = AB - SSA - SSB - SSyint;
        SSAxC    = AC - SSA - SSC - SSyint;
        SSBxC    = sum(sum(sum(perm_data, 5).^2, 3), 4) / n_subs - SSB - SSC - SSyint;
        %SSBxBL   = ABS - AB - AS + A;
        %SSCxBL   = ACS - AC - AS + A;
        SSAxBxC  = ABC - SSAxB - SSAxC - SSBxC - SSA - SSB - SSC - SSyint;
        SSBxCxBL = ABCS - ABC - ABS - ACS + AB + AC + AS - A;
        
        %Calculate F
        if isequal(dims, [3, 4])
            SSBxC(SSBxC < 1e-12) = 0; %Eliminates large F values that result from floating point error 
            F_dist(i, :, :) = (SSBxC/dfBxC) ./ (SSBxCxBL/dfBxCxBL);
        elseif isequal(dims, [3, 4, 5])
            SSAxBxC(SSAxBxC < 1e-12) = 0; %Eliminates large F values that result from floating point error 
            F_dist(i, :, :) = (SSAxBxC/dfAxBxC) ./ (SSBxCxBL/dfBxCxBL);
        else
            error('Something has gone wrong! This design should have been reduced.');
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
    if isequal(dims, [3, 4])
        df_effect = dfBxC;
    elseif isequal(dims, [3, 4, 5])
        df_effect = dfAxBxC;
    end
    df_res = dfBxCxBL;

end

