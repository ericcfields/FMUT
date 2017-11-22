%Calculate F-observed and the empirical F-distribution for a
%within-subjects ANOVA with up to three factors. This function calculates a
%one-way ANOVA, two-way interaction, and three-way interaction.
%
%EXAMPLE USAGE
% >> [F_dist, df_effect, df_res, exact_test] = perm_rbANOVA(data, [3, 4], 1e4)
%
%REQUIRED INPUTS
% data          - An electrode x time points x conditions x subjects array of ERP
%                 data. Array will vary in number of dimensions based on how many
%                 factors there are.
% n_perm        - Number of permutations to conduct
%
%OUTPUT
% Fvals         - F-values at each time point and electrode for each
%                 permutation. The first permutation is F-observed.
% df_effect     - numerator degrees of freedom
% df_res        - denominator degrees of freedom
% exact_test    - Boolean specifying whether the test was an exact test
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

function [F_dist, df_effect, df_res, exact_test] = perm_rbANOVA(data, dims, n_perm, reduce)

    %Eliminate factors not involved in this effect and reduce interactions
    %via subtraction
    if nargin < 4
        reduce = true;
    end
    if reduce
        reduced_data = reduce_data(data, dims);
    else
        reduced_data = data;
    end

    %Calculate appropriate ANOVA
    if ndims(reduced_data) == 4
        [F_dist, df_effect, df_res] = oneway(reduced_data, n_perm);
        exact_test = true;
    elseif ndims(reduced_data) == 5
        [F_dist, df_effect, df_res] = twoway_approx_int(reduced_data, n_perm);
        exact_test = false;
    elseif ndims(reduced_data) == 6
        [F_dist, df_effect, df_res] = threeway_approx_int(reduced_data, n_perm);
        exact_test = false;
    end

end

function [F_dist, df_effect, df_res] = oneway(data, n_perm)
%Perform permutation one-way ANOVA. This is an exact test.

    global VERBLEVEL

    %Make sure there's only one factor
    assert(ndims(data) == 4);
    
    %Some useful numbers
    [n_electrodes, n_time_pts, n_conds, n_subs] = size(data);
    
    %Calculate degrees of freedom
    %(Always the same, so no point calculating in the loop)
    dfA   = n_conds - 1;
    dfBL  = n_subs - 1;
    dfRES = dfA * dfBL;
    
    %Perform n_perm permutations
    F_dist = NaN(n_perm, n_electrodes, n_time_pts);
    for i = 1:n_perm
        
        %Permute the data
        if i ==1
            perm_data = data;
        else
            for n = 1:n_subs
                perm_data(:, :, :, n) = data(:, :, randperm(size(data, 3)), n);
            end
        end
        
        %Calculate sums of squares
        SSyint = (sum(sum(perm_data, 3), 4).^2) / (n_conds * n_subs);
       %SSTO   = sum(sum(perm_data.^2, 3), 4) - SSyint;
        SSA    = (sum(sum(perm_data, 4).^2, 3) / n_subs) - SSyint;
        SSBL   = (sum(sum(perm_data, 3).^2, 4) / n_conds) - SSyint;
        SSRES  = sum(sum(perm_data.^2, 3), 4) - SSA - SSBL - SSyint;
        %assert(all(abs(SSTO(:) - (SSA(:) + SSBL(:) + SSRES(:))) < 1e-9));
        
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


function [F_dist, df_effect, df_res] = twoway_approx_int(data, n_perm)
%Use permutation of residuals method to conduct an approximate test of the
%two-way interaction.

    global VERBLEVEL

    %Make sure we're dealing with a two-way design
    assert(ndims(data) == 5);

    %Some useful numbers
    [n_electrodes, n_time_pts, n_conds_A, n_conds_B, n_subs] = size(data);

    %Subtract main effects within each subject so that the data is 
    %exchangeable under the null hypothesis for the interaction
    int_res = get_int_res(data, [], [3, 4]);
    
    %Calculate degrees of freedom
    %(Always the same, so no point calculating in the loop)
    dfBL     = n_subs - 1;
    dfA      = n_conds_A - 1;
    %dfAerr   = dfA * dfBL;
    dfB      = n_conds_B - 1;
    %dfBerr   = dfB * dfBL;
    dfAxB    = dfA * dfB;
    dfAxBerr = dfAxB * dfBL;
    %dfRES    = (num_subs - 1) * (num_conds_A * num_conds_B - 1);

    %Re-arrange data for permutation
    flat_data = reshape(int_res, [n_electrodes, n_time_pts, n_conds_A*n_conds_B, n_subs]);

    %Perform n_perm permutations
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
            perm_data = reshape(flat_perm_data, n_electrodes, n_time_pts, n_conds_A, n_conds_B, n_subs);
        end
        
        %Calculate sums of squares
        SSyint   = (sum(sum(sum(perm_data, 3), 4), 5).^2)/(n_conds_A*n_conds_B*n_subs);
        %SSTO     = sum(sum(sum(perm_data.^2, 3), 4), 5) - SSyint;
        SSA      = sum(sum(sum(perm_data, 4), 5).^2, 3)/(n_conds_B*n_subs) - SSyint;
        SSB      = sum(sum(sum(perm_data, 3), 5).^2, 4)/(n_conds_A*n_subs) - SSyint;
        SSBL     = sum(sum(sum(perm_data, 3), 4).^2, 5)/(n_conds_A*n_conds_B) - SSyint;
        SSAxB    = sum(sum(sum(perm_data, 5).^2, 3), 4)/n_subs - SSA - SSB -SSyint;
        SSAxBL   = sum(sum(sum(perm_data, 4).^2, 3), 5)/n_conds_B - SSA - SSBL - SSyint;
        SSBxBL   = sum(sum(sum(perm_data, 3).^2, 4), 5)/n_conds_A - SSB - SSBL - SSyint;
        SSAxBxBL = sum(sum(sum(perm_data.^2, 3), 4), 5) - SSA - SSB - SSBL - SSAxB - SSAxBL - SSBxBL - SSyint;
        %SSRES    = sum(sum(sum(perm_data.^2, 3), 4), 5) - SSA - SSB - SSBL - SSAxB - SSyint;

        %Doublechecking that the numbers match up
        %assert(all(SSRES - (SSAxBL + SSBxBL + SSAxBxBL) < 1e-9)); %SSRES is equal to its three subcomponents
        %assert(all(SSTO  - (SSRES + SSBL + SSA + SSB + SSAxB) < 1e-9)); %sums of squares add up

        %Calculate F
        SSAxB(SSAxB < 1e-12) = 0; %Eliminates large F values that result from floating point error 
        F_dist(i, :, :) = (SSAxB/dfAxB) ./ (SSAxBxBL/dfAxBerr);
        
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
    df_effect = dfAxB;
    df_res    = dfAxBerr;

end


function [F_dist, df_effect, df_res] = threeway_approx_int(data, n_perm)
%Use permutation of residuals method to conduct an approximate test of the
%three-way interaction.
    
    global VERBLEVEL

    %Make sure we're dealing with a two-way design
    assert(ndims(data) == 6);

    %Some useful numbers
    [n_electrodes, n_time_pts, n_conds_A, n_conds_B, n_conds_C, n_subs] = size(data);

    %Subtract main effects then two-way effects within each subject so that 
    %the data is exchangeable under the null hypothesis for the three-way
    %interaction
    int_res = get_int_res(data, [], [3, 4, 5]);
    
    %Calculate degrees of freedom
    %(Always the same, so no point calculating in the loop)
    dfBL       = n_subs - 1;
    dfA        = n_conds_A - 1;
   %dfAerr     = dfA * dfBL;
    dfB        = n_conds_B - 1;
   %dfBerr     = dfB * dfBL;
    dfC        = n_conds_C - 1;
   %dfCerr     = dfC * dfBL;
   %dfAxB      = dfA * dfB;
   %dfAxBerr   = dfAxB * dfBL;
   %dfAxC      = dfA * dfC;
   %dfAxCerr   = dfAxC * dfBL;
   %dfBxC      = dfB * dfC;
   %dfBxCerr   = dfBxC * dfBL;
    dfAxBxC    = dfA * dfB * dfC;
    dfAxBxCerr = dfAxBxC * dfBL;
   %dfRES      = (n_subs - 1) * (n_conds_A * n_conds_B * n_conds_C - 1);

    %Re-arrange data for permutation
    flat_data = reshape(int_res, [n_electrodes, n_time_pts, n_conds_A*n_conds_B*n_conds_C, n_subs]);

    %Perform n_perm permutations
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
            perm_data = reshape(flat_perm_data, size(int_res));
        end
        
        %Calculate F at each time point and electrode combination

        %Calculate sums of squares
        SSyint     = (sum(sum(sum(sum(perm_data, 3), 4), 5), 6).^2)/(n_conds_A*n_conds_B*n_conds_C*n_subs);
       %SSTO       = sum(sum(sum(sum(perm_data.^2, 3), 4), 5), 6) - SSyint;
        SSA        = sum(sum(sum(sum(perm_data, 4), 5), 6).^2, 3)/(n_conds_B*n_conds_C*n_subs) - SSyint;
        SSB        = sum(sum(sum(sum(perm_data, 3), 5), 6).^2, 4)/(n_conds_A*n_conds_C*n_subs) - SSyint;
        SSC        = sum(sum(sum(sum(perm_data, 3), 4), 6).^2, 5)/(n_conds_A*n_conds_B*n_subs) - SSyint;
        SSBL       = sum(sum(sum(sum(perm_data, 3), 4), 5).^2, 6)/(n_conds_A*n_conds_B*n_conds_C) - SSyint;
        SSAxB      = sum(sum(sum(sum(perm_data, 5), 6).^2, 3), 4)/(n_conds_C*n_subs) - SSA - SSB -SSyint;
        SSAxC      = sum(sum(sum(sum(perm_data, 4), 6).^2, 3), 5)/(n_conds_B*n_subs) - SSA - SSC -SSyint;
        SSBxC      = sum(sum(sum(sum(perm_data, 3), 6).^2, 4), 5)/(n_conds_A*n_subs) - SSB - SSC -SSyint;
        SSAxBxC    = sum(sum(sum(sum(perm_data, 6).^2, 3), 4), 5)/n_subs - SSA - SSB - SSC - SSAxB - SSAxC - SSBxC - SSyint;
        SSAxBL     = sum(sum(sum(sum(perm_data, 4), 5).^2, 3), 6)/(n_conds_B*n_conds_C) - SSA - SSBL - SSyint;
        SSBxBL     = sum(sum(sum(sum(perm_data, 3), 5).^2, 4), 6)/(n_conds_A*n_conds_C) - SSB - SSBL - SSyint;
        SSCxBL     = sum(sum(sum(sum(perm_data, 3), 4).^2, 5), 6)/(n_conds_A*n_conds_B) - SSC - SSBL - SSyint;
        SSAxBxBL   = sum(sum(sum(sum(perm_data, 5).^2, 3), 4), 6)/n_conds_C - SSA - SSB - SSBL - SSAxB - SSAxBL - SSBxBL - SSyint;
        SSAxCxBL   = sum(sum(sum(sum(perm_data, 4).^2, 3), 5), 6)/n_conds_B - SSA - SSC - SSBL - SSAxC - SSAxBL - SSCxBL - SSyint;
        SSBxCxBL   = sum(sum(sum(sum(perm_data, 3).^2, 4), 5), 6)/n_conds_A - SSB - SSC - SSBL - SSBxC - SSBxBL - SSCxBL - SSyint;
        SSAxBxCxBL = sum(sum(sum(sum(perm_data.^2, 3), 4), 5), 6) - SSA - SSB - SSC - SSBL - SSAxB - SSAxC - SSBxC - SSAxBL - SSBxBL - SSCxBL - SSAxBxC - SSAxBxBL - SSAxCxBL - SSBxCxBL - SSyint;
       %SSRES      = sum(sum(sum(sum(perm_data.^2, 3), 4), 5), 6) - SSA - SSB - SSC - SSBL - SSAxB - SSAxC - SSBxC - SSAxBxC - SSyint;

        %Doublechecking that the numbers match up
        %assert(all(SSRES - (SSAxBL + SSBxBL + SSCxBL + SSAxBxBL + SSAxCxBL + SSBxCxBL + SSAxBxCxBL) < 1e-9)); %SSRES is equal to its three subcomponents
        %assert(all(SSTO  - (SSRES + SSBL + SSA + SSB + SSC + SSAxB + SSAxC + SSBxC + SSAxBxC) < 1e-9)); %sums of squares add up

        SSAxBxC(SSAxBxC < 1e-12) = 0; %Eliminates large F values that result from floating point error 
        F_dist(i, :, :) = (SSAxBxC/dfAxBxC) ./ (SSAxBxCxBL/dfAxBxCerr);
        
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
    df_effect = dfAxBxC;
    df_res    = dfAxBxCerr;

end