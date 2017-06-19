%Calculate parametric ANOVA at each time point and electrode. Results are
%not multiple comparison corrected but can be subjected to FDR correction.
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
% alpha         - Alpha level of the test
%
%OUTPUT
% test_results - A struct with results of the mass univariate ANOVA
%
%
%VERSION DATE: 9 April 2017
%AUTHOR: Eric Fields, Tufts University (Eric.Fields@tufts.edu)
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 
%This is a beta version of this software. It needs additional testing 
%and SHOULD NOT be considered error free.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.
%This function may incorporate some code from the Mass Univariate Toolbox, 
%Copyright (c) 2015, David Groppe: https://github.com/dmgroppe/Mass_Univariate_ERP_Toolbox/blob/master/LICENSE

%%%%%%%%%%%%%%%%%%%  REVISION LOG   %%%%%%%%%%%%%%%%%%%
% 4/9/17  - Added documentation


function test_results = calc_param_ANOVA(data, dims, alpha)
%Reduce data as a necessary and send to the right function for statistical
%test
% 1. Average across any factors not involved in the effect being calculated
% 2. Send the reduced data to the appropriate ANOVA function

    %Average across factors not involved in this effect
    if length(dims) < ndims(data) - 3
        %Put the factors to average across as the initial dimensions
        num_dims = ndims(data);
        all_dims = 1:ndims(data);
        reorder = [all_dims(~ismember(all_dims, [1,2,dims,num_dims])), all_dims(ismember(all_dims, [1,2,dims,num_dims]))];
        reduced_data = permute(data, reorder);
        %Reduce all the factors to average across to a single dimension
        dim_sizes = size(reduced_data);
        num_dims_to_avg = ndims(data) - length(dims) - 3;
        reduced_data = reshape(reduced_data, [prod(dim_sizes(1:num_dims_to_avg)), dim_sizes((num_dims_to_avg+1):end)]);
        %Take the mean across that dimension
        reduced_data = mean(reduced_data, 1);
        %Get rid of the extra singleton dimension
        reduced_data = reshape(reduced_data, dim_sizes((num_dims_to_avg+1):end));
    else
        reduced_data = data;
    end
    
    %Calculate appropriate ANOVA
    if length(dims) == 1
        test_results = oneway(reduced_data, alpha);
    elseif length(dims) == 2
        test_results = twoway_int(reduced_data, alpha);
    elseif length(dims) == 3
        test_results = threeway_int(reduced_data, alpha);
    end

end


function test_results = oneway(data, alpha)

    %Make sure there's only one factor
    assert(ndims(data) == 4);
    
    %Some useful numbers
    [~, ~, n_conds, n_subs] = size(data);
    
    %Calculate degrees of freedom
    %(Always the same, so no point calculating in the loop)
    dfA   = n_conds - 1;
    dfBL  = n_subs - 1;
    dfRES = dfA * dfBL;
        
    %Calculate sums of squares
    SSyint = (sum(sum(data, 3), 4).^2) / (n_conds * n_subs);
   %SSTO   = sum(sum(data.^2, 3), 4) - SSyint;
    SSA    = (sum(sum(data, 4).^2, 3) / n_subs) - SSyint;
    SSBL   = (sum(sum(data, 3).^2, 4) / n_conds) - SSyint;
    SSRES  = sum(sum(data.^2, 3), 4) - SSA - SSBL - SSyint;
    %assert(all(abs(SSTO(:) - (SSA(:) + SSBL(:) + SSRES(:))) < 1e-9));

    %Calculate F
    SSA(SSA < 1e-12) = 0; %Eliminates large F values that result from floating point error 
    F_obs = (SSA/dfA) ./ (SSRES/dfRES);
        
    %Calculate F_crit
    F_crit = finv(.95, dfA, dfRES);

    %Null hypothesis test
    h = F_obs > F_crit;
    
    %Calculate p-value
    p = 1 - fcdf(F_obs, dfA, dfRES);
    assert(isequal(h, p<=alpha));
    
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs;
    test_results.F_crit = F_crit;
    test_results.df = [dfA, dfRES];

end


function test_results = twoway_int(data, alpha)

    %Make sure we're dealing with a two-way design
    assert(ndims(data) == 5);

    %Some useful numbers
    [~, ~, n_conds_A, n_conds_B, n_subs] = size(data);
    
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

    %Calculate sums of squares
    SSyint   = (sum(sum(sum(data, 3), 4), 5).^2)/(n_conds_A*n_conds_B*n_subs);
    %SSTO     = sum(sum(sum(data.^2, 3), 4), 5) - SSyint;
    SSA      = sum(sum(sum(data, 4), 5).^2, 3)/(n_conds_B*n_subs) - SSyint;
    SSB      = sum(sum(sum(data, 3), 5).^2, 4)/(n_conds_A*n_subs) - SSyint;
    SSBL     = sum(sum(sum(data, 3), 4).^2, 5)/(n_conds_A*n_conds_B) - SSyint;
    SSAxB    = sum(sum(sum(data, 5).^2, 3), 4)/n_subs - SSA - SSB -SSyint;
    SSAxBL   = sum(sum(sum(data, 4).^2, 3), 5)/n_conds_B - SSA - SSBL - SSyint;
    SSBxBL   = sum(sum(sum(data, 3).^2, 4), 5)/n_conds_A - SSB - SSBL - SSyint;
    SSAxBxBL = sum(sum(sum(data.^2, 3), 4), 5) - SSA - SSB - SSBL - SSAxB - SSAxBL - SSBxBL - SSyint;
    %SSRES    = sum(sum(sum(data.^2, 3), 4), 5) - SSA - SSB - SSBL - SSAxB - SSyint;

    %Doublechecking that the numbers match up
    %assert(all(SSRES - (SSAxBL + SSBxBL + SSAxBxBL) < 1e-9)); %SSRES is equal to its three subcomponents
    %assert(all(SSTO  - (SSRES + SSBL + SSA + SSB + SSAxB) < 1e-9)); %sums of squares add up

    %Calculate F
    SSAxB(SSAxB < 1e-12) = 0; %Eliminates large F values that result from floating point error 
    F_obs = (SSAxB/dfAxB) ./ (SSAxBxBL/dfAxBerr);
 
    %Calculate F_crit
    F_crit = finv(.95, dfAxB, dfAxBerr);

    %Null hypothesis test
    h = F_obs > F_crit;
    
    %Calculate p-value
    p = 1 - fcdf(F_obs, dfAxB, dfAxBerr);
    assert(isequal(h, p<=alpha));
    
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs;
    test_results.F_crit = F_crit;
    test_results.df = [dfAxB, dfAxBerr];
    
end


function test_results = threeway_int(data, alpha)
    
    %Make sure we're dealing with a two-way design
    assert(ndims(data) == 6);

    %Some useful numbers
    [~, ~, n_conds_A, n_conds_B, n_conds_C, n_subs] = size(data);
    
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

    %Calculate sums of squares
    SSyint     = (sum(sum(sum(sum(data, 3), 4), 5), 6).^2)/(n_conds_A*n_conds_B*n_conds_C*n_subs);
   %SSTO       = sum(sum(sum(sum(data.^2, 3), 4), 5), 6) - SSyint;
    SSA        = sum(sum(sum(sum(data, 4), 5), 6).^2, 3)/(n_conds_B*n_conds_C*n_subs) - SSyint;
    SSB        = sum(sum(sum(sum(data, 3), 5), 6).^2, 4)/(n_conds_A*n_conds_C*n_subs) - SSyint;
    SSC        = sum(sum(sum(sum(data, 3), 4), 6).^2, 5)/(n_conds_A*n_conds_B*n_subs) - SSyint;
    SSBL       = sum(sum(sum(sum(data, 3), 4), 5).^2, 6)/(n_conds_A*n_conds_B*n_conds_C) - SSyint;
    SSAxB      = sum(sum(sum(sum(data, 5), 6).^2, 3), 4)/(n_conds_C*n_subs) - SSA - SSB -SSyint;
    SSAxC      = sum(sum(sum(sum(data, 4), 6).^2, 3), 5)/(n_conds_B*n_subs) - SSA - SSC -SSyint;
    SSBxC      = sum(sum(sum(sum(data, 3), 6).^2, 4), 5)/(n_conds_A*n_subs) - SSB - SSC -SSyint;
    SSAxBxC    = sum(sum(sum(sum(data, 6).^2, 3), 4), 5)/n_subs - SSA - SSB - SSC - SSAxB - SSAxC - SSBxC - SSyint;
    SSAxBL     = sum(sum(sum(sum(data, 4), 5).^2, 3), 6)/(n_conds_B*n_conds_C) - SSA - SSBL - SSyint;
    SSBxBL     = sum(sum(sum(sum(data, 3), 5).^2, 4), 6)/(n_conds_A*n_conds_C) - SSB - SSBL - SSyint;
    SSCxBL     = sum(sum(sum(sum(data, 3), 4).^2, 5), 6)/(n_conds_A*n_conds_B) - SSC - SSBL - SSyint;
    SSAxBxBL   = sum(sum(sum(sum(data, 5).^2, 3), 4), 6)/n_conds_C - SSA - SSB - SSBL - SSAxB - SSAxBL - SSBxBL - SSyint;
    SSAxCxBL   = sum(sum(sum(sum(data, 4).^2, 3), 5), 6)/n_conds_B - SSA - SSC - SSBL - SSAxC - SSAxBL - SSCxBL - SSyint;
    SSBxCxBL   = sum(sum(sum(sum(data, 3).^2, 4), 5), 6)/n_conds_A - SSB - SSC - SSBL - SSBxC - SSBxBL - SSCxBL - SSyint;
    SSAxBxCxBL = sum(sum(sum(sum(data.^2, 3), 4), 5), 6) - SSA - SSB - SSC - SSBL - SSAxB - SSAxC - SSBxC - SSAxBL - SSBxBL - SSCxBL - SSAxBxC - SSAxBxBL - SSAxCxBL - SSBxCxBL - SSyint;
   %SSRES      = sum(sum(sum(sum(data.^2, 3), 4), 5), 6) - SSA - SSB - SSC - SSBL - SSAxB - SSAxC - SSBxC - SSAxBxC - SSyint;

    %Doublechecking that the numbers match up
    %assert(all(SSRES - (SSAxBL + SSBxBL + SSCxBL + SSAxBxBL + SSAxCxBL + SSBxCxBL + SSAxBxCxBL) < 1e-9)); %SSRES is equal to its three subcomponents
    %assert(all(SSTO  - (SSRES + SSBL + SSA + SSB + SSC + SSAxB + SSAxC + SSBxC + SSAxBxC) < 1e-9)); %sums of squares add up

    SSAxBxC(SSAxBxC < 1e-12) = 0; %Eliminates large F values that result from floating point error 
    F_obs = (SSAxBxC/dfAxBxC) ./ (SSAxBxCxBL/dfAxBxCerr);
    
    %Calculate Fobs, Fmax distribution, and Fmax critical value
    F_crit = finv(.95, dfAxBxC, dfAxBxCerr);

    %Null hypothesis test
    h = F_obs > F_crit;
    
    %Calculate p-value
    p = 1 - fcdf(F_obs, dfAxBxC, dfAxBxCerr);
    assert(isequal(h, p<=alpha));
    
    test_results.h = h;
    test_results.p = p;
    test_results.F_obs = F_obs;
    test_results.F_crit = F_crit;
    test_results.df = [dfAxBxC, dfAxBxCerr];

end
