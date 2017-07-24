%Reduce data to the simplest design that can calculate an equivalent
%effect. For example, for a 3 x 2 design, the main effects can be
%calculated by averaging across the other factor and then conducting a
%one-way ANOVA. Similarly, the interaction effect an be calculated by
%subtracting across the factor with two levels and then conducting a
%one-way ANOVA.
%
%reduce_data performs two computations:
%1. Average across any factors not included in the effect to be calculated
%   (as specified by the dims input)
%2. Subtract across any factors involved in an interaction that have only 
%   two levels (unless all factors have two levels, in which case subtract
%   across all but one factor)
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
%
%OUTPUT
% reduced_data  - data reduced for analysis
% new_dims      - the dimensions of reduced_data that are involved in the
%                 effect
%
%VERSION DATE: 13 July 2017
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

%%%%%%%%%%%%%%%%%%%  REVISION LOG   %%%%%%%%%%%%%%%%%%%
% 6/22/17   - First version. Code re-organized from other functions.
% 7/10/17   - Ignore between-subjects factor; return updated dims variable
% 7/13/17   - Now reduces all interactions effects to the maximum extent
%             possible

function [reduced_data, new_dims] = reduce_data(data, dims)

    %Can't reduce across betwee-subjects factors, so extract just the
    %within-subject factors
    wdims = dims(dims ~= ndims(data));
    
    n_electrodes = size(data, 1);
    n_time_pts   = size(data, 2);
    n_subs = size(data, ndims(data));
    
    %Determine if betwee-subjects factors are involved in the effect
    bg = any(dims == ndims(data));
    
    
    %% Average across factors not involved in this effect
    
    if (length(wdims) < ndims(data) - 3) || isempty(wdims)
        %Put the factors to average across as the initial dimensions
        num_dims = ndims(data);
        all_dims = 1:ndims(data);
        reorder = [all_dims(~ismember(all_dims, [1,2,wdims,num_dims])), all_dims(ismember(all_dims, [1,2,wdims,num_dims]))];
        reduced_data = permute(data, reorder);
        %Reduce all the factors to average across to a single dimension
        dim_sizes = size(reduced_data);
        num_dims_to_avg = ndims(data) - length(wdims) - 3;
        reduced_data = reshape(reduced_data, [prod(dim_sizes(1:num_dims_to_avg)), dim_sizes((num_dims_to_avg+1):end)]);
        %Take the mean across that dimension
        reduced_data = mean(reduced_data, 1);
        %Get rid of the extra singleton dimension
        reduced_data = reshape(reduced_data, dim_sizes((num_dims_to_avg+1):end));
    else
        reduced_data = data;
    end
    
    %% For interactions, subtract across factors with two levels
    
    dim_sizes = size(reduced_data);
    factor_levels = dim_sizes(3:(ndims(reduced_data)-1));
    if length(dims)>1 && sum(factor_levels==2)
        %Put the factors to subtract across as the initial dimensions
        if sum(factor_levels>2)
            reorder = [find(factor_levels==2)+2, 1, 2, find(factor_levels~=2)+2, ndims(reduced_data)];
        elseif bg
            reorder = [3:(ndims(reduced_data)-1), 1, 2, ndims(reduced_data)];     
        else
            reorder = [4:(ndims(reduced_data)-1), 1, 2, 3, ndims(reduced_data)];
        end
        reduced_data = permute(reduced_data, reorder);
        %Flatten data and subtract until the data is reduced to the right size
        reduced_data = reshape(reduced_data, 1, []);
        if sum(factor_levels>2)
            end_size = n_electrodes * n_time_pts * prod(factor_levels(factor_levels~=2)) * n_subs;
        elseif bg
            end_size = n_electrodes * n_time_pts * n_subs;
        else
            end_size = n_electrodes * n_time_pts * 2 * n_subs;
        end
        while length(reduced_data) > end_size
            reduced_data = reduced_data(1:2:length(reduced_data)) - reduced_data(2:2:length(reduced_data));
        end
        %Put back in regular format
        if sum(factor_levels>2)
            reduced_data = reshape(reduced_data, [n_electrodes, n_time_pts, factor_levels(factor_levels~=2), size(data,ndims(data))]);
        elseif bg
            reduced_data = reshape(reduced_data, [n_electrodes, n_time_pts, n_subs]);
        else
            reduced_data = reshape(reduced_data, [n_electrodes, n_time_pts, 2, n_subs]);
        end
    end
    
    %% Dimensions of reduced_data involved in effect
    
    new_dims = (3:ndims(reduced_data)-1);
    if bg
        new_dims = [new_dims ndims(reduced_data)];
    end
    
end

