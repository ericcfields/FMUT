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
%
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
% 7/13/17   - Now reduces all interactions effects to the maximum extent
%             possible

function reduced_data = reduce_data(data, dims)

    %% Average across factors not involved in this effect
    
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
    
    %% For interactions, subtract across factors with two levels
    
    dim_sizes = size(reduced_data);
    factor_levels = dim_sizes(3:(ndims(reduced_data)-1));
    if length(dims)>1 && sum(factor_levels==2)
        %Put the factors to subtract across as the initial dimensions
        if sum(factor_levels>2)
            reorder = [find(factor_levels==2)+2, 1, 2, find(factor_levels~=2)+2, ndims(reduced_data)];
        else
            reorder = [4:(ndims(reduced_data)-1), 1, 2, 3, ndims(reduced_data)]; 
        end
        reduced_data = permute(reduced_data, reorder);
        %Flatten data and subtract until the data is reduced to the right size
        reduced_data = reshape(reduced_data, 1, []);
        if sum(factor_levels>2)
            end_size = size(data,1) * size(data,2) * prod(factor_levels(factor_levels~=2)) * size(data,ndims(data));
        else
            end_size = size(data,1) * size(data,2) * 2 * size(data,ndims(data));
        end
        while length(reduced_data) > end_size
            reduced_data = reduced_data(1:2:length(reduced_data)) - reduced_data(2:2:length(reduced_data));
        end
        %Put back in regular format
        if sum(factor_levels>2)
            reduced_data = reshape(reduced_data, [size(data,1), size(data,2), factor_levels(factor_levels~=2), size(data,ndims(data))]);
        else
            reduced_data = reshape(reduced_data, [size(data, 1), size(data, 2), 2, size(data, ndims(data))]);
        end
    end
    
end

