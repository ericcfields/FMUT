%reduce_data performs two computations:
%1. Average across any factors not included in an effect to be calculated. 
%   For example, in a two-way design, the main effect of factor A can be 
%   calculated by averaging across factor B and calculating a one-way ANOVA.
%2. For exact interactions, reduces to a single factor via subtraction. For
%   example a 2x2 interaction can be calculated by subtracting across one 
%   of the factors and calculating a one-way ANOVA.
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
% int_method    - A string that should be either 'exact' or 'approximate'.
%                 If 'exact', the function will attempt to reduce to a 
%                 one-way design by subtraction 
%
%OUTPUT
% reduced_data  - data reduced for analysis
% new_dims      - the dimensions of reduced_data that are involved in the
%                 effect
%
%
%VERSION DATE: 10 July 2017
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

function [reduced_data, new_dims] = reduce_data(data, dims, int_method)

    %Can't reduce across betwee-subjects factors, so extract just the
    %within-subject factors
    wdims = dims(dims ~= ndims(data));
    
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
    
    
    %% For exact interactions, reduce data to one-way design via subtraction
    
    if length(wdims) > 1 && strcmpi(int_method, 'exact')
        %Get/check design structure
        dim_sizes = size(reduced_data);
        factor_levels = dim_sizes(3:(ndims(reduced_data)-1));
        assert(sum(factor_levels > 2) < 2);
        %Put the factors to subtract across as the initial dimensions
        [~, factor_order] = sort(factor_levels);
        reorder = [factor_order(1:end-1)+2, 1, 2, factor_order(end)+2 ndims(reduced_data)];
        reduced_data = permute(reduced_data, reorder);
        %Flatten data and subtract until the data is reduced to the right size
        reduced_data = reshape(reduced_data, 1, []);
        while length(reduced_data) > size(data,1)*size(data,2)*max(factor_levels)*size(data,ndims(data))
            reduced_data = reduced_data(1:2:length(reduced_data)) - reduced_data(2:2:length(reduced_data));
        end
        %Put back in regular format
        reduced_data = reshape(reduced_data, [size(data,1), size(data,2), max(factor_levels), size(data,ndims(data))]);
    end
    
    %% Dimensions of reduced_data involved in effect
    
    new_dims = (3:ndims(reduced_data)-1);
    if any(dims == ndims(data)) && ndims(reduced_data) > 3
        new_dims = [new_dims ndims(reduced_data)];
    end
    
end

