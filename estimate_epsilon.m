%Calculate the Greenhouse-Geisser, Hyunh-Feldt, or lower bound estimate 
%of epsilon at all time points and electrodes
%
%EXAMPLE USAGE
% >> epsilon = estimate_epsilon(data, [], 3, 'gg');
%
%REQUIRED INPUTS
% data          - An electrode x time points x conditions x subjects array 
%                 of ERP data. Array will vary in number of dimensions based 
%                 on how many factors there are.
% cond_subs     - Array giving the number of subjects in each condition of
%                 the between subjects factor. For example, if cond_subs is
%                 [8, 9], then there should be 17 subjects with the first 8
%                 being in condition A and the next 9 being in condition B.
%                 For fully within-subjects designs cond_subs = []
% dims          - Dimensions of the data array involved in the effect to be
%                 calculated. For example, if data is an electrode x time 
%                 points x Factor A x Factor B x subjects array and you want to
%                 calculate the main effect of A, dims = 3. If you want to
%                 calculate the AxB interaciton, dims  = [3, 4].
% method        - 'gg': Greenhouse-Geisser, 'hf':Hyunh-Feldt, or 'lb':
%                 lower bound {default: 'gg'}
%
%OUTPUT
% epsilon      - electrode x time point array of epsilon estimate
%
%
%VERSION DATE: 11 June 2020
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 

%Copyright (c) 2020, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function epsilon = estimate_epsilon(data, cond_subs, dims, method)
    
    %Default to Greenhouse-Geisser method
    if nargin < 4
        method = 'gg';
    end
    
    %Check input
    if ~any(strcmpi(method, {'gg', 'hf', 'lb'}))
        error('method input must be ''gg'', ''hf'' or ''lb''');
    end

    %Get data reduced for analysis
    reduced_data = reduce_data(data, dims);
    
    %This function currently only works for one-way fully repeated measures designs
    if ~isempty(cond_subs) && ~isequal(cond_subs, 0) && length(cond_subs) > 1
        error('estimate_epsilon curently only implemented for fully repetaed measures ANOVA');
    end
    if ndims(reduced_data) ~= 4
        error('estimate_epsilon currently does not support more than one factor with more than 2 levels');
    end
    
    %Get some useful numbers
    [n_electrodes, n_time_pts, n_conds, n_subs] = size(reduced_data);
    df_effect = n_conds - 1;
    
    if strcmpi(method, 'lb')
        
        %Lower bound on epsilon
        epsilon = (1/df_effect) * ones(n_electrodes, n_time_pts);
    
    else
        
        %Greenhouse-Geisser estimate
        epsilon = NaN([n_electrodes, n_time_pts]);
        for e = 1:n_electrodes
            for t = 1:n_time_pts
                epsilon(e,t) = epsGG(squeeze(reduced_data(e, t, :, :))');
            end
        end

        %Hyunh-Feldt
        if strcmpi(method, 'hf')
            epsilon = (n_subs*df_effect*epsilon - 2) ./ (df_effect*(n_subs-1-df_effect*epsilon));
            epsilon(epsilon>1) = 1;
        end
        
    end

end
