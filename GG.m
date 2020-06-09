%Calculate the Greenhouse-Geisser estimate of epsilon at all time points and
%electrodes
%
%EXAMPLE USAGE
% >> epsilon = GG(data, [], 3);
%
%REQUIRED INPUTS
% data          - An electrode x time points x conditions x subjects array of ERP
%                 data. Array will vary in number of dimensions based on how many
%                 factors there are.
% cond_subs     - Array giving the number of subjects in each condition of
%                 the between subjects factor. For example, if cond_subs is
%                 [8, 9], then there should be 17 subjects with the first 8
%                 being in condition A and the next 9 being in condition B.
%                 For fully within-subjects designs cond_subs = []
% dims          - Dimensions of the data array involved in the effect to be
%                 calculated. For example, if data is an electrode x time points
%                 x Factor A x Factor B x subjects array and you want to
%                 calculated the main effect of A, dims = 3. If you want to
%                 calculate the AxB interaciton, dims  = [3, 4].
%
%OUTPUT
% epsilon      - electrode x time point array of GG estimate of epsilon
%
%
%VERSION DATE: 9 June 2020
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 

%Copyright (c) 2020, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function epsilon = GG(data, cond_subs, dims)
    
    %Get data reduced for analysis
    reduced_data = reduce_data(data, dims);
    
    %Get some useful numbers
    n_electrodes = size(reduced_data, 1);
    n_time_pts = size(reduced_data, 2);
    
    %This function currently only works for one-way fully repeated measures
    %designs
    if ~isempty(cond_subs) && ~isequal(cond_subs, 0) && length(cond_subs) > 1
        error('Greenhouse-Geisser correction is currently only implemented for fully repetaed measures ANOVA');
    end
    if ndims(reduced_data) ~= 4
        error('Greenhouse-Geisser correction is currently only implemented for one-way RM ANOVA');
    end
    
    %Calculate epsilon
    epsilon = NaN([n_electrodes, n_time_pts]);
    for e = 1:n_electrodes
        for t = 1:n_time_pts
            epsilon(e,t) = epsGG(squeeze(reduced_data(e, t, :, :))');
        end
    end

end
