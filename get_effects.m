%Find all effects for a given repeated-measures ANOVA design.
%
%EXAMPLE USAGE
% >> [effects, effects_labels] = get_effects({'Frequency', 'Priming'})
%
%REQUIRED INPUTS
% factor_names   - A cell array of strings of the names of all factors in
%                  the design
%
%OUTPUTS
% effects        - An array with each row specifying the dimensions of the
%                  data array that are involved in each effect
% effects labels - A matching cell array with strings describing each effect
%
%VERSION DATE: 28 November 2016
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function [effects, effects_labels] = get_effects(factor_names)
%Given factor names, return dimensions involved and effect names for all effects
    
    %Get numerical coding of all main effects and interactions
    effects = {};
    for i = 1:length(factor_names)
        effects = [effects; num2cell(nchoosek(1:length(factor_names), i), 2)]; %#ok<AGROW>
    end
    
    %Get effect names for all the effects described above
    effects_labels = cell(size(effects));
    for i = 1:length(effects)
        if verLessThan('matlab','8.1')
            %Workaround for the fact that strjoin doesn't exist in older versions of MATLAB
            effects_labels{i} = '';
            for j = 1:length(effects{i})
                if j < length(effects{i})
                    effects_labels{i} = [effects_labels{i} factor_names{effects{i}(j)} 'X'];
                else
                    effects_labels{i} = [effects_labels{i} factor_names{effects{i}(j)}];
                end
            end
        else
            effects_labels{i} = strjoin(factor_names(effects{i}), 'X');            
        end
    end
    
end
