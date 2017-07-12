%Calculate the interaction residuals for the highest order interaction in
%the data. If the data has two factors, the  main effects will be
%subtracted. If the data has three factors, the main effects and two-way
%interactions will be subtracted.
%
%REQUIRED INPUTS
% data          - An electrode x time points x conditions x subjects array of ERP
%                 data. Array will vary in number of dimensions based on how many
%                 factors there are
% cond_subs     - Array giving the number of subjects in each condition of
%                 the between subjects factor. For example, if cond_subs is
%                 [8, 9], then there should be 17 subjects with first 8
%                 being in condition A and the next 9 being in condition B
%                 [] indicates a within-subjects ANOVA
% dims          - Dimensions of the data array involved in the effect
%
%OUTPUT
% int_res       - interaction residuals: the data with all lower order
%                 effects subtracted out
%
%
%VERSION DATE: 12 July 2017
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 
%This is a beta version of this software. It needs additional testing 
%and SHOULD NOT be considered error free.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

%%%%%%%%%%%%%%%%%%%  REVISION LOG   %%%%%%%%%%%%%%%%%%%
% 6/22/17   - First version. Code re-organized from other functions.
% 7/12/17   - Now handles betwee-subjects factors

function int_res = get_int_res(data, cond_subs, dims)

    if isempty(cond_subs) || length(cond_subs)==1 || ~any(dims == ndims(data))
        
        %####################################
        %### WITHIN_SUBJECTS FACTORS ONLY ###
        %####################################
        
        %--- Interaction of two within-subjects factors ---%
        if ndims(data) == 5 && isequal(dims, [3, 4])
            [~, ~, n_conds_A, n_conds_B, ~] = size(data);
            int_res = NaN(size(data));
            %Subtract main effects
            for p = 1:n_conds_A
                for q = 1:n_conds_B
                    int_res(:,:, p, q, :) = data(:,:,p,q,:) ...
                                             - mean(data(:,:,p,:,:), 4) ...
                                             - mean(data(:,:,:,q,:), 3) ...
                                             + mean(mean(data, 3), 4);
                end
            end
         
        %--- Interaction of three within-subjects factors ---%
        elseif ndims(data) == 6 && isequal(dims, [3, 4, 5])
            [~, ~, n_conds_A, n_conds_B, n_conds_C, ~] = size(data);
            int_res1 = NaN(size(data));
            int_res = NaN(size(data));
            %Subtract main effects
            for p = 1:n_conds_A
                for q = 1:n_conds_B
                    for r = 1:n_conds_C
                        int_res1(:,:,p,q,r,:) = data(:,:,p,q,r,:) ... 
                                                - mean(mean(data(:,:,p,:,:,:), 4), 5) ...
                                                - mean(mean(data(:,:,:,q,:,:), 3), 5) ...
                                                - mean(mean(data(:,:,:,:,r,:), 3), 4) ...
                                                + 2*mean(mean(mean(data, 3), 4), 5);
                    end
                end
            end
            %Subtract two-way interaction effects
            for p = 1:n_conds_A
                for q = 1:n_conds_B
                    for r = 1:n_conds_C
                        int_res(:,:,p,q,r,:) = int_res1(:,:,p,q,r,:) ... 
                                                - mean(int_res1(:,:,p,q,:,:), 5) ...
                                                - mean(int_res1(:,:,p,:,r,:), 4) ...
                                                - mean(int_res1(:,:,:,q,r,:), 3) ...
                                                + 2*mean(mean(mean(int_res1, 3), 4), 5);
                    end
                end
            end
        end
        
    else
        
        %###############################################################
        %### INTERACTION EFFECTS INVOLVING A BETWEEN-SUBJECTS FACTOR ###
        %###############################################################
        
        n_conds_A = length(cond_subs);
        
        %--- Interaction of one within-subjects factor and one between-subjects factor ---%
        if ndims(data) == 4
            n_conds_B = size(data, 3);
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
        
        %--- Interaction of two within-subjects factors and on between-subjects factor ---%
        elseif ndims(data) == 5
            [~, ~, n_conds_B, n_conds_C, ~] = size(data);
            int_res = NaN(size(data));
            if isequal(dims, [3, 4])
                int_res = get_int_res(data);
            elseif isequal(dims, [3, 4, 5])
                int_res1 = NaN(size(data));
                int_res  = NaN(size(data));
                %Subtract main effects
                for p = 1:n_conds_A
                    first = sum(cond_subs(1:p)) - cond_subs(p) + 1;
                    for q = 1:n_conds_B
                        for r = 1:n_conds_C
                            for n = 1:cond_subs(p)
                                sub = first+n-1;
                                int_res1(:, :, q, r, sub) = data(:, :, q, r, sub) ...
                                                            - mean(mean(data(:, :, q, :, :), 4), 5) ...
                                                            - mean(mean(data(:, :, :, r, :), 3), 5) ...
                                                            - mean(mean(data(:, :, :, :, sub), 3), 4) ...
                                                            + 2*mean(mean(mean(data, 3), 4), 5);
                            end
                        end
                    end
                end
                %Subtract two-way interactions
                for p = 1:n_conds_A
                    first = sum(cond_subs(1:p)) - cond_subs(p) + 1;
                    for q = 1:n_conds_B
                        for r = 1:n_conds_C
                            for n = 1:cond_subs(p)
                                sub = first+n-1;
                                int_res(:, :, q, r, sub) = int_res1(:, :, q, r, sub) ...
                                                           - mean(int_res1(:, :, q, r, :), 5) ...
                                                           - mean(int_res1(:, :, q, :, sub), 4) ...
                                                           - mean(int_res1(:, :, :, r, sub), 3) ...
                                                           + 2 * mean(mean(mean(int_res1, 3), 4), 5);
                            end
                        end
                    end
                end
            end
        end
        
    end
    
    assert(abs(mean(int_res(:))) < 1e-9)
    
end
