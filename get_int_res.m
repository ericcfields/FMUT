%Calculate the interaction residuals for the highest order interaction in
%the data. If the data has two factors, the  main effects will be
%subtracted. If the data has three factors, the main effects and two-way
%interactions will be subtracted.
%
%EXAMPLE USAGE
% >> int_res = get_int_res(data, [], [3, 4])
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

function int_res = get_int_res(data, cond_subs, dims)

    if isempty(cond_subs) || length(cond_subs)==1 || ~any(dims == ndims(data))
        
        %####################################
        %### WITHIN SUBJECTS FACTORS ONLY ###
        %####################################
        
        %Interaction of two within-subjects factors
        if ndims(data) == 5 && isequal(dims, [3, 4])
            int_res = twoway_rb_main(data);
        %Interaction of three within-subjects factors
        elseif ndims(data) == 6 && isequal(dims, [3, 4, 5])
            int_res1 = threeway_rb_main(data);
            int_res = threeway_rb_int(int_res1);
        end
  
    else
        
        %####################################
        %######## SPLIT PLOT DESIGNS ########
        %####################################
       
        if ndims(data) == 4
            %Interaction of one within-subjects factor and one between-subjects factor
            int_res = twoway_sp_main(data, cond_subs);
        elseif ndims(data) == 5
            if isequal(dims, [3, 4])
                %Interaction of two within-subjects in the presence of a between-subjects factor
                int_res = twoway_rb_main(data, cond_subs);
            elseif isequal(dims, [3, 4, 5])
                %Interaction of two within-subjects factors and one between-subjects factor
                int_res1 = threeway_sp_main(data, cond_subs);
                int_res = threeway_sp_int(int_res1, cond_subs);
            end
        end
        
    end
    
    assert(abs(mean(int_res(:))) < 1e-9)
    
end


function int_res = twoway_rb_main(data)
%Subtract main effects from two-way within-subjects design
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
end

function int_res = threeway_rb_main(data)
%Subtract main effects from three-way within-subjects design
    [~, ~, n_conds_A, n_conds_B, n_conds_C, ~] = size(data);
    int_res = NaN(size(data));
    %Subtract main effects
    for p = 1:n_conds_A
        for q = 1:n_conds_B
            for r = 1:n_conds_C
                int_res(:,:,p,q,r,:) = data(:,:,p,q,r,:) ... 
                                        - mean(mean(data(:,:,p,:,:,:), 4), 5) ...
                                        - mean(mean(data(:,:,:,q,:,:), 3), 5) ...
                                        - mean(mean(data(:,:,:,:,r,:), 3), 4) ...
                                        + 2*mean(mean(mean(data, 3), 4), 5);
            end
        end
    end
end


function int_res = threeway_rb_int(data)
%Subtract two-way interaction effects from three-way within-subjects design
    [~, ~, n_conds_A, n_conds_B, n_conds_C, ~] = size(data);
    int_res = NaN(size(data));
    for p = 1:n_conds_A
        for q = 1:n_conds_B
            for r = 1:n_conds_C
                int_res(:,:,p,q,r,:) = data(:,:,p,q,r,:) ... 
                                        - mean(data(:,:,p,q,:,:), 5) ...
                                        - mean(data(:,:,p,:,r,:), 4) ...
                                        - mean(data(:,:,:,q,r,:), 3) ...
                                        + 2*mean(mean(mean(data, 3), 4), 5);
            end
        end
    end
end


function int_res = twoway_sp_main(data, cond_subs)
%Subtract main effects from two-way split plot design
    n_conds_A = length(cond_subs);
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
end


function int_res = threeway_sp_main(data, cond_subs)
%Subtract main effects from three-way split plot design 
    n_conds_A = length(cond_subs);
    [~, ~, n_conds_B, n_conds_C, ~] = size(data);
    int_res = NaN(size(data));
    for p = 1:n_conds_A
        first = sum(cond_subs(1:p)) - cond_subs(p) + 1;
        for q = 1:n_conds_B
            for r = 1:n_conds_C
                for n = 1:cond_subs(p)
                    sub = first+n-1;
                    int_res(:, :, q, r, sub) = data(:, :, q, r, sub) ...
                                                - mean(mean(data(:, :, q, :, :), 4), 5) ...
                                                - mean(mean(data(:, :, :, r, :), 3), 5) ...
                                                - mean(mean(data(:, :, :, :, sub), 3), 4) ...
                                                + 2*mean(mean(mean(data, 3), 4), 5);
                end
            end
        end
    end
end


function int_res = threeway_sp_int(data, cond_subs)
%Subtract two-way interaction from three-way split plot design
    n_conds_A = length(cond_subs);
    [~, ~, n_conds_B, n_conds_C, ~] = size(data);
    int_res = NaN(size(data));
    for p = 1:n_conds_A
        first = sum(cond_subs(1:p)) - cond_subs(p) + 1;
        for q = 1:n_conds_B
            for r = 1:n_conds_C
                for n = 1:cond_subs(p)
                    sub = first+n-1;
                    int_res(:, :, q, r, sub) = data(:, :, q, r, sub) ...
                                               - mean(data(:, :, q, r, :), 5) ...
                                               - mean(data(:, :, q, :, sub), 4) ...
                                               - mean(data(:, :, :, r, sub), 3) ...
                                               + 2 * mean(mean(mean(data, 3), 4), 5);
                end
            end
        end
    end
end