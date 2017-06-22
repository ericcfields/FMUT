%Calculate the interaction residuals for thei highest order interaction in
%the data. If the data has two factors, the  main effects will be
%subtracted. If the data has three factors, the main effects and two-way
%interactions will be subtracted.
%
%REQUIRED INPUTS
% data          - An electrode x time points x conditions x subjects array of ERP
%                 data. Array will vary in number of dimensions based on how many
%                 factors there are
%OUTPUT
% int-res       - interaction residuals: the data with all lower order
%                 effects subtracted out
%
%
%VERSION DATE: 22 June 2017
%AUTHOR: Eric Fields, Tufts University (Eric.Fields@tufts.edu)
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

function int_res = get_int_res(data)

    if ndims(data) ==5
        [~, ~, n_conds_A, n_conds_B, ~] = size(data);
        %Subtract main effects within each subject so that the data is 
        %exchangeable under the null hypothesis for the interaction
        int_res = NaN(size(data));
        for p = 1:n_conds_A
            for q = 1:n_conds_B
                int_res(:,:, p, q, :) = data(:,:,p,q,:) ...
                                         - mean(data(:,:,p,:,:), 4) ...
                                         - mean(data(:,:,:,q,:), 3) ...
                                         + mean(mean(data, 3), 4);
            end
        end
        assert(abs(mean(int_res(:))) < 1e-9);
    elseif ndims(data) == 6
        [~, ~, n_conds_A, n_conds_B, n_conds_C, ~] = size(data);
        %Subtract main effects then two-way effects within each subject so that 
        %the data is exchangeable under the null hypothesis for the three-way
        %interaction
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
        %Subtract two-way effects
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
        assert(abs(mean(int_res(:))) < 1e-9)
    end
    
end
