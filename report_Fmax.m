%Summarize results of Fmax test at the command window
%
%EXAMPLE USAGE
%
% report_Fmax(GND.F_tests(1))
%
%REQUIRED INPUTS
% results        - An F_tests results struct
%
%OPTIONAL INPUTS
% effects        - An array with each row specifying the dimensions of the
%                  data array that are involved in each effect {default:
%                  calculated from the results struct}
% effects labels - A matching cell array with strings describing each effect
%                  {default: calculated from the results struct}
%
%VERSION DATE: 14 July 2014
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

%%%%%%%%%%%%%%%%%%%  REVISION LOG   %%%%%%%%%%%%%%%%%%%
% 7/14/17  - Moved to separate function from FmaxGND and FmaxGRP

function report_Fmax(results, effects, effects_labels)

    if nargin == 1
        [effects, effects_labels] = get_effects(results.factors);
        assert(all(strcmpi(fieldnames(results.F_obs), effects_labels)))
    end

    fprintf('\n##### RESULTS #####\n');
    for i = 1:length(effects)
        fprintf('\n%s effect\n', effects_labels{i});
        if length(effects) == 1
            if any(results.null_test(:))
                fprintf('Critical F-value: %.3f\n', results.F_crit);
                fprintf('That corresponds to a test-wise alpha level of %.3f.\n', ...
                        fcdf(results.F_crit, results.df(1), results.df(2), 'upper'));
                if strcmpi(p.Results.mean_wind, 'yes') || strcmpi(p.Results.mean_wind, 'y')
                    for t = 1:size(time_wind, 1)
                        fprintf('Significant electrodes for time window %d - %d: ', time_wind(t, 1), time_wind(t, 2));
                        fprintf('%s ', results.include_chans{results.null_test(:, t)});
                        fprintf('\n');
                    end
                else
                    fprintf('Electrodes and time points with significant effects:\n');
                    for t = 1:length(results.used_tpt_ids)
                        if any(results.null_test(:, t))
                            fprintf('%d ms, electrode(s): ', GRP.time_pts(results.used_tpt_ids(t)));
                            fprintf('%s ', results.include_chans{results.null_test(:, t)});
                            fprintf('\n');
                        end
                    end
                end
                fprintf('All significant corrected p-values are between %f and %f.\n', ... 
                         max(results.adj_pval(results.adj_pval <= .05)), ... 
                         min(results.adj_pval(:)));
            else
                fprintf('NO significant time points or electrodes.\n');
            end
        else
            if any(results.null_test.(effects_labels{i})(:))
                fprintf('Critical F-value: %.3f\n', results.F_crit.(effects_labels{i}));
                fprintf('That corresponds to a test-wise alpha level of %.3f.\n', ...
                        fcdf(results.F_crit.(effects_labels{i}), results.df.(effects_labels{i})(1), results.df.(effects_labels{i})(2), 'upper'));
                if strcmpi(p.Results.mean_wind, 'yes') || strcmpi(p.Results.mean_wind, 'y')
                    for t = 1:size(time_wind, 1)
                        fprintf('Significant electrodes for time windonw %d - %d: ', time_wind(t, 1), time_wind(t, 2));
                        fprintf('%s ', results.include_chans{results.null_test.(effects_labels{i})(:, t)});
                        fprintf('\n');
                    end
                else
                    fprintf('Electrodes and time points with significant effects:\n');
                    for t = 1:length(results.used_tpt_ids)
                        if any(results.null_test.(effects_labels{i})(:, t))
                            fprintf('%d ms, electrode(s): ', GRP.time_pts(results.used_tpt_ids(t)));
                            fprintf('%s ', results.include_chans{results.null_test.(effects_labels{i})(:, t)});
                            fprintf('\n');
                        end
                    end
                end
                fprintf('All significant corrected p-values are between %f and %f.\n', ... 
                         max(results.adj_pval.(effects_labels{i})(results.adj_pval.(effects_labels{i}) <= .05)), ... 
                         min(results.adj_pval.(effects_labels{i})(:)));
            else
                fprintf('NO significant time points or electrodes.\n');
            end
        end
    end

end