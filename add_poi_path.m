%Add Apache POI library to MATLAB's static Java class path
%
%EXAMPLE USAGE
% >> add_poi_path
%
%VERSION DATE: 19 November 2017
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function add_poi_path()
    
    %Get full path for POI .jar files not on static Java class path
    [~, missing_poi_files] = get_poi_paths();

    %Append POI locations that are not already on the static Java class
    %path
    if ~isempty(missing_poi_files)
        %Full path for javaclasspath file
        static_file = [prefdir filesep 'javaclasspath.txt'];
        %Append missing POI files to javaclasspath.txt file
        f_out = fopen(static_file, 'a');
        for i = 1:length(missing_poi_files)        
            fprintf(f_out, '%s\n', missing_poi_files{i});
        end
        fclose(f_out);
    end
    
end
    