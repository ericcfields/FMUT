%Add Apache POI library to MATLAB's static Java class path
%
%EXAMPLE USAGE
% >> add_poi_path
%
%VERSION DATE: 18 November 2017
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function add_poi_path()
    
    %POI path
    poi_files = {fullfile(fileparts(which('add_poi_path')), 'poi_library/poi-3.8-20120326.jar')
                 fullfile(fileparts(which('add_poi_path')), 'poi_library/poi-ooxml-3.8-20120326.jar')
                 fullfile(fileparts(which('add_poi_path')), 'poi_library/poi-ooxml-schemas-3.8-20120326.jar')
                 fullfile(fileparts(which('add_poi_path')), 'poi_library/xmlbeans-2.3.0.jar\n')
                 fullfile(fileparts(which('add_poi_path')), 'poi_library/dom4j-1.6.1.jar')
                 fullfile(fileparts(which('add_poi_path')), 'poi_library/stax-api-1.0.1.jar')};
    
    %POI paths not currently on the static Java class path
    missing_poi_files = poi_files(~ismember(poi_files, javaclasspath('-static')));

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
    