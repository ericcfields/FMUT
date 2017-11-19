%Return the full paths of the .jar files of the Apache POI library
%
%EXAMPLE USAGE
% >> add_poi_path
%
%OUTPUT
% poi_files          - Full paths for all .jar files in the Apache POI
%                      library
% missing_poi_files  - Full paths for .jar files in the Apache POI library
%                      that are not on the static Java class path
%
%VERSION DATE: 19 November 2017
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function [poi_files, missing_poi_files] = get_poi_paths()

    poi_files = {fullfile(fileparts(which('add_poi_path')), 'poi_library/poi-3.8-20120326.jar')
                 fullfile(fileparts(which('add_poi_path')), 'poi_library/poi-ooxml-3.8-20120326.jar')
                 fullfile(fileparts(which('add_poi_path')), 'poi_library/poi-ooxml-schemas-3.8-20120326.jar')
                 fullfile(fileparts(which('add_poi_path')), 'poi_library/xmlbeans-2.3.0.jar')
                 fullfile(fileparts(which('add_poi_path')), 'poi_library/dom4j-1.6.1.jar')
                 fullfile(fileparts(which('add_poi_path')), 'poi_library/stax-api-1.0.1.jar')};
             
     missing_poi_files = poi_files(~ismember(poi_files, javaclasspath('-static')));

end