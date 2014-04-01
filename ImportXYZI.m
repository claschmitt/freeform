%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IMPORTPTS
%
% Importing TLS *.pts files into a matrix
% data
% input
% fileName = filename
% pathName = pathname
% varargin = optional parameter bounding Box 2D
%     'boundingbox': Boundingbox with a matrix 5x2, each row is one edge point [x y]
%     'headerline' : number of lines from top which not consists of xyzi
%     'sort2profil': logical if the pc should be sorted into different profiles
%     'delimiter'  : delimiter in the text file
% output
% scan = scan structures
%
% claudius schmitt, 05.2013
% schmitt@gih.uni-hannover.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test



function [scan] = ImportXYZI( absfilePathName, varargin )

% %% test
% clear all;
% varargin{1} = [-2.9 5.05 ; -2.9 5.25 ; 14.0 4.7 ; 14.0 4.5 ; -2.9 5.05 ];
% fileName = '20120412_10_23.xyz.asc';
% fileName = 'profil3.xyz.txt';
% absfilePathName = 'E:\Projekte\20120412_Rethen\Scans\Export\';
% absfilePathName = 'E:\Projekte\Bachelorprojekt2013\PC\Standpunkt1 - Kopie.pts';
%
% load('E:\tmp\newData2.mat');
% numHeaderLines = 0;
% sort2Profil = true;
% varargin{1} = [-2.9 5.05 ; -2.9 5.25 ; 14.5 4.7 ; 14.5 4.5 ; -2.9 5.05 ];

[filePath, fileName, ext] = fileparts(absfilePathName);  
%% define Struct
scan = struct('scanFile' , struct('filePath',[fileName ext] , 'fileName', filePath), 'data', cell(1), 'roiGrid', cell(1), 'boundingBox', []);

%%set default parameters
delimiter = ' ';
headerline = 0;
boundingBox = [];
sort2Profil = false;
%% Read varargin

if isempty(varargin)
    
else
    
    indexBoundingBox = strcmp('boundingbox', varargin);
    if isempty(find(indexBoundingBox,1))
        boundingBox = [];
    else
        boundingBox = varargin{find(indexBoundingBox,1) +1};
    end
    
    indexHeaderlines = strcmp('headerline', varargin);
    if isempty(find(indexHeaderlines,1))
        headerline = 0;
    else
        headerline = varargin{find(indexHeaderlines,1) +1};
    end
    
    indexSort2Profile = strcmp('sort2profil', varargin);
    if isempty(find(indexSort2Profile,1))
        sort2Profil = false;
    else
        sort2Profil = varargin{find(indexSort2Profile,1) +1};
    end
    
    indexDelimiter = strcmp('delimiter', varargin);
    if isempty(find(indexDelimiter,1))
        delimiter = ' ';
    else
        delimiter = varargin{find(indexDelimiter,1) +1};
    end   
    
end
%% Import the file


fid = fopen(absfilePathName,'r+');

startFReading = tic;
newData = textscan(fid, '%f %f %f %f %f %f %f%*[^\n]','headerlines', headerline, 'delimiter', delimiter);
% newData = textscan(fid, '%f %f %*[^\n]','headerlines', headerline, 'delimiter', delimiter);
timeFReading = toc(startFReading)

fclose (fid);

%% BoundixBox
if isempty(boundingBox)
else
    inside = inpolygon(newData{2}, newData{3}, boundingBox(:,1), boundingBox(:,2));
    scan.boundingBox = boundingBox;
    newData{1} = newData{1} (inside,1);
    newData{2} = newData{2} (inside,1);
    newData{3} = newData{3} (inside,1);
    newData{4} = newData{4} (inside,1);
    newData{5} = newData{5} (inside,1);
    newData{6} = newData{6} (inside,1);
    newData{7} = newData{7} (inside,1);
    %     test
    %     plot(newData{2},newData{3},'+');
    %     hold on;
    %     plot(boundingBox(:,1), boundingBox(:,2),'r');
    %     save('E:\tmp\newData.mat', 'newData', '-v7.3');
    
    %         % further binary selection - Messung 1
    %         indexSelect1 = newData{2} > -2.90;
    %         indexSelect2 = newData{2} < 10.00;
    %         indexSelect3 = newData{2} > 11.15;
    %         indexSelect4 = newData{2} < 11.30;
    %         indexSelect5 = newData{2} > 12.45;
    %         indexSelect6 = newData{2} < 14.50;
    %         indexAll = indexSelect1 & indexSelect2 | indexSelect3 & indexSelect4 | indexSelect5 & indexSelect6;
    %
    %         newData{1} = newData{1} (indexAll,1);
    %         newData{2} = newData{2} (indexAll,1);
    %         newData{3} = newData{3} (indexAll,1);
    %         newData{4} = newData{4} (indexAll,1);
    
    %        test
    %          plot(boundingBox(:,1), boundingBox(:,2), newData{2}(indexAll,1),newData{3}(indexAll,1),'+');
    
end


%% Sort Points 2 Profile
if sort2Profil
    
    numPoints = length(newData{2});
    numProfilePoints = 1;
    numProfile = 0;
    pointPrev = 0;
    startProfile = 1;
    endProfile   = 0;
    
    %% reconfigure data
    
    minY = min(newData{2}(:,1));
    minZ = min(newData{3}(:,1));
    scan.roiGrid {1} = [0 abs(minY)  minZ];
    scan.boundingBox = boundingBox;
    
    if minY < 0.0
%         newData{2} = newData{2}(:,1) + abs(minY);
        newData{2} = newData{2}(:,1) + abs(min(boundingBox(:,1)));
    end
    
    startProfilSeparation = tic;
    
    endPotentialProfile = find(newData{2} <= 0.01);
    
    for i=1 : length(endPotentialProfile)
        if (endPotentialProfile(i)+1 <= length(newData{2}))
            if ((newData{2}(endPotentialProfile(i)+1, 1) - newData{2}(endPotentialProfile(i), 1)) >= 1.00)
                endProfile = endPotentialProfile(i);
                numProfile = numProfile + 1;
                
                scan.data{numProfile,1} (:,2) = newData{2}(startProfile:endProfile,1); % Y Koordinate
                scan.data{numProfile,1} (:,3) = newData{3}(startProfile:endProfile,1); % Z Koordinate
                scan.data{numProfile,1} (:,4) = newData{4}(startProfile:endProfile,1); % Intensity
                
                if newData{2}(startProfile,1) > newData{2}(endProfile,1)
                    scan.data{numProfile,1}  = flipud(scan.data{numProfile,1});
                else
                    
                end
                
                startProfile = endProfile + 1;
            end
        else
            endProfile = endPotentialProfile(i);
            numProfile = numProfile + 1;
            
            scan.data{numProfile,1} (:,2) = newData{2}(startProfile:endProfile,1); % Y Koordinate
            scan.data{numProfile,1} (:,3) = newData{3}(startProfile:endProfile,1); % Z Koordinate
            scan.data{numProfile,1} (:,4) = newData{4}(startProfile:endProfile,1); % Intensity
            
            if newData{2}(startProfile,1) > newData{2}(endProfile,1)
                scan.data{numProfile,1}  = flipud(scan.data{numProfile,1});
            else
                
            end
            
        end
        
    end
    
    timeProfileSeparation = toc(startProfilSeparation)
    
    
    
    % for i=1 : numPoints
    %     if ((newData{2}(i,1) - pointPrev) >= 1.000) && (pointPrev >= 0)
    %         numProfile = numProfile + 1;
    %         numProfilePoints = 1;
    %     end
    %
    %     scan.data{numProfile,1} (numProfilePoints,1) = newData{2}(i,1); % Y Koordinate
    %     scan.data{numProfile,1} (numProfilePoints,2) = newData{3}(i,1); % Z Koordinate
    %     scan.data{numProfile,1} (numProfilePoints,3) = newData{4}(i,1); % Intensity
    %     numProfilePoints = numProfilePoints + 1;
    %     pointPrev = newData{2}(i,1);
    %
    % end
    
else
    minX = min(newData{1}(:,1));
    minY = min(newData{2}(:,1));
    minZ = min(newData{3}(:,1));
    scan.roiGrid {1} = [minX abs(minY)  minZ];
    scan.data{1}(:,1) = newData{1};
    scan.data{1}(:,2) = newData{2};
    scan.data{1}(:,3) = newData{3};
    scan.data{1}(:,4) = newData{4};
    scan.data{1}(:,5) = newData{5};
    scan.data{1}(:,6) = newData{6};
    scan.data{1}(:,7) = newData{7};
    scan.data{1}      = scan.data{1}(:,~isnan(scan.data{1}(1,:))); 
end



