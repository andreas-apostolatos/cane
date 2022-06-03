function writeOutMultipatchBSplineSurface4GiD ...
    (BSplinePatches, pathToOutput, caseName)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Writes out a single/multipatch NURBS surface to be read and visualized by 
% GiD.
%
%            Input :
%   BSplinePatches : A cell array of B-Spline patches each of which
%                    containing the following information
%                       .p,.q : Polynomial orders in both parametric
%                               directions
%                    .Xi,.Eta : Knot vectors in both parametric directions
%                         .CP : Set of Control Point coordinates and
%                               weights
%                    .isNURBS : Flag of whether the basis is a NURBS or a
%                               B-Spline
%     pathToOutput : Absolute path to the folder where to write the GiD
%                    file
%         caseName : The name of the GiD file
%
%           Output :
%                     Writes out information in a GiD gile at
%                     pathToOutput\caseName.gid
%
% Function layout :
%
% 0. Read input
%
% 1. Write the preamble onto the file
%
% 2. Loop over all the patches of the multipatch geometry
% ->
%    2i. Get patch parameters
%
%   2ii. Loop over the trimming curves in the outer loop
%   ->
%        2ii.1. Get fixed and running parameters on the trimming curve
%
%        2ii.2. Initialize parameters and counters
%
%        2ii.3. Write all the points onto the trimming curve
%        ->
%               2ii.3i. Write the parametric coordinates onto the file
%
%              2ii.3ii. Update point counter
%
%             2ii.3iii. Update surface coordinates
%               
%              2ii.3iv. Compute and write the Cartesian coordinates on the B-Spline surface patch
%        <-
%       
%        2ii.4. Write all line segments onto the trimming curves
%        ->
%               2ii.4i. Write line segment
%
%              2ii.4ii. Update segment's counter
%        <-
%
%        2ii.5. Write the last segment of the outer loop
%   <-
%
%  2iii. Write out the surface header and update the surface counter
%
%   2iv. Write out the number of the trimming curves for the patch, their IDs and their orientation
%
%    2v. Write the center point of the patch and the corresponding normal to the surface vector at the center of the patch point
%
%   2vi. Write number of Control Points in xi- and eta- directions as well as polynomial orders in both parametric directions
%
%  2viii. Write the knot vectors
%
%    2ix. Write the Control Point weights
%
%     2x. Write the epilogue of the file
% <-
%
%% Function main body

%% 0. Read input

% No. patches
noPatches = length(BSplinePatches);

% No. points for the linearization of the trimming curves
noPntsOnTrCrvs = 49;

% Make directory to write out the results of the analysis
isExistent = exist(strcat(pathToOutput,caseName,'.gid'),'dir');
if ~isExistent
    mkdir(strcat(pathToOutput,caseName,'.gid'));
end
fileHandle = fopen(strcat(pathToOutput,caseName,'.gid','/',caseName,'.geo'),'w');

% Initialization of general counters
counterPoint = 1;
counterCurve = 1;
counterSurface = 1;

%% 1. Write the preamble onto the file
fprintf(fileHandle,'RAMSAN-ASCII-gid-v7.6\n');
fprintf(fileHandle,'UNKNOWN 0\n');
fprintf(fileHandle,'0\n');
fprintf(fileHandle,'1 Layer0 0 1 153 153 153\n');
fprintf(fileHandle,'0\n');
fprintf(fileHandle,'0\n\n');

%% 2. Loop over all the patches of the multipatch geometry
for counterPatches = 1:noPatches
    %% 2i. Get patch parameters
    p = BSplinePatches{counterPatches}.p;
    q = BSplinePatches{counterPatches}.q;
    Xi = BSplinePatches{counterPatches}.Xi;
    Eta = BSplinePatches{counterPatches}.Eta;
    CP = BSplinePatches{counterPatches}.CP;
    isNURBS = BSplinePatches{counterPatches}.isNURBS;
    nxi = length(CP(:,1,1));
    neta = length(CP(1,:,1));
    counterCurveBeginning = counterCurve;
    firstLastPnt = counterCurve;
    
    %% 2ii. Loop over the trimming curves in the outer loop
    for counterTrCrvs = 1:4
        %% 2ii.1. Get fixed and running parameters on the trimming curve
        if counterTrCrvs == 1
            isOnXi = true;
            xiEtaStart = Xi(1);
            xiEtaEnd = Xi(end);
            xiEtaFixed = Eta(1);
            xi = xiEtaStart;
            eta = xiEtaFixed;
        elseif counterTrCrvs == 2
            isOnXi = false;
            xiEtaStart = Eta(1);
            xiEtaEnd = Eta(end);
            xiEtaFixed = Xi(end);
            xi = xiEtaFixed;
            eta = xiEtaStart;
        elseif counterTrCrvs == 3
            isOnXi = true;
            xiEtaStart = Xi(end);
            xiEtaEnd = Xi(1);
            xiEtaFixed = Eta(end);
            xi = xiEtaStart;
            eta = xiEtaFixed;
        elseif counterTrCrvs == 4
            isOnXi = false;
            xiEtaStart = Eta(end);
            xiEtaEnd = Eta(1);
            xiEtaFixed = Xi(1);
            xi = xiEtaFixed;
            eta = xiEtaStart;
        end
        
        %% 2ii.2. Initialize parameters and counters
        dXiEtaRunning = (xiEtaEnd - xiEtaStart)/noPntsOnTrCrvs;
        startIndex = 1;
        endIndex = noPntsOnTrCrvs;
        if counterTrCrvs == 1
            startIndex = 0;
        end
        if counterTrCrvs == 4
            endIndex = noPntsOnTrCrvs - 1;
        end
%         firstLastPnt = 1;
        
        %% 2ii.3. Write all the points onto the trimming curve
        for i = startIndex:endIndex
            %% 2ii.3i. Write the parametric coordinates onto the file
            fprintf(fileHandle,['1',' ',num2str(counterPoint),' ','0 0 2 0 0 2 0\n']);
            
            %% 2ii.3ii. Update point counter
            counterPoint = counterPoint + 1;
            
            %% 2ii.3iii. Update surface coordinates
            if isOnXi
                xi = xiEtaStart + i * dXiEtaRunning;
            else
                eta = xiEtaStart + i * dXiEtaRunning;
            end
            
            %% 2ii.3iv. Compute and write the Cartesian coordinates on the B-Spline surface patch
            xiSpan = findKnotSpan(xi,Xi,nxi);
            etaSpan = findKnotSpan(eta,Eta,neta);
            R = computeIGABasisFunctionsAndDerivativesForSurface...
                (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,0);
            X = computeCartesianCoordinatesOfAPointOnBSplineSurface...
                (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,R);
            fprintf(fileHandle,'%.15f %.15f %.15f\n',X(1,1),X(2,1),X(3,1));
            fprintf(fileHandle,'\n');
        end
        
        %% 2ii.4. Write all line segments onto the trimming curves
        for i = 1:noPntsOnTrCrvs - 1
            %% 2ii.4i. Write line segment
            fprintf(fileHandle,['2',' ',num2str(counterCurve),' ','0 0 1 0 0 2 0\n']);
            fprintf(fileHandle,[num2str(counterCurve),' ',num2str(counterCurve + 1),'\n']);
            
            %% 2ii.4ii. Update segment's counter
            counterCurve = counterCurve + 1;
        end
        
        %% 2ii.5. Write the last segment of the outer loop
        if counterTrCrvs ~= 4
            fprintf(fileHandle,['2',' ',num2str(counterCurve),' ','0 0 1 0 0 2 0\n']);
            fprintf(fileHandle,[num2str(counterCurve),' ',num2str(counterCurve + 1),'\n']);
            counterCurve = counterCurve + 1;
        else
            fprintf(fileHandle,['2',' ',num2str(counterCurve),' ','0 0 1 0 0 2 0\n']);
            fprintf(fileHandle,[num2str(counterCurve),' ',num2str(firstLastPnt),'\n']);
            counterCurve = counterCurve + 1;
        end
        fprintf(fileHandle,'\n');
    end
    
    %% 2iii. Write out the surface header and update the surface counter
    fprintf(fileHandle,['14' ' ' num2str(counterSurface) ' ' '0 0 0 0 0 2 0\n']);
    counterSurface = counterSurface + 1;
    
    %% 2iv. Write out the number of the trimming curves for the patch, their IDs and their orientation
    
    % Write the number per patch
    noCurvesPerPatch = counterCurve - counterCurveBeginning;
    fprintf(fileHandle,[num2str(noCurvesPerPatch) '\n']);
    
    % Write the IDs of the trimming curves
    for iCurve = counterCurveBeginning:counterCurve-1
        fprintf(fileHandle,[num2str(iCurve) ' ']);
    end
    fprintf(fileHandle,'\n');
    
    % write the orientation of the curves
    for iCurve = 1:noCurvesPerPatch
        fprintf(fileHandle,['0' ' ']);
    end
    fprintf(fileHandle,'\n');
    
    %% 2v. Write the center point of the patch and the corresponding normal to the surface vector at the center of the patch point
    xiMiddle = (Xi(end) + Xi(1))/2;
    etaMiddle = (Eta(end) + Eta(1))/2;
    xiSpanMiddle = findKnotSpan(xiMiddle,Xi,nxi);
    etaSpanMiddle = findKnotSpan(etaMiddle,Eta,neta);
    dRMiddle = computeIGABasisFunctionsAndDerivativesForSurface...
        (xiSpanMiddle,p,xiMiddle,Xi,etaSpanMiddle,q,etaMiddle,Eta,CP,...
        isNURBS,1);
    XMiddle = computeCartesianCoordinatesOfAPointOnBSplineSurface...
        (xiSpanMiddle,p,xiMiddle,Xi,etaSpanMiddle,q,etaMiddle,Eta,CP,dRMiddle(:,1));
    surfaceNormalMiddle = computeBaseVectorsAndDerivativesForBSplineSurface...
        (xiSpanMiddle,p,etaSpanMiddle,q,CP,0,dRMiddle);
    fprintf(fileHandle,'%.15f %.15f %.15f\n',XMiddle(1,1),XMiddle(2,1),XMiddle(3,1));
    fprintf(fileHandle,'%.15f %.15f %.15f\n',surfaceNormalMiddle(1,1),surfaceNormalMiddle(2,1),surfaceNormalMiddle(3,1));
    
    %% 2vi. Write number of Control Points in xi- and eta- directions as well as polynomial orders in both parametric directions
    fprintf(fileHandle,['0' ' ' num2str(nxi) ' ' num2str(neta) ' ' num2str(p) ' ' num2str(q) '\n']);
    
    %% 2vii. Write the coordinates of the Control Points belonging to the patch
    for iEtaCounter = 1:neta
        for iXiCounter = 1:nxi
            xCP = CP(iXiCounter,iEtaCounter,1);
            yCP = CP(iXiCounter,iEtaCounter,2);
            zCP = CP(iXiCounter,iEtaCounter,3);
            fprintf(fileHandle,'%.15f %.15f %.15f\n',xCP,yCP,zCP);
        end
    end
    
    %% 2viii. Write the knot vectors
    for iXiCounter = 1:length(Xi)
        fprintf(fileHandle,[num2str(Xi(iXiCounter)) ' ']);
    end
    fprintf(fileHandle,'\n');
    for iEtaCounter = 1:length(Eta)
        fprintf(fileHandle,[num2str(Eta(iEtaCounter)) ' ']);
    end
    fprintf(fileHandle,'\n');
    
    %% 2ix. Write the Control Point weights
    fprintf(fileHandle,['1' ' ']);
    for iEtaCounter = 1:neta
        for iXiCounter = 1:nxi
            fprintf(fileHandle,[num2str(CP(iXiCounter,iEtaCounter,4)) ' ']);
        end
    end
    
    %% 2x. Write the epilogue of the file
    fprintf(fileHandle,'\n\n');
    fprintf(fileHandle,'0');
end

end
