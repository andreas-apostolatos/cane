function is_loaded = load_dependencies(path_prefix)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
% 
% Loads the dependecies for the isogeometric membrane analysis
%
%       Input :
% path_prefix : Prefix of the paths to be added
%
%      Output :
%   is_loaded : Flag on whether some path was not loaded appropriately
%
% Function layout :
%
% 1. Add paths
%
% 2. Check if there was any warning message when adding the paths
%
%% 1. Add paths

path_prefix = char(path_prefix);

% Add general math functions
addpath([path_prefix 'generalMath/']);

% Add general auxiliary functions
addpath([path_prefix 'auxiliary/']);

% Add system solvers
addpath([path_prefix 'equationSystemSolvers/']);

% Add efficient computation functions
addpath([path_prefix 'efficientComputation/']);

% Add transient analysis solvers
addpath([path_prefix 'transientAnalysis/']);

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath([path_prefix 'CAGDKernel/CAGDKernel_basisFunctions'],...
        [path_prefix 'CAGDKernel/CAGDKernel_geometryResolutionRefinement/'],...
        [path_prefix 'CAGDKernel/CAGDKernel_baseVectors/'],...
        [path_prefix 'CAGDKernel/CAGDKernel_graphics/'],...
        [path_prefix 'CAGDKernel/CAGDKernel_BSplineCurve/'],...
        [path_prefix 'CAGDKernel/CAGDKernel_BSplineSurface/']);
    
% Add all functions related to the isogeometric Kirchhoff-Love shell formulation
addpath([path_prefix 'isogeometricThinStructureAnalysis/graphicsSinglePatch/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/graphicsMultipatches/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/loads/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/solutionMatricesAndVectors/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/solvers/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/metrics/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/auxiliary/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/postprocessing/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/BOperatorMatrices/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/penaltyDecompositionKLShell/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/penaltyDecompositionMembrane/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/lagrangeMultipliersDecompositionKLShell/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/lagrangeMultipliersDecompositionMembrane/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/nitscheDecompositionMembrane/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/errorComputation/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/output/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/transientAnalysis/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/initialConditions/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/weakDBCMembrane/'],...
        [path_prefix 'isogeometricThinStructureAnalysis/formFindingAnalysis/']);

%% 2. Check if there was any warning message when adding the paths
is_loaded = true;
[warnMsg, ~] = lastwarn;
if ~isempty(warnMsg)
    is_loaded = false;
end

end