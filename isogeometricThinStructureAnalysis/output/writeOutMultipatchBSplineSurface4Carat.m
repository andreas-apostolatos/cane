function writeOutMultipatchBSplineSurface4Carat ...
    (BSplinePatches, strongDBC, connections, pathToOutput, caseName)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Writes out an csv file with the multipatch B-Spline information in a
% Carat++ input file.
%
%          Input :
% BSplinePatches : Array of BSpline patches each of which contains,
%                        .p,.q : Polynomial orders in xi-,eta- directions
%                     .Xi,.Eta : Knot vectors in xi-,eta- directions
%                          .CP : Control Points in xi-,eta- directions
%                     .isNURBS : Boolean on whether the B-Spline is a NURBS
%                                or not
%      strongDBC : Structure array containing for each patch the extension
%                  of its boundary where strongly imposed conditions are
%                  applied
%    connections : Define the connection between the patches:
%                        .No : Number of connections
%                 .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                 ...      ...    ...   ...  ...   ...]
%   pathToOutput : The directory on which the results are written out
%       caseName : The name of the case which is assigned to the file name
%
%         Output :
%                  No output but writting the results out into a file with
%                  the name extension caseName.csv under the directory
%                  pathToOutput
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over the patches
% ->
%    1i. Get the patch data
%
%   1ii. Write out the data into file
% <-
%
%% Function main body

%% 0. Read input

% No patches
noPatches = length(BSplinePatches);

% Number of coodinates of the Control points
noCoord = 4;

% Initialize counters
counterCPs = 1;
counterWeakDBC = 1;
counterStrongDBC = 1;

% Total number of weak Dirichlet boundary conditions
noWeakDBC = 0;
for iPatches = 1:noPatches
    if ~isempty(BSplinePatches{iPatches}.weakDBC)
        noWeakDBC = noWeakDBC + BSplinePatches{iPatches}.weakDBC.noCnd;
    end
end
idWeakDBC = zeros(noWeakDBC,2);

% Array containing the global ID of the connections
if isfield(connections,'No')
    if connections.No > 1
        idConn = zeros(connections.No,1);
    end
end

% Make directory to write out the results of the analysis
isExistent = exist(strcat(pathToOutput,caseName),'dir');
if ~isExistent
    mkdir(strcat(pathToOutput,caseName));
end
fileHandle = fopen(strcat(pathToOutput,caseName,'/',caseName,'_GiD.georhino.txt'),'w');

% Material properties
parameters = BSplinePatches{1}.parameters;

%% 1. Write the analysis information
fprintf(fileHandle,'!###################################################################\n');
fprintf(fileHandle,'!####                          PC-BLOCK                         ####\n');
fprintf(fileHandle,'!###################################################################\n');
fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-PROBLEM\n');
fprintf(fileHandle,'  MASTERJOB = PC-ANALYSIS 1\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-ANALYSIS 1: EMPIRE_CoSimulation\n');
fprintf(fileHandle,'  CARAT_ANALYSIS = PC-ANALYSIS 3\n');
fprintf(fileHandle,'  COSIMULATION_INTERFACE = DESIGN_ELEMENT2D ');
for iPatches = 1:noPatches
    if iPatches < noPatches
        fprintf(fileHandle,[num2str(iPatches) blanks(1)]);
    else
        fprintf(fileHandle,[num2str(iPatches) '\n']);
    end
end
fprintf(fileHandle,'  EMPIRE_INPUT_FILE = empireCarat.xml\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-ANALYSIS 2: STA_GEO_NONLIN\n');
fprintf(fileHandle,'  PATHCONTROL = ARCLENGTH      !Options: FORCE or DISPLACEMENT or ARCLENGTH\n');
fprintf(fileHandle,'  SOLVER = PC-SOLVER 1\n');
fprintf(fileHandle,'  OUTPUT = PC-OUT 1\n');
fprintf(fileHandle,'  COMPCASE = LD-COM 1\n');
fprintf(fileHandle,'  DOMAIN = EL-DOMAIN 1\n');
fprintf(fileHandle,'  NUM_STEP = 1\n');
fprintf(fileHandle,'  MAX_ITER_EQUILIBRIUM = 100\n');
fprintf(fileHandle,'  EQUILIBRIUM_ACCURACY = 1e-8\n');
fprintf(fileHandle,'  CURVE = LD-CURVE 1\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-ANALYSIS 4: DYNAMIC\n');
fprintf(fileHandle,'  SOLVER = PC-SOLVER 1\n');
fprintf(fileHandle,'  STARTTIME = 0.0\n');
fprintf(fileHandle,'  ENDTIME   = 30.0\n');
fprintf(fileHandle,'  TIMESTEP  = 0.01\n');
fprintf(fileHandle,'  ALGORITHM = GENALPHA_NLN  !NEWMARK_NLN !GENALPHA_LIN,NEWMARK_LIN,CENTRALDIFFERENCES_LIN,NEWMARK_NLN,GENALPHA_NLN\n');
fprintf(fileHandle,'  BETA     = 0.309\n');
fprintf(fileHandle,'  GAMMA    = 0.611\n');
fprintf(fileHandle,'  ALPHA_M  = 0.333\n');
fprintf(fileHandle,'  ALPHA_F  = 0.444\n');
fprintf(fileHandle,'  OUTPUT   = PC-OUT 1\n');
fprintf(fileHandle,'  COMPCASE = LD-COM 1\n');
fprintf(fileHandle,'  DOMAIN   = EL-DOMAIN 1\n');
fprintf(fileHandle,'  MAX_ITER_EQUILIBRIUM = 50\n');
fprintf(fileHandle,'  EQUILIBRIUM_ACCURACY = 1e-7\n');
fprintf(fileHandle,'  DAMPING = 0         ! 0=off, 1=on\n');
fprintf(fileHandle,'  A1 = 2.709094\n');
fprintf(fileHandle,'  A2 = 0.000441\n');
fprintf(fileHandle,'  RESTARTRUN = 0           ! 0=off, 1=on\n');
fprintf(fileHandle,'  RESTARTOUTPUT = 0        ! 0=off, 1=on\n');
fprintf(fileHandle,'  RESTARTFREQUENCY = 10\n');
fprintf(fileHandle,'  RESTARTFILEPREFIX = restartfile\n');
fprintf(fileHandle,'  RESTARTINFOINSTANCES = 3\n');
fprintf(fileHandle,'  INITIAL_STATIC_CALCULATION = TRUE\n');
fprintf(fileHandle,'  INITIAL_LOAD_STEPS=1\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-SOLVER 1: CROUT_SKYLINE\n');
fprintf(fileHandle,'  BANDWITH = CUTHILL_MCKEE\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-OUT 1 : RHINO\n');
fprintf(fileHandle,'   GEOM=1\n');
fprintf(fileHandle,'   DISP=1\n');
fprintf(fileHandle,'   !STRESS=1\n');
fprintf(fileHandle,'   PREC=7\n');
fprintf(fileHandle,'   FPN=1\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'LD-CURVE 1 TYPE=DISCRETE\n');
fprintf(fileHandle,'  TIME = 0  VAL = 0\n');
fprintf(fileHandle,'  TIME = 30  VAL = 1\n');

fprintf(fileHandle,'!###################################################################\n');
fprintf(fileHandle,'!####                          EL-BLOCK                         ####\n');
fprintf(fileHandle,'!###################################################################\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'EL-PART 1\n');
fprintf(fileHandle,'  NAME=Support\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'EL-MAT 1 : LIN_ELAST_ISOTROPIC \n');
fprintf(fileHandle,['  EMOD = ' num2str(parameters.E) ' \n']);
fprintf(fileHandle,'  ALPHAT =  0.0\n');
fprintf(fileHandle,['DENS = ' num2str(parameters.rho) '\n']);
fprintf(fileHandle,['NUE = ' num2str(parameters.nue) '\n']);

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'EL-PROP 1 : MEMBRANE_NURBS\n');
fprintf(fileHandle,'  MAT= EL-MAT 1\n');
fprintf(fileHandle,['  THICKNESS =  ' num2str(parameters.t) ' \n']);
fprintf(fileHandle,'  INT_TYPE_MEMBRANE_NURBS = FULL\n');
if isfield(parameters,'prestress')
    if ~isa(parameters.prestress.voigtVector,'function_handle')
        fprintf(fileHandle,['  PROJECTED_PRESTRESS=1 SIG11 = ' num2str(parameters.prestress.voigtVector(1,1)) '     SIG22 = ' num2str(parameters.prestress.voigtVector(2,1)) '      SIG12 = ' num2str(parameters.prestress.voigtVector(3,1)) '\n']);
    end
end
fprintf(fileHandle,'  A_X= 1     A_Y= 0    A_Z= 0   !AREA FOR THE DEFINITION OF THE PRESTRESS\n');
fprintf(fileHandle,'  B_X= 0     B_Y= 1    B_Z= 0   !VECTOR A AND VECTOR B DEFINE THE AREA\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'EL-PROP 2 : NURBS_BREP_EDGE_COUPLING\n');
fprintf(fileHandle,'  BREP_TYPE = COUPLING_DISP\n');
fprintf(fileHandle,'  BL_FACTOR_DISP = 1e7\n');
fprintf(fileHandle,'  INT_TYPE_NURBS_BREP_EDGE_COUPLING = FULL\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'EL-PROP 3 : NURBS_BREP_EDGE_COUPLING\n');
fprintf(fileHandle,'  BREP_TYPE = COUPLING_DISP_ROT\n');
fprintf(fileHandle,'  BL_FACTOR_DISP = 1e7\n');
fprintf(fileHandle,'  BL_FACTOR_ROT = 1e7\n');
fprintf(fileHandle,'  INT_TYPE_NURBS_BREP_EDGE_COUPLING = FULL\n');

%% 2. Loop over all patches and write the geometrical information of each patch out
for iPatches = 1:noPatches
    %% 2i. Get the patch data
    p = BSplinePatches{iPatches}.p;
    q = BSplinePatches{iPatches}.q;
    Xi = BSplinePatches{iPatches}.Xi;
    Eta = BSplinePatches{iPatches}.Eta;
    CP = BSplinePatches{iPatches}.CP;
    mxi = length(Xi);
    meta = length(Eta);
    nxi = length(CP(:,1,1));
    neta = length(CP(1,:,1));
    
    %% 2ii. Write out the block with the NURBS data
    fprintf(fileHandle,'!###################################################################\n');
    fprintf(fileHandle,'!####                        NURBS-BLOCK                        ####\n');
    fprintf(fileHandle,'!###################################################################\n');
    fprintf(fileHandle,strcat('NURBS_PATCH',char([blanks(2) num2str(iPatches)]),' : NURBS_2D\n'));
    fprintf(fileHandle,strcat(char([blanks(2) 'CTRL_PTS' ' = ' 'CTRL_PTS_DEF' blanks(2) num2str(iPatches)]),'\n'));
    fprintf(fileHandle,strcat(char([blanks(2) 'NCTRL' ' =' blanks(2) num2str(nxi - 1)]),'\n'));
    fprintf(fileHandle,strcat(char([blanks(2) 'MCTRL' ' =' blanks(2) num2str(neta - 1)]),'\n'));
    fprintf(fileHandle,strcat(char([blanks(2) 'PDEG' blanks(1) '=' blanks(3) num2str(p)]),'\n'));
    fprintf(fileHandle,strcat(char([blanks(2) 'QDEG' blanks(1) '=' blanks(3) num2str(q)]),'\n'));
    fprintf(fileHandle,char([blanks(2) 'UKNOT' blanks(1) '=' blanks(2)]));
    for iXi = 1:mxi
        if iXi < mxi
            fprintf(fileHandle,'%.16f, ',Xi(iXi));
        else
            fprintf(fileHandle,'%.16f\n',Xi(iXi));
        end
    end
    fprintf(fileHandle,char([blanks(2) 'VKNOT' blanks(1) '=' blanks(2)]));
    for iEta = 1:meta
        if iEta < meta
            fprintf(fileHandle,'%.16f, ',Eta(iEta));
        else
            fprintf(fileHandle,'%.16f\n',Eta(iEta));
        end
    end
    fprintf(fileHandle,strcat(char([blanks(2) 'TRIMMING' ' = ' 'B_REP ' num2str(iPatches)]),'\n'));
    fprintf(fileHandle,'!===================================================================\n');
    fprintf(fileHandle,strcat('CTRL_PTS_DEF',char([blanks(2) num2str(iPatches)]),'\n'));
    counter = 1;
    for iEta = 1:neta
        for iXi = 1:nxi
            fprintf(fileHandle,char([blanks(1) 'CTRLPT'  blanks(2) num2str(counter) blanks(2)]));
            for iCoord = 1:noCoord
                if iCoord < noCoord
                    fprintf(fileHandle,'%.16f  ',CP(iXi,iEta,iCoord));
                else
                    fprintf(fileHandle,'%.16f\n',CP(iXi,iEta,iCoord));
                end
            end
                counter = counter + 1;
        end
    end
    
    %% 2ii. Write out all boundary edges of the patch in its parameter space
    verticesMap = [Xi(1)    Eta(1)   0 1
                   Xi(end)  Eta(1)   0 1
                   Xi(end)  Eta(end) 0 1
                   Xi(1)    Eta(end) 0 1];
               
	%% 2iii. Loop over all the boundary edges of the patch
    for iBEdges = 1:length(verticesMap(:,1))
        %% 2iii.1. Write the preamble out
        fprintf(fileHandle,'!###################################################################\n');
        fprintf(fileHandle,'!####                  NURBS-BLOCK-PARAMETER                    ####\n');
        fprintf(fileHandle,'!###################################################################\n');
        
        %% 2iii.2. Write the geometrical information of the curve out
        strBRepEdge = [num2str(iPatches) '00' num2str(iBEdges)];
        fprintf(fileHandle,['NURBS_PATCH_PAR' blanks(2) strBRepEdge ' : NURBS_1D\n']);
        fprintf(fileHandle,['  CTRL_PTS = CTRL_PTS_PAR' blanks(2) strBRepEdge '\n']);
        fprintf(fileHandle,'  NCTRL =  1\n');
        fprintf(fileHandle,'  PDEG  =  1\n');
        fprintf(fileHandle,char([blanks(2) 'UKNOT' blanks(1) '=' blanks(2)]));
        if iBEdges == 1 || iBEdges == 3
            knotVct = [Xi(1) Xi(1) Xi(end) Xi(end)];
        elseif iBEdges == 2 || iBEdges == 4
            knotVct = [Eta(1) Eta(1) Eta(end) Eta(end)];
        end
        for iKnotVct = 1:length(knotVct)
            if iKnotVct < length(knotVct)
                fprintf(fileHandle,'%.16f, ',knotVct(iKnotVct));
            else
                fprintf(fileHandle,'%.16f\n',knotVct(iKnotVct));
            end
        end
        fprintf(fileHandle,'!===================================================================\n');
        
        %% 2iii.3. Write the corresponding Control Point IDs out
        fprintf(fileHandle,['CTRL_PTS_PAR' blanks(2) strBRepEdge '\n']);
        for iCPs = 1:length(verticesMap(1,1:2))
            fprintf(fileHandle,['CTRLPT_PAR_ID' blanks(1) num2str(counterCPs) '\n']);
            counterCPs = counterCPs + 1;
        end
    end
    
    %% 2iv. Write the trimming information of the patch out
    fprintf(fileHandle,'!###################################################################\n');
    fprintf(fileHandle,'!####                        TRIMMING                           ####\n');
    fprintf(fileHandle,'!###################################################################\n');
    fprintf(fileHandle,['B_REP ' num2str(iPatches) '\n']);
    fprintf(fileHandle,[' B_REP_LOOP' blanks(3)  num2str(iPatches) blanks(2) '0' blanks(2)]);
    for iBEdges = 1:length(verticesMap(:,1))
        strBRepEdge = [num2str(iPatches) '00' num2str(iBEdges)];
        if iBEdges < length(verticesMap(:,1))
            fprintf(fileHandle,[strBRepEdge blanks(2) 'TRUE' blanks(2)]);
        else
            fprintf(fileHandle,[strBRepEdge blanks(2) 'TRUE\n']);
        end
    end
    
    %% 2v. Define the Control Points of the boundary curves in the parameter space of the patch
    fprintf(fileHandle,'!==================================================================\n');
    fprintf(fileHandle,'CTRL_PTS_PAR_DEF\n');
	counterIndexCPs = counterCPs - 2*length(verticesMap(:,1)) : counterCPs - 1;
    counter = 1;
    for iCPs = 1:length(verticesMap(:,1))
        if iCPs < length(verticesMap(:,1))
            range = [iCPs iCPs + 1];
        else
            range = [iCPs 1];
        end
        for iCPsInner = range
            indexCP = counterCPs - 2*iPatches*length(verticesMap(:,1)) - 1 + iCPsInner;
            fprintf(fileHandle,[blanks(2) 'CTRLPT_PAR' blanks(2) num2str(counterIndexCPs(counter)) blanks(2)]);
            counter = counter + 1;
            for iCoord = 1:noCoord
                if iCoord < noCoord
                    fprintf(fileHandle,'%.16f  ',verticesMap(indexCP,iCoord));
                else
                    fprintf(fileHandle,'%.16f\n',verticesMap(indexCP,iCoord));
                end
            end
        end
    end
end

%% 3. Loop ove all patches and write the design block information
fprintf(fileHandle,'!###################################################################\n');
fprintf(fileHandle,'!####                        DESIGN-BLOCK                       ####\n');
fprintf(fileHandle,'!###################################################################\n');
fprintf(fileHandle,'!         ID  PART  PROP   NURBS_TOP\n');
fprintf(fileHandle,'DE-ELTOP-NURBS\n');
for iPatches = 1:noPatches
    fprintf(fileHandle,['DE-EL' blanks(3) num2str(iPatches) blanks(3) '1' blanks(3) '1' blanks(3) num2str(iPatches) '\n']);
end
fprintf(fileHandle,'!==================================================================\n');

%% 4. Write out the elements
fprintf(fileHandle,'EL-DOMAIN 1\n');
fprintf(fileHandle,'  ELEMENTS = EL-TOP 1\n');
fprintf(fileHandle,'!==================================================================\n');

%% 5. Write preamble for the weak application of conditions
fprintf(fileHandle,'!                  ID  PART  PROP   NURBS_TOP\n');
fprintf(fileHandle,'DE-BREP-NURBS\n');

%% 6. Loop over all patches and add the weak Dirichlet boundary conditions
for iPatches = 1:noPatches
    %% 6i. Get the weak Dirichlet boundary conditions of the patch
    weakDBC = BSplinePatches{iPatches}.weakDBC;
    
    %% 6i. Loop over all weak Dirichlet boundary conditions of the patch
    if ~isempty(weakDBC)
        for iWeakDBC = 1:weakDBC.noCnd
            xiExtension = weakDBC.xiExtension{iWeakDBC};
            etaExtension = weakDBC.etaExtension{iWeakDBC};
            if etaExtension(1,1) == etaExtension(1,2)
                if etaExtension(1,1) == 0
                    idCurve =  [num2str(iPatches) '001'];
                elseif etaExtension(1,1) == 1
                    idCurve =  [num2str(iPatches) '003'];
                else
                    error('Eta extension should run from 0 to 1 but not %d',etaExtension(1,1));
                end
            elseif xiExtension(1,1) == xiExtension(1,2)
                if xiExtension(1,1) == 0
                    idCurve =  [num2str(iPatches) '004'];
                elseif xiExtension(1,1) == 1
                    idCurve =  [num2str(iPatches) '002'];
                else
                    error('Xi extension should run from 0 to 1 but not %d',xiExtension(1,1));
                end
            else
                error('Wrong definition of the weak Dirichlet conditions %d at patch %d',iWeakDBC,iPatches);
            end
            idWeakDBC(counterWeakDBC,1) = str2double(idCurve);
            idWeakDBC(counterWeakDBC,2) = iPatches;
            fprintf(fileHandle,['  DE-BREP-EL' blanks(3) idCurve blanks(3) '1' blanks(3) '2' blanks(3) 'PATCH_PAR' blanks(2) idCurve blanks(2) 'OF' blanks(3) 'DE-EL' blanks(2) num2str(iPatches) '\n']);
            counterWeakDBC = counterWeakDBC + 1;
        end
    end
end

%% 7. Loop over all the connections and write out the patch connectivity information
if isfield(connections,'No')
    if connections.No > 1
        for iConn = 1:connections.No
            idI = connections.xiEtaCoup(iConn,1);
            idJ = connections.xiEtaCoup(iConn,2);
            XiI = BSplinePatches{idI}.Xi;
            EtaI = BSplinePatches{idI}.Eta;
            XiJ = BSplinePatches{idJ}.Xi;
            EtaJ = BSplinePatches{idJ}.Eta;
            xicoupI = connections.xiEtaCoup(iConn,3:4);
            etacoupI = connections.xiEtaCoup(iConn,5:6);
            xicoupJ = connections.xiEtaCoup(iConn,7:8);
            etacoupJ = connections.xiEtaCoup(iConn,9:10);    
            if etacoupI(1,1) == etacoupI(1,2)
                if etacoupI(1,1) == 0
                    idCurveI =  [num2str(idI) '001'];
                    knotStartI = XiI(1);
                    knotEndI = XiI(end);
                elseif etacoupI(1,1) == 1
                    idCurveI =  [num2str(idI) '003'];
                    knotStartI = XiI(end);
                    knotEndI = XiI(1);
                else
                    error('Eta coupling extension of patch %d should run from 0 to 1 but not %d',idI,etacoupI(1,1));
                end
            elseif xicoupI(1,1) == xicoupI(1,2)
                if xicoupI(1,1) == 0
                    idCurveI =  [num2str(idI) '004'];
                    knotStartI = EtaI(end);
                    knotEndI = EtaI(1);
                elseif xicoupI(1,1) == 1
                    idCurveI =  [num2str(idI) '002'];
                    knotStartI = EtaI(1);
                    knotEndI = EtaI(end);
                else
                    error('Xi coupling extension of patch %d should run from 0 to 1 but not %d',idI,xicoupI(1,1));
                end
            else
                error('Wrong definition of connection %d at patch %d',iConn,idI);
            end
            if etacoupJ(1,1) == etacoupJ(1,2)
                if etacoupJ(1,1) == 0
                    idCurveJ =  [num2str(idJ) '001'];
                    knotStartJ = XiJ(1);
                    knotEndJ = XiJ(end);
                elseif etacoupJ(1,1) == 1
                    idCurveJ =  [num2str(idJ) '003'];
                    knotStartJ = XiJ(end);
                    knotEndJ = XiJ(1);
                else
                    error('Eta coupling extension of patch %d should run from 0 to 1 but not %d',idJ,etacoupJ(1,1));
                end
            elseif xicoupJ(1,1) == xicoupJ(1,2)
                if xicoupJ(1,1) == 0
                    idCurveJ =  [num2str(idJ) '004'];
                    knotStartJ = EtaJ(end);
                    knotEndJ = EtaJ(1);
                elseif xicoupJ(1,1) == 1
                    idCurveJ =  [num2str(idJ) '002'];
                    knotStartJ = EtaJ(1);
                    knotEndJ = EtaJ(end);
                else
                    error('Xi coupling extension of patch %d should run from 0 to 1 but not %d',idJ,xicoupJ(1,1));
                end
            else
                error('Wrong definition of connection %d at patch %d',iConn,idJ);
            end
            idConn(iConn,1) = str2double(idCurveI);
            fprintf(fileHandle,['  DE-BREP-EL' blanks(3) idCurveI blanks(3) '1' blanks(3) '3' blanks(3) 'PATCH_PAR' blanks(2) idCurveI blanks(2) 'OF' blanks(3) 'DE-EL' blanks(2) num2str(idI) blanks(1)...
                'PATCH_PAR' blanks(2) idCurveJ blanks(2) 'OF' blanks(3) 'DE-EL' blanks(2) num2str(idJ) blanks(1) 'TOL=0.01' blanks(1) ...
                '[' num2str(knotStartI) blanks(1) num2str(knotEndI) ']' '[' num2str(knotStartJ) blanks(1) num2str(knotEndJ) ']' '\n']);
        end
    end
end

%% 8. Write out the preamble for the refinement
fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'DE-REFINEMENT\n');

%% 9. Define the refinement for each patch
for iPatches = 1:noPatches
    fprintf(fileHandle,[' DE-EL' blanks(3) num2str(iPatches) blanks(3) 'ep=0' blanks(3) 'eq=0' blanks(3) 'ru=0' blanks(3) 'rv=0\n']);
end

%% 10. Define the refinement for each weak Dirichlet condition
for iWeakDBC = 1:noWeakDBC
    fprintf(fileHandle,[' DE-BREP-EL' blanks(3) num2str(idWeakDBC(iWeakDBC,1)) blanks(3) 'dp=auto' blanks(3) 'ru=auto\n']);
end

%% 11. Define the refinement for each b-rep element associated to multipatch coupling
if isfield(connections,'No')
    if connections.No > 1
        for iConn = 1:connections.No
            fprintf(fileHandle,[' DE-BREP-EL' blanks(3) num2str(idConn(iConn,1)) blanks(3) 'dp=auto' blanks(3) 'ru=auto\n']);
        end
    end
end

%% 12. Preamble for the support conditions
fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'!        ID  DE-EL     LOC COORD  BC\n');

%% 13. Loop over all weak Dirichlet boundary conditions
for iWeakDBC = 1:noWeakDBC
    idCurve = idWeakDBC(iWeakDBC,1);
    idPatch = idWeakDBC(iWeakDBC,2);
    fprintf(fileHandle,['DE-SUP' blanks(3) num2str(iWeakDBC) blanks(3) num2str(idPatch) blanks(11) 'DE-BREP' blanks(2) num2str(idCurve) blanks(5) 'DISP_X,DISP_Y,DISP_Z !CLAMPED' '\n']);
end

%% 14. Loop over all the strong Dirichlet boundary conditions
if ~isempty(strongDBC)
    for iPatches = 1:noPatches
        strongdbc = strongDBC{iPatches};
        if ~isempty(strongdbc)
            for iStrongDBC = 1:length(strongdbc.xiExtensionHom)
                xiExtension = strongdbc.xiExtensionHom{iStrongDBC};
                etaExtension = strongdbc.etaExtensionHom{iStrongDBC};
                if etaExtension(1,1) == etaExtension(1,2) && xiExtension(1,1) == xiExtension(1,2)
                    xiKnot = xiExtension(1,1);function writeOutMultipatchBSplineSurface4Carat ...
    (BSplinePatches, strongDBC, connections, pathToOutput, caseName)
%% Licensing
%
% Copyright (c)  2020, by Technische Universitaet Muenchen
%                Dr.-Ing. Andreas Apostolatos and Dr.-Ing. Roland Wuechner
%                Lehrstuhl fuer Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger
%
% All rights reserved.
%
% AFEM default license: afem/license.txt
%
%% Function documentation
%
% Writes out an csv file with the multipatch B-Spline information in a
% Carat++ input file.
%
%          Input :
% BSplinePatches : Array of BSpline patches each of which contains,
%                        .p,.q : Polynomial orders in xi-,eta- directions
%                     .Xi,.Eta : Knot vectors in xi-,eta- directions
%                          .CP : Control Points in xi-,eta- directions
%                     .isNURBS : Boolean on whether the B-Spline is a NURBS
%                                or not
%      strongDBC : Structure array containing for each patch the extension
%                  of its boundary where strongly imposed conditions are
%                  applied
%    connections : Define the connection between the patches:
%                        .No : Number of connections
%                 .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                 ...      ...    ...   ...  ...   ...]
%   pathToOutput : The directory on which the results are written out
%       caseName : The name of the case which is assigned to the file name
%
%         Output :
%                  No output but writting the results out into a file with
%                  the name extension caseName.csv under the directory
%                  pathToOutput
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over the patches
% ->
%    1i. Get the patch data
%
%   1ii. Write out the data into file
% <-
%
%% Function main body

%% 0. Read input

% No patches
noPatches = length(BSplinePatches);

% Number of coodinates of the Control points
noCoord = 4;

% Initialize counters
counterCPs = 1;
counterWeakDBC = 1;
counterStrongDBC = 1;

% Total number of weak Dirichlet boundary conditions
noWeakDBC = 0;
for iPatches = 1:noPatches
    if ~isempty(BSplinePatches{iPatches}.weakDBC)
        noWeakDBC = noWeakDBC + BSplinePatches{iPatches}.weakDBC.noCnd;
    end
end
idWeakDBC = zeros(noWeakDBC,2);

% Array containing the global ID of the connections
if isfield(connections,'No')
    if connections.No > 1
        idConn = zeros(connections.No,1);
    end
end

% Make directory to write out the results of the analysis
isExistent = exist(strcat(pathToOutput,caseName),'dir');
if ~isExistent
    mkdir(strcat(pathToOutput,caseName));
end
fileHandle = fopen(strcat(pathToOutput,caseName,'/',caseName,'_GiD.georhino.txt'),'w');

% Material properties
parameters = BSplinePatches{1}.parameters;

%% 1. Write the analysis information
fprintf(fileHandle,'!###################################################################\n');
fprintf(fileHandle,'!####                          PC-BLOCK                         ####\n');
fprintf(fileHandle,'!###################################################################\n');
fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-PROBLEM\n');
fprintf(fileHandle,'  MASTERJOB = PC-ANALYSIS 1\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-ANALYSIS 1: EMPIRE_CoSimulation\n');
fprintf(fileHandle,'  CARAT_ANALYSIS = PC-ANALYSIS 3\n');
fprintf(fileHandle,'  COSIMULATION_INTERFACE = DESIGN_ELEMENT2D ');
for iPatches = 1:noPatches
    if iPatches < noPatches
        fprintf(fileHandle,[num2str(iPatches) blanks(1)]);
    else
        fprintf(fileHandle,[num2str(iPatches) '\n']);
    end
end
fprintf(fileHandle,'  EMPIRE_INPUT_FILE = empireCarat.xml\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-ANALYSIS 2: STA_GEO_NONLIN\n');
fprintf(fileHandle,'  PATHCONTROL = ARCLENGTH      !Options: FORCE or DISPLACEMENT or ARCLENGTH\n');
fprintf(fileHandle,'  SOLVER = PC-SOLVER 1\n');
fprintf(fileHandle,'  OUTPUT = PC-OUT 1\n');
fprintf(fileHandle,'  COMPCASE = LD-COM 1\n');
fprintf(fileHandle,'  DOMAIN = EL-DOMAIN 1\n');
fprintf(fileHandle,'  NUM_STEP = 1\n');
fprintf(fileHandle,'  MAX_ITER_EQUILIBRIUM = 100\n');
fprintf(fileHandle,'  EQUILIBRIUM_ACCURACY = 1e-8\n');
fprintf(fileHandle,'  CURVE = LD-CURVE 1\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-ANALYSIS 4: DYNAMIC\n');
fprintf(fileHandle,'  SOLVER = PC-SOLVER 1\n');
fprintf(fileHandle,'  STARTTIME = 0.0\n');
fprintf(fileHandle,'  ENDTIME   = 30.0\n');
fprintf(fileHandle,'  TIMESTEP  = 0.01\n');
fprintf(fileHandle,'  ALGORITHM = GENALPHA_NLN  !NEWMARK_NLN !GENALPHA_LIN,NEWMARK_LIN,CENTRALDIFFERENCES_LIN,NEWMARK_NLN,GENALPHA_NLN\n');
fprintf(fileHandle,'  BETA     = 0.309\n');
fprintf(fileHandle,'  GAMMA    = 0.611\n');
fprintf(fileHandle,'  ALPHA_M  = 0.333\n');
fprintf(fileHandle,'  ALPHA_F  = 0.444\n');
fprintf(fileHandle,'  OUTPUT   = PC-OUT 1\n');
fprintf(fileHandle,'  COMPCASE = LD-COM 1\n');
fprintf(fileHandle,'  DOMAIN   = EL-DOMAIN 1\n');
fprintf(fileHandle,'  MAX_ITER_EQUILIBRIUM = 50\n');
fprintf(fileHandle,'  EQUILIBRIUM_ACCURACY = 1e-7\n');
fprintf(fileHandle,'  DAMPING = 0         ! 0=off, 1=on\n');
fprintf(fileHandle,'  A1 = 2.709094\n');
fprintf(fileHandle,'  A2 = 0.000441\n');
fprintf(fileHandle,'  RESTARTRUN = 0           ! 0=off, 1=on\n');
fprintf(fileHandle,'  RESTARTOUTPUT = 0        ! 0=off, 1=on\n');
fprintf(fileHandle,'  RESTARTFREQUENCY = 10\n');
fprintf(fileHandle,'  RESTARTFILEPREFIX = restartfile\n');
fprintf(fileHandle,'  RESTARTINFOINSTANCES = 3\n');
fprintf(fileHandle,'  INITIAL_STATIC_CALCULATION = TRUE\n');
fprintf(fileHandle,'  INITIAL_LOAD_STEPS=1\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-SOLVER 1: CROUT_SKYLINE\n');
fprintf(fileHandle,'  BANDWITH = CUTHILL_MCKEE\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-OUT 1 : RHINO\n');
fprintf(fileHandle,'   GEOM=1\n');
fprintf(fileHandle,'   DISP=1\n');
fprintf(fileHandle,'   !STRESS=1\n');
fprintf(fileHandle,'   PREC=7\n');
fprintf(fileHandle,'   FPN=1\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'LD-CURVE 1 TYPE=DISCRETE\n');
fprintf(fileHandle,'  TIME = 0  VAL = 0\n');
fprintf(fileHandle,'  TIME = 30  VAL = 1\n');

fprintf(fileHandle,'!###################################################################\n');
fprintf(fileHandle,'!####                          EL-BLOCK                         ####\n');
fprintf(fileHandle,'!###################################################################\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'EL-PART 1\n');
fprintf(fileHandle,'  NAME=Support\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'EL-MAT 1 : LIN_ELAST_ISOTROPIC \n');
fprintf(fileHandle,['  EMOD = ' num2str(parameters.E) ' \n']);
fprintf(fileHandle,'  ALPHAT =  0.0\n');
fprintf(fileHandle,['DENS = ' num2str(parameters.rho) '\n']);
fprintf(fileHandle,['NUE = ' num2str(parameters.nue) '\n']);

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'EL-PROP 1 : MEMBRANE_NURBS\n');
fprintf(fileHandle,'  MAT= EL-MAT 1\n');
fprintf(fileHandle,['  THICKNESS =  ' num2str(parameters.t) ' \n']);
fprintf(fileHandle,'  INT_TYPE_MEMBRANE_NURBS = FULL\n');
if isfield(parameters,'prestress')
    if ~isa(parameters.prestress.voigtVector,'function_handle')
        fprintf(fileHandle,['  PROJECTED_PRESTRESS=1 SIG11 = ' num2str(parameters.prestress.voigtVector(1,1)) '     SIG22 = ' num2str(parameters.prestress.voigtVector(2,1)) '      SIG12 = ' num2str(parameters.prestress.voigtVector(3,1)) '\n']);
    end
end
fprintf(fileHandle,'  A_X= 1     A_Y= 0    A_Z= 0   !AREA FOR THE DEFINITION OF THE PRESTRESS\n');
fprintf(fileHandle,'  B_X= 0     B_Y= 1    B_Z= 0   !VECTOR A AND VECTOR B DEFINE THE AREA\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'EL-PROP 2 : NURBS_BREP_EDGE_COUPLING\n');
fprintf(fileHandle,'  BREP_TYPE = COUPLING_DISP\n');
fprintf(fileHandle,'  BL_FACTOR_DISP = 1e7\n');
fprintf(fileHandle,'  INT_TYPE_NURBS_BREP_EDGEfunction writeOutMultipatchBSplineSurface4Carat ...
    (BSplinePatches, strongDBC, connections, pathToOutput, caseName)
%% Licensing
%
% Copyright (c)  2020, by Technische Universitaet Muenchen
%                Dr.-Ing. Andreas Apostolatos and Dr.-Ing. Roland Wuechner
%                Lehrstuhl fuer Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger
%
% All rights reserved.
%
% AFEM default license: afem/license.txt
%
%% Function documentation
%
% Writes out an csv file with the multipatch B-Spline information in a
% Carat++ input file.
%
%          Input :
% BSplinePatches : Array of BSpline patches each of which contains,
%                        .p,.q : Polynomial orders in xi-,eta- directions
%                     .Xi,.Eta : Knot vectors in xi-,eta- directions
%                          .CP : Control Points in xi-,eta- directions
%                     .isNURBS : Boolean on whether the B-Spline is a NURBS
%                                or not
%      strongDBC : Structure array containing for each patch the extension
%                  of its boundary where strongly imposed conditions are
%                  applied
%    connections : Define the connection between the patches:
%                        .No : Number of connections
%                 .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                 ...      ...    ...   ...  ...   ...]
%   pathToOutput : The directory on which the results are written out
%       caseName : The name of the case which is assigned to the file name
%
%         Output :
%                  No output but writting the results out into a file with
%                  the name extension caseName.csv under the directory
%                  pathToOutput
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over the patches
% ->
%    1i. Get the patch data
%
%   1ii. Write out the data into file
% <-
%
%% Function main body

%% 0. Read input

% No patches
noPatches = length(BSplinePatches);

% Number of coodinates of the Control points
noCoord = 4;

% Initialize counters
counterCPs = 1;
counterWeakDBC = 1;
counterStrongDBC = 1;

% Total number of weak Dirichlet boundary conditions
noWeakDBC = 0;
for iPatches = 1:noPatches
    if ~isempty(BSplinePatches{iPatches}.weakDBC)
        noWeakDBC = noWeakDBC + BSplinePatches{iPatches}.weakDBC.noCnd;
    end
end
idWeakDBC = zeros(noWeakDBC,2);

% Array containing the global ID of the connections
if isfield(connections,'No')
    if connections.No > 1
        idConn = zeros(connections.No,1);
    end
end

% Make directory to write out the results of the analysis
isExistent = exist(strcat(pathToOutput,caseName),'dir');
if ~isExistent
    mkdir(strcat(pathToOutput,caseName));
end
fileHandle = fopen(strcat(pathToOutput,caseName,'/',caseName,'_GiD.georhino.txt'),'w');

% Material properties
parameters = BSplinePatches{1}.parameters;

%% 1. Write the analysis information
fprintf(fileHandle,'!###################################################################\n');
fprintf(fileHandle,'!####                          PC-BLOCK                         ####\n');
fprintf(fileHandle,'!###################################################################\n');
fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-PROBLEM\n');
fprintf(fileHandle,'  MASTERJOB = PC-ANALYSIS 1\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-ANALYSIS 1: EMPIRE_CoSimulation\n');
fprintf(fileHandle,'  CARAT_ANALYSIS = PC-ANALYSIS 3\n');
fprintf(fileHandle,'  COSIMULATION_INTERFACE = DESIGN_ELEMENT2D ');
for iPatches = 1:noPatches
    if iPatches < noPatches
        fprintf(fileHandle,[num2str(iPatches) blanks(1)]);
    else
        fprintf(fileHandle,[num2str(iPatches) '\n']);
    end
end
fprintf(fileHandle,'  EMPIRE_INPUT_FILE = empireCarat.xml\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-ANALYSIS 2: STA_GEO_NONLIN\n');
fprintf(fileHandle,'  PATHCONTROL = ARCLENGTH      !Options: FORCE or DISPLACEMENT or ARCLENGTH\n');
fprintf(fileHandle,'  SOLVER = PC-SOLVER 1\n');
fprintf(fileHandle,'  OUTPUT = PC-OUT 1\n');
fprintf(fileHandle,'  COMPCASE = LD-COM 1\n');
fprintf(fileHandle,'  DOMAIN = EL-DOMAIN 1\n');
fprintf(fileHandle,'  NUM_STEP = 1\n');
fprintf(fileHandle,'  MAX_ITER_EQUILIBRIUM = 100\n');
fprintf(fileHandle,'  EQUILIBRIUM_ACCURACY = 1e-8\n');
fprintf(fileHandle,'  CURVE = LD-CURVE 1\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-ANALYSIS 4: DYNAMIC\n');
fprintf(fileHandle,'  SOLVER = PC-SOLVER 1\n');
fprintf(fileHandle,'  STARTTIME = 0.0\n');
fprintf(fileHandle,'  ENDTIME   = 30.0\n');
fprintf(fileHandle,'  TIMESTEP  = 0.01\n');
fprintf(fileHandle,'  ALGORITHM = GENALPHA_NLN  !NEWMARK_NLN !GENALPHA_LIN,NEWMARK_LIN,CENTRALDIFFERENCES_LIN,NEWMARK_NLN,GENALPHA_NLN\n');
fprintf(fileHandle,'  BETA     = 0.309\n');
fprintf(fileHandle,'  GAMMA    = 0.611\n');
fprintf(fileHandle,'  ALPHA_M  = 0.333\n');
fprintf(fileHandle,'  ALPHA_F  = 0.444\n');
fprintf(fileHandle,'  OUTPUT   = PC-OUT 1\n');
fprintf(fileHandle,'  COMPCASE = LD-COM 1\n');
fprintf(fileHandle,'  DOMAIN   = EL-DOMAIN 1\n');
fprintf(fileHandle,'  MAX_ITER_EQUILIBRIUM = 50\n');
fprintf(fileHandle,'  EQUILIBRIUM_ACCURACY = 1e-7\n');
fprintf(fileHandle,'  DAMPING = 0         ! 0=off, 1=on\n');
fprintf(fileHandle,'  A1 = 2.709094\n');
fprintf(fileHandle,'  A2 = 0.000441\n');
fprintf(fileHandle,'  RESTARTRUN = 0           ! 0=off, 1=on\n');
fprintf(fileHandle,'  RESTARTOUTPUT = 0        ! 0=off, 1=on\n');
fprintf(fileHandle,'  RESTARTFREQUENCY = 10\n');
fprintf(fileHandle,'  RESTARTFILEPREFIX = restartfile\n');
fprintf(fileHandle,'  RESTARTINFOINSTANCES = 3\n');
fprintf(fileHandle,'  INITIAL_STATIC_CALCULATION = TRUE\n');
fprintf(fileHandle,'  INITIAL_LOAD_STEPS=1\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-SOLVER 1: CROUT_SKYLINE\n');
fprintf(fileHandle,'  BANDWITH = CUTHILL_MCKEE\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'PC-OUT 1 : RHINO\n');
fprintf(fileHandle,'   GEOM=1\n');
fprintf(fileHandle,'   DISP=1\n');
fprintf(fileHandle,'   !STRESS=1\n');
fprintf(fileHandle,'   PREC=7\n');
fprintf(fileHandle,'   FPN=1\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'LD-CURVE 1 TYPE=DISCRETE\n');
fprintf(fileHandle,'  TIME = 0  VAL = 0\n');
fprintf(fileHandle,'  TIME = 30  VAL = 1\n');

fprintf(fileHandle,'!###################################################################\n');
fprintf(fileHandle,'!####                          EL-BLOCK                         ####\n');
fprintf(fileHandle,'!###################################################################\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'EL-PART 1\n');
fprintf(fileHandle,'  NAME=Support\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'EL-MAT 1 : LIN_ELAST_ISOTROPIC \n');
fprintf(fileHandle,['  EMOD = ' num2str(parameters.E) ' \n']);
fprintf(fileHandle,'  ALPHAT =  0.0\n');
fprintf(fileHandle,['DENS = ' num2str(parameters.rho) '\n']);
fprintf(fileHandle,['NUE = ' num2str(parameters.nue) '\n']);

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'EL-PROP 1 : MEMBRANE_NURBS\n');
fprintf(fileHandle,'  MAT= EL-MAT 1\n');
fprintf(fileHandle,['  THICKNESS =  ' num2str(parameters.t) ' \n']);
fprintf(fileHandle,'  INT_TYPE_MEMBRANE_NURBS = FULL\n');
if isfield(parameters,'prestress')
    if ~isa(parameters.prestress.voigtVector,'function_handle')
        fprintf(fileHandle,['  PROJECTED_PRESTRESS=1 SIG11 = ' num2str(parameters.prestress.voigtVector(1,1)) '     SIG22 = ' num2str(parameters.prestress.voigtVector(2,1)) '      SIG12 = ' num2str(parameters.prestress.voigtVector(3,1)) '\n']);
    end
end
fprintf(fileHandle,'  A_X= 1     A_Y= 0    A_Z= 0   !AREA FOR THE DEFINITION OF THE PRESTRESS\n');
fprintf(fileHandle,'  B_X= 0     B_Y= 1    B_Z= 0   !VECTOR A AND VECTOR B DEFINE THE AREA\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'EL-PROP 2 : NURBS_BREP_EDGE_COUPLING\n');
fprintf(fileHandle,'  BREP_TYPE = COUPLING_DISP\n');
fprintf(fileHandle,'  BL_FACTOR_DISP = 1e7\n');
fprintf(fileHandle,'  INT_TYPE_NURBS_BREP_EDGE_COUPLING = FULL\n');

fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'EL-PROP 3 : NURBS_BREP_EDGE_COUPLING\n');
fprintf(fileHandle,'  BREP_TYPE = COUPLING_DISP_ROT\n');
fprintf(fileHandle,'  BL_FACTOR_DISP = 1e7\n');
fprintf(fileHandle,'  BL_FACTOR_ROT = 1e7\n');
fprintf(fileHandle,'  INT_TYPE_NURBS_BREP_EDGE_COUPLING = FULL\n');

%% 2. Loop over all patches and write the geometrical information of each patch out
for iPatches = 1:noPatches
    %% 2i. Get the patch data
    p = BSplinePatches{iPatches}.p;
    q = BSplinePatches{iPatches}.q;
    Xi = BSplinePatches{iPatches}.Xi;
    Eta = BSplinePatches{iPatches}.Eta;
    CP = BSplinePatches{iPatches}.CP;
    mxi = length(Xi);
    meta = length(Eta);
    nxi = length(CP(:,1,1));
    neta = length(CP(1,:,1));
    
    %% 2ii. Write out the block with the NURBS data
    fprintf(fileHandle,'!###################################################################\n');
    fprintf(fileHandle,'!####                        NURBS-BLOCK                        ####\n');
    fprintf(fileHandle,'!###################################################################\n');
    fprintf(fileHandle,strcat('NURBS_PATCH',char([blanks(2) num2str(iPatches)]),' : NURBS_2D\n'));
    fprintf(fileHandle,strcat(char([blanks(2) 'CTRL_PTS' ' = ' 'CTRL_PTS_DEF' blanks(2) num2str(iPatches)]),'\n'));
    fprintf(fileHandle,strcat(char([blanks(2) 'NCTRL' ' =' blanks(2) num2str(nxi - 1)]),'\n'));
    fprintf(fileHandle,strcat(char([blanks(2) 'MCTRL' ' =' blanks(2) num2str(neta - 1)]),'\n'));
    fprintf(fileHandle,strcat(char([blanks(2) 'PDEG' blanks(1) '=' blanks(3) num2str(p)]),'\n'));
    fprintf(fileHandle,strcat(char([blanks(2) 'QDEG' blanks(1) '=' blanks(3) num2str(q)]),'\n'));
    fprintf(fileHandle,char([blanks(2) 'UKNOT' blanks(1) '=' blanks(2)]));
    for iXi = 1:mxi
        if iXi < mxi
            fprintf(fileHandle,'%.16f, ',Xi(iXi));
        else
            fprintf(fileHandle,'%.16f\n',Xi(iXi));
        end
    end
    fprintf(fileHandle,char([blanks(2) 'VKNOT' blanks(1) '=' blanks(2)]));
    for iEta = 1:meta
        if iEta < meta
            fprintf(fileHandle,'%.16f, ',Eta(iEta));
        else
            fprintf(fileHandle,'%.16f\n',Eta(iEta));
        end
    end
    fprintf(fileHandle,strcat(char([blanks(2) 'TRIMMING' ' = ' 'B_REP ' num2str(iPatches)]),'\n'));
    fprintf(fileHandle,'!===================================================================\n');
    fprintf(fileHandle,strcat('CTRL_PTS_DEF',char([blanks(2) num2str(iPatches)]),'\n'));
    counter = 1;
    for iEta = 1:neta
        for iXi = 1:nxi
            fprintf(fileHandle,char([blanks(1) 'CTRLPT'  blanks(2) num2str(counter) blanks(2)]));
            for iCoord = 1:noCoord
                if iCoord < noCoord
                    fprintf(fileHandle,'%.16f  ',CP(iXi,iEta,iCoord));
                else
                    fprintf(fileHandle,'%.16f\n',CP(iXi,iEta,iCoord));
                end
            end
                counter = counter + 1;
        end
    end
    
    %% 2ii. Write out all boundary edges of the patch in its parameter space
    verticesMap = [Xi(1)    Eta(1)   0 1
                   Xi(end)  Eta(1)   0 1
                   Xi(end)  Eta(end) 0 1
                   Xi(1)    Eta(end) 0 1];
               
	%% 2iii. Loop over all the boundary edges of the patch
    for iBEdges = 1:length(verticesMap(:,1))
        %% 2iii.1. Write the preamble out
        fprintf(fileHandle,'!###################################################################\n');
        fprintf(fileHandle,'!####                  NURBS-BLOCK-PARAMETER                    ####\n');
        fprintf(fileHandle,'!###################################################################\n');
        
        %% 2iii.2. Write the geometrical information of the curve out
        strBRepEdge = [num2str(iPatches) '00' num2str(iBEdges)];
        fprintf(fileHandle,['NURBS_PATCH_PAR' blanks(2) strBRepEdge ' : NURBS_1D\n']);
        fprintf(fileHandle,['  CTRL_PTS = CTRL_PTS_PAR' blanks(2) strBRepEdge '\n']);
        fprintf(fileHandle,'  NCTRL =  1\n');
        fprintf(fileHandle,'  PDEG  =  1\n');
        fprintf(fileHandle,char([blanks(2) 'UKNOT' blanks(1) '=' blanks(2)]));
        if iBEdges == 1 || iBEdges == 3
            knotVct = [Xi(1) Xi(1) Xi(end) Xi(end)];
        elseif iBEdges == 2 || iBEdges == 4
            knotVct = [Eta(1) Eta(1) Eta(end) Eta(end)];
        end
        for iKnotVct = 1:length(knotVct)
            if iKnotVct < length(knotVct)
                fprintf(fileHandle,'%.16f, ',knotVct(iKnotVct));
            else
                fprintf(fileHandle,'%.16f\n',knotVct(iKnotVct));
            end
        end
        fprintf(fileHandle,'!===================================================================\n');
        
        %% 2iii.3. Write the corresponding Control Point IDs out
        fprintf(fileHandle,['CTRL_PTS_PAR' blanks(2) strBRepEdge '\n']);
        for iCPs = 1:length(verticesMap(1,1:2))
            fprintf(fileHandle,['CTRLPT_PAR_ID' blanks(1) num2str(counterCPs) '\n']);
            counterCPs = counterCPs + 1;
        end
    end
    
    %% 2iv. Write the trimming information of the patch out
    fprintf(fileHandle,'!###################################################################\n');
    fprintf(fileHandle,'!####                        TRIMMING                           ####\n');
    fprintf(fileHandle,'!###################################################################\n');
    fprintf(fileHandle,['B_REP ' num2str(iPatches) '\n']);
    fprintf(fileHandle,[' B_REP_LOOP' blanks(3)  num2str(iPatches) blanks(2) '0' blanks(2)]);
    for iBEdges = 1:length(verticesMap(:,1))
        strBRepEdge = [num2str(iPatches) '00' num2str(iBEdges)];
        if iBEdges < length(verticesMap(:,1))
            fprintf(fileHandle,[strBRepEdge blanks(2) 'TRUE' blanks(2)]);
        else
            fprintf(fileHandle,[strBRepEdge blanks(2) 'TRUE\n']);
        end
    end
    
    %% 2v. Define the Control Points of the boundary curves in the parameter space of the patch
    fprintf(fileHandle,'!==================================================================\n');
    fprintf(fileHandle,'CTRL_PTS_PAR_DEF\n');
	counterIndexCPs = counterCPs - 2*length(verticesMap(:,1)) : counterCPs - 1;
    counter = 1;
    for iCPs = 1:length(verticesMap(:,1))
        if iCPs < length(verticesMap(:,1))
            range = [iCPs iCPs + 1];
        else
            range = [iCPs 1];
        end
        for iCPsInner = range
            indexCP = counterCPs - 2*iPatches*length(verticesMap(:,1)) - 1 + iCPsInner;
            fprintf(fileHandle,[blanks(2) 'CTRLPT_PAR' blanks(2) num2str(counterIndexCPs(counter)) blanks(2)]);
            counter = counter + 1;
            for iCoord = 1:noCoord
                if iCoord < noCoord
                    fprintf(fileHandle,'%.16f  ',verticesMap(indexCP,iCoord));
                else
                    fprintf(fileHandle,'%.16f\n',verticesMap(indexCP,iCoord));
                end
            end
        end
    end
end

%% 3. Loop ove all patches and write the design block information
fprintf(fileHandle,'!###################################################################\n');
fprintf(fileHandle,'!####                        DESIGN-BLOCK                       ####\n');
fprintf(fileHandle,'!###################################################################\n');
fprintf(fileHandle,'!         ID  PART  PROP   NURBS_TOP\n');
fprintf(fileHandle,'DE-ELTOP-NURBS\n');
for iPatches = 1:noPatches
    fprintf(fileHandle,['DE-EL' blanks(3) num2str(iPatches) blanks(3) '1' blanks(3) '1' blanks(3) num2str(iPatches) '\n']);
end
fprintf(fileHandle,'!==================================================================\n');

%% 4. Write out the elements
fprintf(fileHandle,'EL-DOMAIN 1\n');
fprintf(fileHandle,'  ELEMENTS = EL-TOP 1\n');
fprintf(fileHandle,'!==================================================================\n');

%% 5. Write preamble for the weak application of conditions
fprintf(fileHandle,'!                  ID  PART  PROP   NURBS_TOP\n');
fprintf(fileHandle,'DE-BREP-NURBS\n');

%% 6. Loop over all patches and add the weak Dirichlet boundary conditions
for iPatches = 1:noPatches
    %% 6i. Get the weak Dirichlet boundary conditions of the patch
    weakDBC = BSplinePatches{iPatches}.weakDBC;
    
    %% 6i. Loop over all weak Dirichlet boundary conditions of the patch
    if ~isempty(weakDBC)
        for iWeakDBC = 1:weakDBC.noCnd
            xiExtension = weakDBC.xiExtension{iWeakDBC};
            etaExtension = weakDBC.etaExtension{iWeakDBC};
            if etaExtension(1,1) == etaExtension(1,2)
                if etaExtension(1,1) == 0
                    idCurve =  [num2str(iPatches) '001'];
                elseif etaExtension(1,1) == 1
                    idCurve =  [num2str(iPatches) '003'];
                else
                    error('Eta extension should run from 0 to 1 but not %d',etaExtension(1,1));
                end
            elseif xiExtension(1,1) == xiExtension(1,2)
                if xiExtension(1,1) == 0
                    idCurve =  [num2str(iPatches) '004'];
                elseif xiExtension(1,1) == 1
                    idCurve =  [num2str(iPatches) '002'];
                else
                    error('Xi extension should run from 0 to 1 but not %d',xiExtension(1,1));
                end
            else
                error('Wrong definition of the weak Dirichlet conditions %d at patch %d',iWeakDBC,iPatches);
            end
            idWeakDBC(counterWeakDBC,1) = str2double(idCurve);
            idWeakDBC(counterWeakDBC,2) = iPatches;
            fprintf(fileHandle,['  DE-BREP-EL' blanks(3) idCurve blanks(3) '1' blanks(3) '2' blanks(3) 'PATCH_PAR' blanks(2) idCurve blanks(2) 'OF' blanks(3) 'DE-EL' blanks(2) num2str(iPatches) '\n']);
            counterWeakDBC = counterWeakDBC + 1;
        end
    end
end

%% 7. Loop over all the connections and write out the patch connectivity information
if isfield(connections,'No')
    if connections.No > 1
        for iConn = 1:connections.No
            idI = connections.xiEtaCoup(iConn,1);
            idJ = connections.xiEtaCoup(iConn,2);
            XiI = BSplinePatches{idI}.Xi;
            EtaI = BSplinePatches{idI}.Eta;
            XiJ = BSplinePatches{idJ}.Xi;
            EtaJ = BSplinePatches{idJ}.Eta;
            xicoupI = connections.xiEtaCoup(iConn,3:4);
            etacoupI = connections.xiEtaCoup(iConn,5:6);
            xicoupJ = connections.xiEtaCoup(iConn,7:8);
            etacoupJ = connections.xiEtaCoup(iConn,9:10);    
            if etacoupI(1,1) == etacoupI(1,2)
                if etacoupI(1,1) == 0
                    idCurveI =  [num2str(idI) '001'];
                    knotStartI = XiI(1);
                    knotEndI = XiI(end);
                elseif etacoupI(1,1) == 1
                    idCurveI =  [num2str(idI) '003'];
                    knotStartI = XiI(end);
                    knotEndI = XiI(1);
                else
                    error('Eta coupling extension of patch %d should run from 0 to 1 but not %d',idI,etacoupI(1,1));
                end
            elseif xicoupI(1,1) == xicoupI(1,2)
                if xicoupI(1,1) == 0
                    idCurveI =  [num2str(idI) '004'];
                    knotStartI = EtaI(end);
                    knotEndI = EtaI(1);
                elseif xicoupI(1,1) == 1
                    idCurveI =  [num2str(idI) '002'];
                    knotStartI = EtaI(1);
                    knotEndI = EtaI(end);
                else
                    error('Xi coupling extension of patch %d should run from 0 to 1 but not %d',idI,xicoupI(1,1));
                end
            else
                error('Wrong definition of connection %d at patch %d',iConn,idI);
            end
            if etacoupJ(1,1) == etacoupJ(1,2)
                if etacoupJ(1,1) == 0
                    idCurveJ =  [num2str(idJ) '001'];
                    knotStartJ = XiJ(1);
                    knotEndJ = XiJ(end);
                elseif etacoupJ(1,1) == 1
                    idCurveJ =  [num2str(idJ) '003'];
                    knotStartJ = XiJ(end);
                    knotEndJ = XiJ(1);
                else
                    error('Eta coupling extension of patch %d should run from 0 to 1 but not %d',idJ,etacoupJ(1,1));
                end
            elseif xicoupJ(1,1) == xicoupJ(1,2)
                if xicoupJ(1,1) == 0
                    idCurveJ =  [num2str(idJ) '004'];
                    knotStartJ = EtaJ(end);
                    knotEndJ = EtaJ(1);
                elseif xicoupJ(1,1) == 1
                    idCurveJ =  [num2str(idJ) '002'];
                    knotStartJ = EtaJ(1);
                    knotEndJ = EtaJ(end);
                else
                    error('Xi coupling extension of patch %d should run from 0 to 1 but not %d',idJ,xicoupJ(1,1));
                end
            else
                error('Wrong definition of connection %d at patch %d',iConn,idJ);
            end
            idConn(iConn,1) = str2double(idCurveI);
            fprintf(fileHandle,['  DE-BREP-EL' blanks(3) idCurveI blanks(3) '1' blanks(3) '3' blanks(3) 'PATCH_PAR' blanks(2) idCurveI blanks(2) 'OF' blanks(3) 'DE-EL' blanks(2) num2str(idI) blanks(1)...
                'PATCH_PAR' blanks(2) idCurveJ blanks(2) 'OF' blanks(3) 'DE-EL' blanks(2) num2str(idJ) blanks(1) 'TOL=0.01' blanks(1) ...
                '[' num2str(knotStartI) blanks(1) num2str(knotEndI) ']' '[' num2str(knotStartJ) blanks(1) num2str(knotEndJ) ']' '\n']);
        end
    end
end

%% 8. Write out the preamble for the refinement
fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'DE-REFINEMENT\n');

%% 9. Define the refinement for each patch
for iPatches = 1:noPatches
    fprintf(fileHandle,[' DE-EL' blanks(3) num2str(iPatches) blanks(3) 'ep=0' blanks(3) 'eq=0' blanks(3) 'ru=0' blanks(3) 'rv=0\n']);
end

%% 10. Define the refinement for each weak Dirichlet condition
for iWeakDBC = 1:noWeakDBC
    fprintf(fileHandle,[' DE-BREP-EL' blanks(3) num2str(idWeakDBC(iWeakDBC,1)) blanks(3) 'dp=auto' blanks(3) 'ru=auto\n']);
end

%% 11. Define the refinement for each b-rep element associated to multipatch coupling
if isfield(connections,'No')
    if connections.No > 1
        for iConn = 1:connections.No
            fprintf(fileHandle,[' DE-BREP-EL' blanks(3) num2str(idConn(iConn,1)) blanks(3) 'dp=auto' blanks(3) 'ru=auto\n']);
        end
    end
end

%% 12. Preamble for the support conditions
fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'!        ID  DE-EL     LOC COORD  BC\n');

%% 13. Loop over all weak Dirichlet boundary conditions
for iWeakDBC = 1:noWeakDBC
    idCurve = idWeakDBC(iWeakDBC,1);
    idPatch = idWeakDBC(iWeakDBC,2);
    fprintf(fileHandle,['DE-SUP' blanks(3) num2str(iWeakDBC) blanks(3) num2str(idPatch) blanks(11) 'DE-BREP' blanks(2) num2str(idCurve) blanks(5) 'DISP_X,DISP_Y,DISP_Z !CLAMPED' '\n']);
end

%% 14. Loop over all the strong Dirichlet boundary conditions
if ~isempty(strongDBC)
    for iPatches = 1:noPatches
        strongdbc = strongDBC{iPatches};
        if ~isempty(strongdbc)
            for iStrongDBC = 1:length(strongdbc.xiExtensionHom)
                xiExtension = strongdbc.xiExtensionHom{iStrongDBC};
                etaExtension = strongdbc.etaExtensionHom{iStrongDBC};
                if etaExtension(1,1) == etaExtension(1,2) && xiExtension(1,1) == xiExtension(1,2)
                    xiKnot = xiExtension(1,1);
                    etaKnot = etaExtension(1,1);
                    fprintf(fileHandle,['DE-SUP' blanks(3) num2str(noWeakDBC + counterStrongDBC) blanks(3) num2str(iPatches) blanks(11)  'u = ' num2str(xiKnot) blanks(1) 'v = ' num2str(etaKnot) blanks(5) 'DISP_X,DISP_Y,DISP_Z !CLAMPED' '\n']);
                elseif etaExtension(1,1) == etaExtension(1,2) && xiExtension(1,1) ~= xiExtension(1,2)
                    etaKnot = etaExtension(1,1);
                    fprintf(fileHandle,['DE-SUP' blanks(3) num2str(noWeakDBC + counterStrongDBC) blanks(3) num2str(iPatches) blanks(11)  'v = ' num2str(etaKnot) blanks(5) 'DISP_X,DISP_Y,DISP_Z !CLAMPED' '\n']);
                elseif etaExtension(1,1) ~= etaExtension(1,2) && xiExtension(1,1) == xiExtension(1,2)
                    xiKnot = xiExtension(1,1);
                    fprintf(fileHandle,['DE-SUP' blanks(3) num2str(noWeakDBC + counterStrongDBC) blanks(3) num2str(iPatches) blanks(11)  'u = ' num2str(xiKnot) blanks(5) 'DISP_X,DISP_Y,DISP_Z !CLAMPED' '\n']);
                else % Strong condition applied along whole patch
                    fprintf(fileHandle,['DE-SUP' blanks(3) num2str(noWeakDBC + counterStrongDBC) blanks(3) num2str(iPatches) blanks(5) 'DISP_X,DISP_Y,DISP_Z !CLAMPED' '\n']);
                end
                counterStrongDBC = counterStrongDBC + 1;
            end
        end
    end
end
counterStrongDBC = 1;

%% 15. Define the preamble for the loads
fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'!         ID  TYPE    DE-EL     LOC COOR        D1   D2   D3     VAL\n');
fprintf(fileHandle,'!===================================================================\n');

%% 16. Loop over all patches and define the load
for iPatches = 1:noPatches
    fprintf(fileHandle,['DE-LOAD ' num2str(iPatches) ' PRES_FL ' num2str(iPatches) ' VAL=10\n']);
end

%% 17. Enable both the strong Dirichlet boundary conditions and the loads
fprintf(fileHandle,'!===================================================================\n');
fprintf(fileHandle,'LD-COM 1\n');
if ~isempty(strongDBC)
    for iPatches = 1:noPatches
        strongdbc = strongDBC{iPatches};
        if ~isempty(strongdbc)
            for iStrongDBC = 1:length(strongdbc.xiExtensionHom)
                fprintf(fileHandle,['  TYPE=BC-DIRICHLET ' num2str(noWeakDBC + counterStrongDBC) '\n']);
                counterStrongDBC = counterStrongDBC + 1;
            end
        end
    end
end
for iPatches = 1:noPatches
    fprintf(fileHandle,['  TYPE=LD-ELEM   ' num2str(iPatches) '   FAC= 1.0\n']);
end

%% 18. Close the file
fclose(fileHandle);

end
