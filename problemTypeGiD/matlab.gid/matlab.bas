%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   MSc.-Ing. Aditya Ghantasala        (aditya.ghantasala@tum.de)         %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Matlab Input File                                                     %
%   _________________                                                     %
%                                                                         %
%   FiniteElementAnalysisProgramStructuralAnalysisInstituteTUM            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Structural Boundary Value Problem                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STRUCTURE_ANALYSIS
 ANALYSIS_TYPE,*GenData(STR_ANA-TYPE)

STRUCTURE_MATERIAL_PROPERTIES
*Loop materials
*if(strcmp(MatProp(0),"Structure")==0)
  DENSITY,*MatProp(Density)
  YOUNGS_MODULUS,*MatProp(Young_Modulus)
  POISSON_RATIO,*MatProp(Poisson_Ratio)
*endif
*end loop

STRUCTURE_NLINEAR_SCHEME
 NLINEAR_SCHEME,*GenData(STR_NL_SOLVER-TYPE)
 TOLERANCE,*GenData(STR_TOL)
 MAX_ITERATIONS,*GenData(STR_MAX_IT)

STRUCTURE_TRANSIENT_ANALYSIS
 SOLVER *GenData(STR_TIME_ANA-TYPE)
 TIME_INTEGRATION *GenData(STR_TIME_INTEGRATION-SCHEME)
 START_TIME *GenData(STR_START_TIME)
 END_TIME *GenData(STR_END_TIME)
 NUMBER_OF_TIME_STEPS *GenData(STR_NUMBER_OF_TIME_STEPS)

STRUCTURE_NODES*\
*set Cond Structure-Nodes *nodes
*loop nodes OnlyInCond
*format "%8i%10.5f%10.5f%10.5f"

*NodesNum *NodesCoord(1,real) *NodesCoord(2,real) *NodesCoord(3,real) *\
*end loop

STRUCTURE_ELEMENTS*\
*set Cond Structure-Elements *elems
*loop elems OnlyInCond
*format "%8i%6i%6i%8i%8i%8i%8i%8i%8i%8i%8i"

*ElemsNum *ElemsConec*\
*end loop

STRUCTURE_DIRICHLET_NODES*\
*set Cond Structure-Dirichlet *nodes
*loop nodes OnlyInCond
*format "%8i"

*NodesNum *\
*if(cond(X-Constraint,int)==1)
*cond(X-Value)  *\
*else
NaN  *\
*endif
*if(cond(Y-Constraint,int)==1)
*cond(Y-Value)  *\
*else
NaN  *\
*endif
*if(cond(Z-Constraint,int)==1)
*cond(Z-Value)  *\
*else
NaN  *\
*endif
*end loop

STRUCTURE_FORCE_NODES*\
*set Cond Structure-Force *nodes
*loop nodes OnlyInCond
*format "%8i"

*if(strcmp(cond(ForceType),"boundaryLoad")==0)
*NodesNum *cond(ForceType) *cond(FunctionHandleToForceComputation)*\
*endif
*end loop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Fluid Boundary Value Problem                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FLUID_ANALYSIS
 ANALYSIS_TYPE,*GenData(CFD_ANA-TYPE)

FLUID_MATERIAL_PROPERTIES
*Loop materials
*if(strcmp(MatProp(0),"Fluid")==0)
  DENSITY,*MatProp(Density)
  DYNAMIC_VISCOSITY,*MatProp(Dynamic_Viscosity)
*endif
*end loop

FLUID_NLINEAR_SCHEME
 NLINEAR_SCHEME,*GenData(CFD_NL_SOLVER-TYPE)
 NO_LOAD_STEPS,*GenData(CFD_N_STEPS)
 TOLERANCE,*GenData(CFD_TOL)
 MAX_ITERATIONS,*GenData(CFD_MAX_IT)

FLUID_TRANSIENT_ANALYSIS
 SOLVER *GenData(CFD_TIME_ANA-TYPE)
 TIME_INTEGRATION *GenData(CFD_TIME_INTEGRATION-SCHEME)
 ALPHA_BETA *GenData(CFD_ALPHA/BETA)
 GAMMA *GenData(CFD_GAMMA)
 START_TIME *GenData(CFD_START_TIME)
 END_TIME *GenData(CFD_END_TIME)
 NUMBER_OF_TIME_STEPS *GenData(CFD_NUMBER_OF_TIME_STEPS)
 ADAPTIVE_TIME_STEPPING *GenData(CFD_ADAPTIVE_TIME_STEPPING)
 
FLUID_INTEGRATION
 DOMAIN *GenData(CFD_DOMAIN_TYPE-TYPE)
 domainNoGP *GenData(CFD_DOMAIN_NO_GP)
 boundaryNoGP *GenData(CFD_BOUNDARY_NO_GP)
	
FLUID_ELEMENTS
*set Cond Fluid-Elements-Over-Surfaces *elems
*loop elems OnlyInCond
*format "%8i%6i%6i%8i%8i%8i%8i%8i%8i%8i%8i"
*ElemsNum *ElemsConec
*end loop

*set Cond Fluid-Elements-Over-Volumes *elems
*loop elems OnlyInCond
*format "%8i%6i%6i%8i%8i%8i%8i%8i%8i%8i%8i"
*ElemsNum *ElemsConec
*end loop

FLUID_DIRICHLET_NODES*\
*set Cond Fluid-Dirichlet-Over-Lines *nodes
*loop nodes OnlyInCond
*format "%8i"

*NodesNum *\
*if(cond(X-Constraint,int)==1)
*cond(X-Value)  *\
*else
NaN  *\
*endif
*if(cond(Y-Constraint,int)==1)
*cond(Y-Value)  *\
*else
NaN  *\
*endif
*if(cond(Z-Constraint,int)==1)
*cond(Z-Value)  *\
*else
NaN  *\
*endif
*if(cond(P-Constraint,int)==1)
*cond(P-Value)  *\
*else
NaN  *\
*endif
*end loop

*set Cond Fluid-Dirichlet-Over-Surfaces *nodes
*loop nodes OnlyInCond
*format "%8i"

*NodesNum *\
*if(cond(X-Constraint,int)==1)
*cond(X-Value)  *\
*else
NaN  *\
*endif
*if(cond(Y-Constraint,int)==1)
*cond(Y-Value)  *\
*else
NaN  *\
*endif
*if(cond(Z-Constraint,int)==1)
*cond(Z-Value)  *\
*else
NaN  *\
*endif
*if(cond(P-Constraint,int)==1)
*cond(P-Value)  *\
*else
NaN  *\
*endif
*end loop

FLUID_DIRICHLET_ALE_NODES*\
*set Cond Fluid-Dirichlet-ALE *nodes
*loop nodes OnlyInCond
*format "%8i"

*NodesNum *cond(FunctionHandleToALEMotion)*\
*end loop

FLUID_POST_PROC_NODES*\
*set Cond Fluid-PostProcessing-Over-Lines *nodes
*loop nodes OnlyInCond
*format "%8i"

*NodesNum *cond(BodyNumber) *cond(FunctionHandleToPostProcessing)*\
*end loop

*set Cond Fluid-PostProcessing-Over-Surfaces *nodes
*loop nodes OnlyInCond
*format "%8i"

*NodesNum *cond(BodyNumber) *cond(FunctionHandleToPostProcessing)*\
*end loop

FLUID_NODES
*set Cond Fluid-Nodes-Over-Surfaces *nodes
*loop nodes OnlyInCond
*format "%8i%10.5f%10.5f%10.5f"

*NodesNum *NodesCoord(1,real) *NodesCoord(2,real) *NodesCoord(3,real) *\
*end loop

*set Cond Fluid-Nodes-Over-Volumes *nodes
*loop nodes OnlyInCond
*format "%8i%10.5f%10.5f%10.5f"

*NodesNum *NodesCoord(1,real) *NodesCoord(2,real) *NodesCoord(3,real) *\
*end loop