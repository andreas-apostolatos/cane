%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   License:        BSD License                                           %
%                   cane Multiphysics default license: cane/license.txt   %
%                                                                         %
%   Main authors:   Andreas Apostolatos                                   %
%                   Marko Leskovar                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Matlab Input File                                                     %
%   _________________                                                     %
%                                                                         %
%   cane Multiphysics                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Structural Boundary Value Problem                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STRUCTURE_ANALYSIS
 ANALYSIS_TYPE,*GenData(STR_Analysis_Type)

STRUCTURE_MATERIAL_PROPERTIES
*Loop materials
*if(strcmp(MatProp(0),"Steel")==0 || strcmp(MatProp(0),"Aluminium")==0)
 DENSITY,*MatProp(Density)
 YOUNGS_MODULUS,*MatProp(Young_Modulus)
 POISSON_RATIO,*MatProp(Poisson_Ratio)
 THICKNESS,*MatProp(Thickness)
*endif
*end loop

STRUCTURE_NLINEAR_SCHEME
 NLINEAR_SCHEME,*GenData(STR_Non-Linear_Solver_Type)
 NO_LOAD_STEPS,*GenData(STR_Number_of_Steps)
 TOLERANCE,*GenData(STR_Tolerance)
 MAX_ITERATIONS,*GenData(STR_Max_Iterations)

STRUCTURE_TRANSIENT_ANALYSIS
 SOLVER *GenData(STR_Time_Analysis_Type)
 TIME_INTEGRATION *GenData(STR_Time_Integration_Scheme)
 ALPHA_BETA *GenData(STR_AlphaBeta)
 GAMMA *GenData(STR_Gamma)
 START_TIME *GenData(STR_Start_Time)
 END_TIME *GenData(STR_End_Time)
 NUMBER_OF_TIME_STEPS *GenData(STR_Number_of_Time_Steps)
 ADAPTIVE_TIME_STEPPING *GenData(STR_Adaptive_Time_Stepping)
 
STRUCTURE_INTEGRATION
 DOMAIN *GenData(STR_Gauss_Integration_Type)
 domainNoGP *GenData(STR_Domain_NO_GP)
 boundaryNoGP *GenData(STR_Boundary_NO_GP)

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

STRUCTURE_CONTACT_NODES*\
*set Cond Structure-Contact *nodes
*loop nodes OnlyInCond
*format "%8i"

*NodesNum *\
*end loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Fluid Boundary Value Problem                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FLUID_ANALYSIS
 ANALYSIS_TYPE,*GenData(CFD_Analysis_Type)

FLUID_MATERIAL_PROPERTIES
*Loop materials
*if(strcmp(MatProp(0),"Water")==0)
 DENSITY,*MatProp(Density)
 DYNAMIC_VISCOSITY,*MatProp(Dynamic_Viscosity)
*endif
*end loop

FLUID_NLINEAR_SCHEME
 NLINEAR_SCHEME,*GenData(CFD_Non-Linear_Solver_Type)
 NO_LOAD_STEPS,*GenData(CFD_Number_of_Steps)
 TOLERANCE,*GenData(CFD_Tolerance)
 MAX_ITERATIONS,*GenData(CFD_Max_Iterations)

FLUID_TRANSIENT_ANALYSIS
 SOLVER *GenData(CFD_Time_Analysis_Type)
 TIME_INTEGRATION *GenData(CFD_Time_Integration_Scheme)
 ALPHA_BETA *GenData(CFD_AlphaBeta)
 GAMMA *GenData(CFD_Gamma)
 START_TIME *GenData(CFD_Start_Time)
 END_TIME *GenData(CFD_End_Time)
 NUMBER_OF_TIME_STEPS *GenData(CFD_Number_of_Time_Steps)
 ADAPTIVE_TIME_STEPPING *GenData(CFD_Adaptive_Time_Stepping)
 
FLUID_INTEGRATION
 DOMAIN *GenData(CFD_Gauss_Integration_Type)
 domainNoGP *GenData(CFD_Domain_NO_GP)
 boundaryNoGP *GenData(CFD_Boundary_NO_GP)
	
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
*add Cond Fluid-Dirichlet-Over-Points *nodes
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
*set Cond Fluid-Dirichlet-ALE-Lines *nodes
*add Cond Fluid-Dirichlet-ALE-Points *nodes
*loop nodes OnlyInCond
*format "%8i"

*NodesNum *cond(FunctionHandleToALEMotion) *cond(FreeBoundary)*\
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Heat Boundary Value Problem                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HEAT_ANALYSIS
 ANALYSIS_TYPE,*GenData(HT_Analysis_Type)

HEAT_MATERIAL_PROPERTIES
*Loop materials
*if(strcmp(MatProp(0),"Steel")==0 || strcmp(MatProp(0),"Aluminium")==0)
 DENSITY,*MatProp(Density)
 THERMAL_CONDUCTIVITY,*MatProp(Thermal_Conductivity)
 SPECIFIC_HEAT,*MatProp(Specific_Heat)
*endif
*end loop

HEAT_TRANSIENT_ANALYSIS
 SOLVER *GenData(HT_Time_Analysis_Type)
 TIME_INTEGRATION *GenData(HT_Time_Integration_Scheme)
 START_TIME *GenData(HT_Start_Time)
 END_TIME *GenData(HT_End_Time)
 NUMBER_OF_TIME_STEPS *GenData(HT_Number_of_Time_Steps)
 ADAPTIVE_TIME_STEPPING *GenData(HT_Adaptive_Time_Stepping)
 
HEAT_INTEGRATION
 DOMAIN *GenData(HT_Gauss_Integration_Type)
 domainNoGP *GenData(HT_Domain_NO_GP)
 boundaryNoGP *GenData(HT_Boundary_NO_GP)

HEAT_NODES*\
*set Cond Heat-Nodes *nodes
*loop nodes OnlyInCond
*format "%8i%10.5f%10.5f%10.5f"

*NodesNum *NodesCoord(1,real) *NodesCoord(2,real) *NodesCoord(3,real) *\
*end loop

HEAT_ELEMENTS*\
*set Cond Heat-Elements *elems
*loop elems OnlyInCond
*format "%8i%6i%6i%8i%8i%8i%8i%8i%8i%8i%8i"

*ElemsNum *ElemsConec*\
*end loop

HEAT_DIRICHLET_NODES*\
*set Cond Heat-Dirichlet *nodes
*loop nodes OnlyInCond
*format "%8i"

*NodesNum *\
*if(cond(T-Constraint,int)==1)
*cond(T-Value)  *\
*else
NaN  *\
*endif
*end loop

HEAT_FLUX_NODES*\
*set Cond Heat-Flux *nodes
*loop nodes OnlyInCond
*format "%8i"

*if(strcmp(cond(FluxType),"boundaryFlux")==0)
*NodesNum *cond(FluxType) *cond(FunctionHandleToFluxComputation)*\
*endif
*end loop
 
 