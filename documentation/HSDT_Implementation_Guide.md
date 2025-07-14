# HSDT Implementation Guide for cane Multiphysics

## Overview

This guide documents the implementation of **Higher-order Shear Deformation Theory (HSDT)** in the cane Multiphysics framework with **5 degrees of freedom per control point**.

### Traditional vs. HSDT Formulations

| Formulation | DOF per Control Point | Variables | Theory |
|-------------|----------------------|-----------|--------|
| **Kirchhoff-Love** | 3 | [u, v, w] | Thin shell, no transverse shear |
| **HSDT** | 5 | [u, v, w, θx, θy] | Thick shell, includes transverse shear |

## Theory Background

### HSDT Kinematics

HSDT extends classical shell theory by:

1. **Independent rotation variables**: θx, θy are independent of w derivatives
2. **Transverse shear deformation**: γxz, γyz are computed from displacement gradients and rotations
3. **Shear correction factor**: Accounts for non-uniform shear distribution through thickness

### Strain-Displacement Relations

#### Membrane Strains
```
εxx = ∂u/∂x
εyy = ∂v/∂y  
γxy = ∂u/∂y + ∂v/∂x
```

#### Bending Strains (Curvatures)
```
κxx = -∂θx/∂x
κyy = -∂θy/∂y
κxy = -(∂θx/∂y + ∂θy/∂x)
```

#### Shear Strains
```
γxz = ∂w/∂x + θx
γyz = ∂w/∂y + θy
```

### Material Matrices

#### Membrane Stiffness
```
Dm = (Et)/(1-ν²) * [1   ν   0  ]
                   [ν   1   0  ]
                   [0   0  (1-ν)/2]
```

#### Bending Stiffness  
```
Db = (Et³)/(12(1-ν²)) * [1   ν   0  ]
                        [ν   1   0  ]
                        [0   0  (1-ν)/2]
```

#### Shear Stiffness
```
Ds = κ*G*t * [1  0]    where κ = shear correction factor (typically 5/6)
             [0  1]           G = E/(2(1+ν))
```

## Implementation Structure

### New Files Created

#### Core Functions
- `computeElStiffMtxHSDTShellLinear.m` - Element stiffness matrix computation
- `computeStiffMtxAndLoadVctIGAHSDTShellLinear.m` - Global assembly
- `solve_IGAHSDTShellLinear.m` - Main solver

#### Supporting Functions  
- `findDofs5D_HSDT.m` - DOF identification for boundary conditions
- `fillUpPatch_HSDT.m` - Patch structure initialization
- `computeLoadVctAreaIGAHSDTShell.m` - Load vector computation

#### Analysis Scripts
- `main_scordelisLoRoof_HSDT.m` - Example analysis using HSDT

## Usage Instructions

### 1. Basic HSDT Analysis Setup

```matlab
%% Material Parameters
parameters.E = 2.1e11;           % Young's modulus [Pa]
parameters.nue = 0.3;            % Poisson's ratio
parameters.t = 0.01;             % Thickness [m]  
parameters.shearCorrection = 5/6; % Shear correction factor

%% Analysis Type
analysis.type = 'isogeometricHSDTShellAnalysis';

%% Create HSDT Patch
BSplinePatch = fillUpPatch_HSDT...
    (analysis, p, Xi, q, Eta, CP, isNURBS, parameters, ...
     homDOFs, inhomDOFs, valuesInhomDOFs, weakDBC, cables, NBC, int);

%% Solve
[dHat, F, minElArea] = solve_IGAHSDTShellLinear...
    (BSplinePatch, solve_LinearSystem, 'outputEnabled');
```

### 2. Boundary Conditions for HSDT

#### Constraint Types (5 DOF per control point)
```matlab
% DOF ordering: [u, v, w, θx, θy]
% dirSupp values:
%   1 = u-displacement
%   2 = v-displacement  
%   3 = w-displacement
%   4 = θx-rotation
%   5 = θy-rotation

% Example: Fix all displacements at boundary
xiSup = [0 0]; etaSup = [0 1];
for dirSupp = [1 2 3]  % constrain u, v, w
    homDOFs = findDofs5D_HSDT(homDOFs, xiSup, etaSup, dirSupp, CP);
end

% Example: Fix rotations at supports
for dirSupp = [4 5]  % constrain θx, θy  
    homDOFs = findDofs5D_HSDT(homDOFs, xiSup, etaSup, dirSupp, CP);
end
```

### 3. Result Interpretation

#### DOF Extraction
```matlab
% Extract displacement components
numCPs = nxi * neta;
u_displ = dHat(1:5:end);      % u-displacements
v_displ = dHat(2:5:end);      % v-displacements
w_displ = dHat(3:5:end);      % w-displacements  
theta_x = dHat(4:5:end);      % θx-rotations [rad]
theta_y = dHat(5:5:end);      % θy-rotations [rad]
```

#### Result Analysis
```matlab
% Maximum values
max_w_displacement = max(abs(w_displ));
max_rotation_x = max(abs(theta_x));
max_rotation_y = max(abs(theta_y));

fprintf('Max |w|: %.6e\n', max_w_displacement);
fprintf('Max |θx|: %.6e rad (%.3f deg)\n', max_rotation_x, rad2deg(max_rotation_x));
fprintf('Max |θy|: %.6e rad (%.3f deg)\n', max_rotation_y, rad2deg(max_rotation_y));
```

## Key Differences from Kirchhoff-Love

### DOF Structure
- **Kirchhoff-Love**: 3 DOF per CP → Total DOFs = 3 × numCPs
- **HSDT**: 5 DOF per CP → Total DOFs = 5 × numCPs

### Strain Computation
- **Kirchhoff-Love**: Curvatures from w second derivatives
- **HSDT**: Curvatures from rotation gradients, includes shear strains

### Material Properties
- **Kirchhoff-Love**: Membrane + bending stiffness
- **HSDT**: Membrane + bending + shear stiffness

### Boundary Conditions
- **Kirchhoff-Love**: 3 constraints possible per CP
- **HSDT**: 5 constraints possible per CP (including rotations)

## Shear Correction Factor Guidelines

| Shell Type | Recommended κ | Reference |
|------------|---------------|-----------|
| **Reissner-Mindlin** | 5/6 | Classical value |
| **Thick shells** | π²/12 ≈ 0.822 | Refined theory |
| **Very thick shells** | 1.0 | No correction |
| **Composite shells** | Variable | Material dependent |

## Performance Considerations

### Computational Cost
- **Memory**: 67% increase in DOFs (5/3 ratio)
- **Stiffness matrix**: 178% increase in size ((5/3)² ratio)  
- **Assembly time**: ~67% increase
- **Solution time**: ~178% increase (depends on solver)

### Convergence
- **HSDT generally converges faster** for thick shells
- **Better behavior** in bending-dominated problems
- **More accurate** transverse shear stresses

## Validation and Testing

### Test Cases
1. **Scordelis-Lo Roof**: Compare HSDT vs. Kirchhoff-Love
2. **Thick plate bending**: Validate shear deformation effects
3. **Simply supported plate**: Compare with analytical solutions

### Expected Results
- **Thin shells (t/L < 1/100)**: HSDT ≈ Kirchhoff-Love
- **Thick shells (t/L > 1/20)**: HSDT shows softer response
- **Shear-dominated cases**: Significant differences

## Troubleshooting

### Common Issues

#### 1. DOF Mismatch Errors
```
Error: Load vector size mismatch. Expected XXXXX DOFs for HSDT
```
**Solution**: Ensure all functions use 5 DOF per control point consistently

#### 2. Overconstrained System
```
Error: System is overconstrained
```
**Solution**: Review boundary conditions - avoid fixing all 5 DOFs unnecessarily

#### 3. Shear Locking
**Symptoms**: Overly stiff response, poor convergence
**Solution**: 
- Reduce shear correction factor
- Use selective reduced integration
- Increase element order

#### 4. Ill-Conditioned Matrices
**Symptoms**: Poor convergence, numerical instabilities
**Solution**:
- Check material parameters (avoid ν → 0.5)
- Verify element quality
- Use appropriate boundary conditions

### Debugging Tips

1. **Check DOF ordering**: Verify [u, v, w, θx, θy] sequence
2. **Validate material matrices**: Ensure positive definiteness
3. **Monitor rotations**: Large rotations may indicate modeling issues
4. **Compare with Kirchhoff-Love**: For validation in thin shell limit

## Future Extensions

### Planned Enhancements
1. **Nonlinear HSDT**: Geometric and material nonlinearity
2. **Dynamic analysis**: Mass matrix for transient problems
3. **Multilayer shells**: Layer-wise HSDT formulation
4. **Postprocessing**: Stress recovery and visualization

### Research Directions
1. **Locking mitigation**: Advanced integration schemes
2. **Adaptive refinement**: Error-driven mesh adaptation
3. **Multiphysics coupling**: Thermal-structural interaction
4. **Optimization**: Shape and topology optimization with HSDT

## References

1. J.N. Reddy, "Mechanics of Laminated Composite Plates and Shells"
2. T.J.R. Hughes, "The Finite Element Method"
3. J. Kiendl, "Isogeometric Analysis and Shape Optimal Design for Shell Structures"
4. F. Auricchio et al., "Isogeometric collocation for elastostatics and explicit dynamics"

---

**Note**: This implementation provides a solid foundation for HSDT analysis in the cane Multiphysics framework. Users should validate results against known benchmarks and consult literature for specific application guidelines.