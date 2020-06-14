function [tangMtxElPage, resVctElPage] = ...
    computeIGAElTangentStiffMtxResVctMembraneNLinear ...
    (dRdxiMatrixGPPage, dRdetaMatrixGPPage, GCovariantGPPage, ...
    g1CovariantGPPage, g2CovariantGPPage, thickness, DmGPPage, ...
    prestressActPage)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the pagewise tangent stiffness matrices and residual vectors for
% all the elements corresponding to the isogeometric membrane analysis.
%
%                        Input :
%            dRdxiMatrixGPPage : The matrix containing the derivatives of 
%                                the basis functions with respect to xi at 
%                                the Gauss point pagewise for all the 
%                                elements
%           dRdetaMatrixGPPage : The matrix containing the derivatives of 
%                                the basis functions with respect to eta at 
%                                the Gauss point pagewise for all the 
%                                elements
%             GCovariantGPPage : The base vectors of the reference 
%                                configuration at the Gauss point pagewise 
%                                for all the elements
%            g1CovariantGPPage : The base vector g1 of the current 
%                                configuration at the Gauss point pagewise 
%                                for all the elements
%            g2CovariantGPPage : The base vector g2 of the current 
%                                configuration at the Gauss point pagewise 
%                                for all the elements
%                    thickness : The thickness of the membrane
%                     DmGPPage : The membrane material Voigt-matrix at the 
%                                Gauss point pagewise for all the elements
%             prestressActPage : The Voigt-vector of the actual prestress, 
%                                that is the prestress scaled by the 
%                                thickess of the membrane at the Gauss 
%                                point pagewise for all the elements
%
%                       Output :
%                tangMtxElPage : The tangent matrices at the Gauss point 
%                                pagewise for all the elements
%                 resVctElPage : The residual vectors at the Gauss point 
%                                pagewise for all the elements
%
% Function layout :
%
% 1. Compute the metric tensors of the reference  configuration pagewise
%
% 2. Compute the metric tensors of the reference  configuration pagewise
%
% 3. Compute the contravariant base vectors pagewise
%
% 4. Compute the local Cartesian basis pagewise
%
% 5. Compute the transformation matrix from the contravariant basis to the local Cartesian one pagewise
%
% 6. Compute the strain Voigt-vector in the local Cartesian basis pagewise
%
% 7. Compute the strain Voigt-vector in the local Cartesian basis pagewise
%
% 8. Compute the force Voigt-vector in the local Cartesian basis pagewise
%
% 9. Compute the first variation of the strain Voigt-vector with respect to the DOFs in the contravariant basis pagewise
%
% 10. Compute the first variation of the strain Voigt-vector with respect to the DOFs in the local Cartesian basis pagewise
%
% 11. Compute the product NCartesian*TFromContraToLC pagewise
%
% 12. Compute the tangent stiffness matrix pagewise
%
% 13. Compute the residual vector pagewise
%
%% Function main body

%% 1. Compute the metric tensors of the reference  configuration pagewise
%
% GabCovariant = GCovariant'*GCovariant;
%
GabCovariantPage = ....
	pmtimes(ptranspose(permute(GCovariantGPPage,[1 3 2])),permute(GCovariantGPPage,[1 3 2]));

%% 2. Compute the metric tensors of the reference  configuration pagewise
%
% gabCovariant = gCovariant'*gCovariant;
%
gGPCovariant = cat(3,g1CovariantGPPage,g2CovariantGPPage);
gabCovariantPage = ....
    pmtimes(ptranspose(gGPCovariant),gGPCovariant);

%% 3. Compute the contravariant base vectors pagewise
%
% GContravariant = GabCovariant\GCovariant';
% GContravariant = GContravariant';
%
GContravariantPage = solve6x6LinearSystemPagewise...
    (GabCovariantPage,GCovariantGPPage);
GContravariantPage = ptranspose(GContravariantPage);

%% 4. Compute the local Cartesian basis pagewise
%
% eLC = computeLocalCartesianBasis4BSplineSurface(GCovariant,GContravariant);
%
eCL1Page = normr(squeeze(GCovariantGPPage(:,1,:)));
eCL2Page = normr(squeeze(GContravariantPage(:,:,2)));

%% 5. Compute the transformation matrix from the contravariant basis to the local Cartesian one pagewise
%
% TFromContraToLC = computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell...
%    (eLC,GContravariant);
%
% TContra2LC4VoigtStrain11 = pmtimes(ptranspose(GContravariantPage(:,:,1)),eCL1Page(:,:)).^2;
% TContra2LC4VoigtStrain12 = pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL1Page(:,:)).^2;
% TContra2LC4VoigtStrain13 = 2.*pmtimes(ptranspose(GContravariantPage(:,:,1)),eCL1Page(:,:)).*...
%     pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL1Page(:,:));
% TContra2LC4VoigtStrain21 = pmtimes(ptranspose(GContravariantPage(:,:,1)),eCL2Page(:,:)).^2;
% TContra2LC4VoigtStrain22 = pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL2Page(:,:)).^2;
% TContra2LC4VoigtStrain23 = 2.*pmtimes(ptranspose(GContravariantPage(:,:,1)),eCL2Page(:,:)).*...
%     pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL2Page(:,:));
% TContra2LC4VoigtStrain31 = 2.*pmtimes(ptranspose(GContravariantPage(:,:,1)),eCL1Page(:,:)).*...
%     pmtimes(ptranspose(GContravariantPage(:,:,1)),eCL2Page(:,:));
% TContra2LC4VoigtStrain32 = 2.*pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL1Page(:,:)).*...
%     pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL2Page(:,:));
% TContra2LC4VoigtStrain33 = 2.*pmtimes(ptranspose(GContravariantPage(:,:,1)),eCL1Page(:,:)).*...
%     pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL2Page(:,:)) + ...
%     2.*pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL1Page(:,:)).*...
%     pmtimes(ptranspose(GContravariantPage(:,:,1)),eCL2Page(:,:));
%
TContra2LC4VoigtStrainPage = cat(3,cat(2,pmtimes(ptranspose(GContravariantPage(:,:,1)),eCL1Page(:,:)).^2,...
    pmtimes(ptranspose(GContravariantPage(:,:,1)),eCL2Page(:,:)).^2,2.*pmtimes(ptranspose(GContravariantPage(:,:,1)),...
    eCL1Page(:,:)).*pmtimes(ptranspose(GContravariantPage(:,:,1)),eCL2Page(:,:))),cat(2,pmtimes(ptranspose(GContravariantPage(:,:,2)),...
    eCL1Page(:,:)).^2,pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL2Page(:,:)).^2,2.*pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL1Page(:,:)).*...
    pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL2Page(:,:))),cat(2,2.*pmtimes(ptranspose(GContravariantPage(:,:,1)),...
    eCL1Page(:,:)).*pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL1Page(:,:)),2.*pmtimes(ptranspose(GContravariantPage(:,:,1)),eCL2Page(:,:)).*...
    pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL2Page(:,:)),2.*pmtimes(ptranspose(GContravariantPage(:,:,1)),eCL1Page(:,:)).*...
    pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL2Page(:,:)) + 2.*pmtimes(ptranspose(GContravariantPage(:,:,2)),eCL1Page(:,:)).*...
    pmtimes(ptranspose(GContravariantPage(:,:,1)),eCL2Page(:,:))));

%% 6. Compute the strain Voigt-vector in the local Cartesian basis pagewise
%
% EpsilonContravariant = .5*[gabCovariant(1,1)-GabCovariant(1,1)
%                            gabCovariant(2,2)-GabCovariant(2,2)
%                            gabCovariant(1,2)-GabCovariant(1,2)];
%
EpsilonContravariantPage = .5*cat(2,gabCovariantPage(:,1,1) - GabCovariantPage(:,1,1),...
    gabCovariantPage(:,2,2) - GabCovariantPage(:,2,2),...
    gabCovariantPage(:,1,2) - GabCovariantPage(:,1,2));

%% 7. Compute the strain Voigt-vector in the local Cartesian basis pagewise
%
% ELocalCartesian = TFromContraToLC*EpsilonContravariant;
%
ELocalCartesianPage = pmtimes(TContra2LC4VoigtStrainPage,EpsilonContravariantPage);

%% 8. Compute the force Voigt-vector in the local Cartesian basis pagewise
%
% NCartesian = prestress + Dm*ELocalCartesian;
%
NCartesianPage = pstimes(prestressActPage,thickness) + pmtimes(DmGPPage,ELocalCartesianPage);

%% 9. Compute the first variation of the strain Voigt-vector with respect to the DOFs in the contravariant basis pagewise
%
% dEContravariant = [gCovariant(:,1)'*dRdxiMatrix
%                    gCovariant(:,2)'*dRdetaMatrix
%                    .5*(gCovariant(:,2)'*dRdxiMatrix + gCovariant(:,1)'*dRdetaMatrix)];
%
dEContravariantPage = cat(2,pmtimes(ptranspose(g1CovariantGPPage),dRdxiMatrixGPPage),...
    pmtimes(ptranspose(g2CovariantGPPage),dRdetaMatrixGPPage),...
    .5*(pmtimes(ptranspose(g2CovariantGPPage),dRdxiMatrixGPPage) + pmtimes(ptranspose(g1CovariantGPPage),dRdetaMatrixGPPage)));

%% 10. Compute the first variation of the strain Voigt-vector with respect to the DOFs in the local Cartesian basis pagewise
%
% dECartesian = TFromContraToLC*dEContravariant;
%
dECartesianPage = pmtimes(TContra2LC4VoigtStrainPage,dEContravariantPage);

%% 11. Compute the product NCartesian*TFromContraToLC pagewise
%
% NCartesianTimesTFromContraToLC = NCartesian'*TFromContraToLC;
%
NCartesianTimesTFromContraToLCPage = pmtimes(ptranspose(NCartesianPage),TContra2LC4VoigtStrainPage);

%% 12. Compute the tangent stiffness matrix pagewise
%
% KTEl = dECartesian'*Dm*dECartesian + ...
%     NCartesianTimesTFromContraToLC(1,1)*(dRdxiMatrix'*dRdxiMatrix) + ...
%     NCartesianTimesTFromContraToLC(1,2)*(dRdetaMatrix'*dRdetaMatrix) + ...
%     NCartesianTimesTFromContraToLC(1,3)*.5*(dRdxiMatrix'*dRdetaMatrix + dRdetaMatrix'*dRdxiMatrix);
%
tangMtxElPage = pmtimes(pmtimes(ptranspose(dECartesianPage),DmGPPage),dECartesianPage) + ...
    pstimes(pmtimes(ptranspose(dRdxiMatrixGPPage),dRdxiMatrixGPPage),NCartesianTimesTFromContraToLCPage(:,1,1)) + ...
    pstimes(pmtimes(ptranspose(dRdetaMatrixGPPage),dRdetaMatrixGPPage),NCartesianTimesTFromContraToLCPage(:,1,2)) + ...
    pstimes(.5*(pmtimes(ptranspose(dRdxiMatrixGPPage),dRdetaMatrixGPPage) + pmtimes(ptranspose(dRdetaMatrixGPPage),dRdxiMatrixGPPage)),...
    NCartesianTimesTFromContraToLCPage(:,1,3));
   
%% 13. Compute the residual vector pagewise
%
% FIEl = dECartesian'*NCartesian;
%
resVctElPage = pmtimes(ptranspose(dECartesianPage),NCartesianPage);

end
