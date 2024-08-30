function stiffMtxElPage = ...
    computeElStiffMtxAndLoadVctPlateInMembraneActionLinearPagewise ...
    (dRdxiMatrixGPPage, dRdetaMatrixGPPage, GCovariantGPPage, ...
    G1CovariantGPPage, G2CovariantGPPage, DmGPPage)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the pagewise element stiffness matrices corresponding to the
% isogeometric plate in membrane action analysis.
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
%            G1CovariantGPPage : The base vector g1 of the reference 
%                                configuration at the Gauss point pagewise 
%                                for all the elements
%            G2CovariantGPPage : The base vector g2 of the reference 
%                                configuration at the Gauss point pagewise 
%                                for all the elements
%                     DmGPPage : The membrane material Voigt-matrix at the 
%                                Gauss point pagewise for all the elements
%
%                       Output :
%               stiffMtxElPage : The element stiffness matrices at the 
%                                Gauss point pagewise
%
% Function layout :
%
% 1. Compute the metric tensors of the reference  configuration pagewise
%
% 2. Compute the contravariant base vectors pagewise
%
% 3. Compute the local Cartesian basis pagewise
%
% 4. Compute the transformation matrix from the contravariant basis to the local Cartesian one pagewise
%
% 5. Compute the first variation of the strain Voigt-vector with respect to the DOFs in the contravariant basis pagewise
%
% 6. Compute the first variation of the strain Voigt-vector with respect to the DOFs in the local Cartesian basis pagewise
%
% 7. Compute the tangent stiffness matrix pagewise
%
%% Function main body

%% 1. Compute the metric tensors of the reference  configuration pagewise
%
% GabCovariant = GCovariant'*GCovariant;
%
GabCovariantPage = ....
	pmtimes(ptranspose(permute(GCovariantGPPage, [1 3 2])), ...
    permute(GCovariantGPPage, [1 3 2]));

%% 2. Compute the contravariant base vectors pagewise
%
% GContravariant = GabCovariant\GCovariant';
% GContravariant = GContravariant';
%
GContravariantPage = solve6x6LinearSystemPagewise...
    (GabCovariantPage, GCovariantGPPage);
GContravariantPage = ptranspose(GContravariantPage);

%% 3. Compute the local Cartesian basis pagewise
%
% eLC = computeLocalCartesianBasis4BSplineSurface(GCovariant,GContravariant);
%
eCL1Page = normr(squeeze(GCovariantGPPage(:, 1, :)));
eCL2Page = normr(squeeze(GContravariantPage(:, :, 2)));

%% 4. Compute the transformation matrix from the contravariant basis to the local Cartesian one pagewise
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
TContra2LC4VoigtStrainPage = cat(3, cat(2, pmtimes(ptranspose(GContravariantPage(:, :, 1)), eCL1Page(:, :)).^2, ...
    pmtimes(ptranspose(GContravariantPage(:, :, 1)), eCL2Page(:, :)).^2, 2.*pmtimes(ptranspose(GContravariantPage(:, :, 1)), ...
    eCL1Page(:, :)).*pmtimes(ptranspose(GContravariantPage(:, :, 1)), eCL2Page(:, :))), cat(2, pmtimes(ptranspose(GContravariantPage(:, :, 2)), ...
    eCL1Page(:, :)).^2, pmtimes(ptranspose(GContravariantPage(:, :, 2)), eCL2Page(:, :)).^2, 2.*pmtimes(ptranspose(GContravariantPage(:, :, 2)), eCL1Page(:, :)).*...
    pmtimes(ptranspose(GContravariantPage(:, :, 2)), eCL2Page(:, :))), cat(2, 2.*pmtimes(ptranspose(GContravariantPage(:, :, 1)), ...
    eCL1Page(:, :)).*pmtimes(ptranspose(GContravariantPage(:, :, 2)), eCL1Page(:, :)), 2.*pmtimes(ptranspose(GContravariantPage(:, :, 1)), eCL2Page(:, :)).*...
    pmtimes(ptranspose(GContravariantPage(:, :, 2)), eCL2Page(:, :)), 2.*pmtimes(ptranspose(GContravariantPage(:, :, 1)), eCL1Page(:, :)).* ...
    pmtimes(ptranspose(GContravariantPage(:, :, 2)), eCL2Page(:, :)) + 2.*pmtimes(ptranspose(GContravariantPage(:, :, 2)), eCL1Page(:, :)).* ...
    pmtimes(ptranspose(GContravariantPage(:, :, 1)), eCL2Page(:, :))));

%% 5. Compute the first variation of the strain Voigt-vector with respect to the DOFs in the contravariant basis pagewise
%
% dEContravariant = [gCovariant(:,1)'*dRdxiMatrix
%                    gCovariant(:,2)'*dRdetaMatrix
%                    .5*(gCovariant(:,2)'*dRdxiMatrix + gCovariant(:,1)'*dRdetaMatrix)];
%
dEContravariantPage = cat(2,pmtimes(ptranspose(G1CovariantGPPage), dRdxiMatrixGPPage), ...
    pmtimes(ptranspose(G2CovariantGPPage), dRdetaMatrixGPPage), ...
    .5*(pmtimes(ptranspose(G2CovariantGPPage), dRdxiMatrixGPPage) + pmtimes(ptranspose(G1CovariantGPPage), dRdetaMatrixGPPage)));

%% 6. Compute the first variation of the strain Voigt-vector with respect to the DOFs in the local Cartesian basis pagewise
%
% dECartesian = TFromContraToLC*dEContravariant;
%
dECartesianPage = pmtimes(TContra2LC4VoigtStrainPage,dEContravariantPage);

%% 7. Compute the tangent stiffness matrix pagewise
%
% KTEl = dECartesian'*Dm*dECartesian + ...
%     NCartesianTimesTFromContraToLC(1,1)*(dRdxiMatrix'*dRdxiMatrix) + ...
%     NCartesianTimesTFromContraToLC(1,2)*(dRdetaMatrix'*dRdetaMatrix) + ...
%     NCartesianTimesTFromContraToLC(1,3)*.5*(dRdxiMatrix'*dRdetaMatrix + dRdetaMatrix'*dRdxiMatrix);
%
stiffMtxElPage = pmtimes(pmtimes(ptranspose(dECartesianPage), DmGPPage), dECartesianPage);

end