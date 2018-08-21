function result = Eval2dTriBasis(X, Y, vE, basis_idx, basis_type, dxdy)

der_order = sum(dxdy);

Ex = vE(:,1);
Ey = vE(:,2);


Ex21 = Ex(2) - Ex(1);
Ex31 = Ex(3) - Ex(1);
Ey21 = Ey(2) - Ey(1);
Ey31 = Ey(3) - Ey(1);
 
detJ = Ex21*Ey31-Ex31*Ey21;
inv_detJ = 1/detJ;



% basisIdxE = MException('FEM2D:invalidBasisIdx', 'basis index %d is not valid', basis_idx);
% basisTypeE = MException('FEM2D:unsptBasisType', 'basis type %d is currently not supported', basis_type);
% derWithE = MException('FEM2D:invalidDerWith', 'cannot take order %d derivative with %s', der_order, der_with);

switch basis_type
	case 1
		switch der_order
			case 0
				switch basis_idx
					case 1
						ref_x = ((Ey31)*(X-Ex(1)) - (Ex31)*(Y-Ey(1)))*inv_detJ;
						ref_y = ((-Ey21)*(X-Ex(1)) + (Ex21)*(Y-Ey(1)))*inv_detJ;
						result = -ref_x-ref_y+1;
					case 2
						ref_x = ((Ey31)*(X-Ex(1)) - (Ex31)*(Y-Ey(1)))*inv_detJ;
						result = ref_x;
					case 3
						ref_y = ((-Ey21)*(X-Ex(1)) + (Ex21)*(Y-Ey(1)))*inv_detJ;
						result = ref_y;
					otherwise
						throw(MException('FEM2D:invalidBasisIdx', 'basis index %d is not valid', basis_idx))
				end
			case 1
				switch dxdy(1)
					case 1
						switch basis_idx
							case 1
								result = (Ey(2)-Ey(3))*inv_detJ;
							case 2
								result = (Ey(3)-Ey(1))*inv_detJ;
							case 3
								result = -(Ey(2)-Ey(1))*inv_detJ;
							otherwise
								throw(MException('FEM2D:invalidBasisIdx', 'basis index %d is not valid', basis_idx))
						end
					case 0
						switch basis_idx
							case 1
								result = (Ex(3)-Ex(2))*inv_detJ;
							case 2
								result = -(Ex(3)-Ex(1))*inv_detJ;
							case 3
								result = (Ex(2)-Ex(1))*inv_detJ;
							otherwise
								throw(MException('FEM2D:invalidBasisIdx', 'basis index %d is not valid', basis_idx))
						end
					otherwise
						throw(MException('FEM2D:invalidDerWith', 'cannot take order %d derivative with %s', der_order, der_with))
				end
			otherwise
				result = zeros(size(pts,1),1);
		end
	case 2
		switch der_order
			case 0
				switch basis_idx
					case 1
						ref_x = ((Ey31)*(X-Ex(1)) - (Ex31)*(Y-Ey(1)))*inv_detJ;
						ref_y = ((-Ey21)*(X-Ex(1)) + (Ex21)*(Y-Ey(1)))*inv_detJ;
						result = 2*ref_x.^2+2*ref_y.^2+4*ref_x.*ref_y- ...
							3*ref_y-3*ref_x+1;
					case 2
						ref_x = ((Ey31)*(X-Ex(1)) - (Ex31)*(Y-Ey(1)))*inv_detJ;
						result = 2*ref_x.^2-ref_x;
					case 3
						ref_y = ((-Ey21)*(X-Ex(1)) + (Ex21)*(Y-Ey(1)))*inv_detJ;
						result = 2*ref_y.^2-ref_y;
					case 4
						ref_x = ((Ey31)*(X-Ex(1)) - (Ex31)*(Y-Ey(1)))*inv_detJ;
						ref_y = ((-Ey21)*(X-Ex(1)) + (Ex21)*(Y-Ey(1)))*inv_detJ;
						result = -4*ref_x.^2-4*ref_x.*ref_y+4*ref_x;
					case 5
						ref_x = ((Ey31)*(X-Ex(1)) - (Ex31)*(Y-Ey(1)))*inv_detJ;
						ref_y = ((-Ey21)*(X-Ex(1)) + (Ex21)*(Y-Ey(1)))*inv_detJ;
						result = 4*ref_x.*ref_y;
					case 6
						ref_x = ((Ey31)*(X-Ex(1)) - (Ex31)*(Y-Ey(1)))*inv_detJ;
						ref_y = ((-Ey21)*(X-Ex(1)) + (Ex21)*(Y-Ey(1)))*inv_detJ;
						result = -4*ref_y.^2-4*ref_x.*ref_y+4*ref_y;
					otherwise
						throw(MException('FEM2D:invalidBasisIdx', 'basis index %d is not valid', basis_idx))
				end
			case 1
				switch basis_idx
					case 1
						ref_x = ((Ey31)*(X-Ex(1)) - (Ex31)*(Y-Ey(1)))*inv_detJ;
						ref_y = ((-Ey21)*(X-Ex(1)) + (Ex21)*(Y-Ey(1)))*inv_detJ;
						ref_dfdx = 4*ref_x+4*ref_y-3;
						ref_dfdy = ref_dfdx;
					case 2
						ref_x = ((Ey31)*(X-Ex(1)) - (Ex31)*(Y-Ey(1)))*inv_detJ;
						ref_dfdx = 4*ref_x-1;
						ref_dfdy = 0;
					case 3
						ref_dfdx = 0;
						ref_y = ((-Ey21)*(X-Ex(1)) + (Ex21)*(Y-Ey(1)))*inv_detJ;
						ref_dfdy = 4*ref_y-1;
					case 4
						ref_x = ((Ey31)*(X-Ex(1)) - (Ex31)*(Y-Ey(1)))*inv_detJ;
						ref_y = ((-Ey21)*(X-Ex(1)) + (Ex21)*(Y-Ey(1)))*inv_detJ;
						ref_dfdx = -8*ref_x-4*ref_y+4;
						ref_dfdy = -4*ref_x;
					case 5
						ref_x = ((Ey31)*(X-Ex(1)) - (Ex31)*(Y-Ey(1)))*inv_detJ;
						ref_y = ((-Ey21)*(X-Ex(1)) + (Ex21)*(Y-Ey(1)))*inv_detJ;
						ref_dfdx = 4*ref_y;
						ref_dfdy = 4*ref_x;
					case 6
						ref_x = ((Ey31)*(X-Ex(1)) - (Ex31)*(Y-Ey(1)))*inv_detJ;
						ref_y = ((-Ey21)*(X-Ex(1)) + (Ex21)*(Y-Ey(1)))*inv_detJ;
						ref_dfdx = -4*ref_y;
						ref_dfdy = -8*ref_y-4*ref_x+4;
					otherwise
						throw(MException('FEM2D:invalidBasisIdx', 'basis index %d is not valid', basis_idx))
				end
				
				switch dxdy(1)
					case 1
						result = ref_dfdx*(Ey(3)-Ey(1))*inv_detJ+ref_dfdy*(Ey(1)-Ey(2))*inv_detJ;
					case 0
						result = ref_dfdx*(Ex(1)-Ex(3))*inv_detJ+ref_dfdy*(Ex(2)-Ex(1))*inv_detJ;
					otherwise
						throw(MException('FEM2D:invalidDerWith', 'cannot take order %d derivative with %s', der_order, der_with))
				end
			case 2
				switch basis_idx
					case 1
						ref_dfdxx = 4;
						ref_dfdxy = 4;
						ref_dfdyy = 4;
					case 2
						ref_dfdxx = 4;
						ref_dfdxy = 0;
						ref_dfdyy = 0;
					case 3
						ref_dfdxx = 0;
						ref_dfdxy = 0;
						ref_dfdyy = 4;
					case 4
						ref_dfdxx = -8;
						ref_dfdxy = -4;
						ref_dfdyy = 0;
					case 5
						ref_dfdxx = 0;
						ref_dfdxy = 4;
						ref_dfdyy = 0;
					case 6
						ref_dfdxx = 0;
						ref_dfdxy = -4;
						ref_dfdyy = -8;
					otherwise
						throw(MException('FEM2D:invalidBasisIdx', 'basis index %d is not valid', basis_idx))
				end
				
				switch dxdy(1)
					case 2
						k1 = (Ey(3)-Ey(1))^2;
						k2 = 2*(Ey(3)-Ey(1))*(Ey(1)-Ey(2));
						k3 = (Ey(1)-Ey(2))^2;
					case 1
						k1 = (Ex(1)-Ex(3))*(Ey(3)-Ey(1));
						k2 = (Ex(1)-Ex(3))*(Ey(1)-Ey(2))+(Ex(1)-Ex(2))*(Ey(3)-Ey(1));
						k3 = (Ex(2)-Ex(1))*(Ey(1)-Ey(2));
					case 0
						k1 = (Ex(1)-Ex(3))^2;
						k2 = 2*(Ex(1)-Ex(3))*(Ex(2)-Ex(1));
						k3 = (Ex(2)-Ex(1))^2;
					otherwise
						throw(MException('FEM2D:invalidDerWith', 'cannot take order %d derivative with %s', der_order, der_with))
				end
				
				result = k1*ref_dfdxx+k2*ref_dfdxy+k3*ref_dfdyy;
				result = result*(inv_detJ^2);
				
			otherwise
				result = zeros(size(pts,1),1);
		end
	otherwise
		throw(MException('FEM2D:unsptBasisType', 'basis type %d is currently not supported', basis_type))
end