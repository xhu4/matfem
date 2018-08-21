function result = evalTriElem(X, Y, vE, basisIdx, basisType, ders)

% take care of line segment evaluation by adding an artificial point to
% form an element
if size(vE,2)==4
	pt1 = vE(:,1:2);
	pt2 = vE(:,3:4);
	dx = pt1(:,1) - pt2(:,1);
	dy = pt1(:,2) - pt2(:,2);
	artiPt = pt1+[dy,dx(dy==0)];
	vE = [vE artiPt];
	if basisIdx==3
		basisIdx = 4;
	end
end

der_order = sum(ders);

Ex = vE(:,1:2:end);
Ey = vE(:,2:2:end);


Ex21 = Ex(:,2) - Ex(:,1);
Ex31 = Ex(:,3) - Ex(:,1);
Ey21 = Ey(:,2) - Ey(:,1);
Ey31 = Ey(:,3) - Ey(:,1);

detJ = Ex21.*Ey31-Ex31.*Ey21;
inv_detJ = 1./detJ;


switch basisType
	case 1
		switch der_order
			case 0
				switch basisIdx
					case 1
						ref_x = ((Ey31).*(X-Ex(:,1)) - (Ex31).*(Y-Ey(:,1))).*inv_detJ;
						ref_y = ((-Ey21).*(X-Ex(:,1)) + (Ex21).*(Y-Ey(:,1))).*inv_detJ;
						result = -ref_x-ref_y+1;
					case 2
						ref_x = ((Ey31).*(X-Ex(:,1)) - (Ex31).*(Y-Ey(:,1))).*inv_detJ;
						result = ref_x;
					case 3
						ref_y = ((-Ey21).*(X-Ex(:,1)) + (Ex21).*(Y-Ey(:,1))).*inv_detJ;
						result = ref_y;
					otherwise
						error('MatFem:invalidBasisIdx', 'basis index %d is not valid', basisIdx)
				end
			case 1
				switch ders(1)
					case 1
						switch basisIdx
							case 1
								result = (Ey(:,2)-Ey(:,3)).*inv_detJ;
							case 2
								result = (Ey(:,3)-Ey(:,1)).*inv_detJ;
							case 3
								result = -(Ey(:,2)-Ey(:,1)).*inv_detJ;
							otherwise
								error('MatFem:invalidBasisIdx', 'basis index %d is not valid', basisIdx)
						end
					case 0
						switch basisIdx
							case 1
								result = (Ex(:,3)-Ex(:,2)).*inv_detJ;
							case 2
								result = -(Ex(:,3)-Ex(:,1)).*inv_detJ;
							case 3
								result = (Ex(:,2)-Ex(:,1)).*inv_detJ;
							otherwise
								error('MatFem:invalidBasisIdx', 'basis index %d is not valid', basisIdx)
						end
					otherwise
						error('MatFem:invalidDerWith', 'cannot take order %d derivative with %s', der_order, der_with)
				end
			otherwise
				result = zeros(size(X));
		end
	case 2
		switch der_order
			case 0
				switch basisIdx
					case 1
						ref_x = ((Ey31).*(X-Ex(:,1)) - (Ex31).*(Y-Ey(:,1))).*inv_detJ;
						ref_y = ((-Ey21).*(X-Ex(:,1)) + (Ex21).*(Y-Ey(:,1))).*inv_detJ;
						result = 2*ref_x.^2+2*ref_y.^2+4*ref_x.*ref_y- ...
							3*ref_y-3*ref_x+1;
					case 2
						ref_x = ((Ey31).*(X-Ex(:,1)) - (Ex31).*(Y-Ey(:,1))).*inv_detJ;
						result = 2*ref_x.^2-ref_x;
					case 3
						ref_y = ((-Ey21).*(X-Ex(:,1)) + (Ex21).*(Y-Ey(:,1))).*inv_detJ;
						result = 2*ref_y.^2-ref_y;
					case 4
						ref_x = ((Ey31).*(X-Ex(:,1)) - (Ex31).*(Y-Ey(:,1))).*inv_detJ;
						ref_y = ((-Ey21).*(X-Ex(:,1)) + (Ex21).*(Y-Ey(:,1))).*inv_detJ;
						result = -4*ref_x.^2-4*ref_x.*ref_y+4*ref_x;
					case 5
						ref_x = ((Ey31).*(X-Ex(:,1)) - (Ex31).*(Y-Ey(:,1))).*inv_detJ;
						ref_y = ((-Ey21).*(X-Ex(:,1)) + (Ex21).*(Y-Ey(:,1))).*inv_detJ;
						result = 4*ref_x.*ref_y;
					case 6
						ref_x = ((Ey31).*(X-Ex(:,1)) - (Ex31).*(Y-Ey(:,1))).*inv_detJ;
						ref_y = ((-Ey21).*(X-Ex(:,1)) + (Ex21).*(Y-Ey(:,1))).*inv_detJ;
						result = -4*ref_y.^2-4*ref_x.*ref_y+4*ref_y;
					otherwise
						error('MatFem:invalidBasisIdx', 'basis index %d is not valid', basisIdx)
				end
			case 1
				switch basisIdx
					case 1
						ref_x = ((Ey31).*(X-Ex(:,1)) - (Ex31).*(Y-Ey(:,1))).*inv_detJ;
						ref_y = ((-Ey21).*(X-Ex(:,1)) + (Ex21).*(Y-Ey(:,1))).*inv_detJ;
						ref_dfdx = 4*ref_x+4*ref_y-3;
						ref_dfdy = ref_dfdx;
					case 2
						ref_x = ((Ey31).*(X-Ex(:,1)) - (Ex31).*(Y-Ey(:,1))).*inv_detJ;
						ref_dfdx = 4*ref_x-1;
						ref_dfdy = 0;
					case 3
						ref_dfdx = 0;
						ref_y = ((-Ey21).*(X-Ex(:,1)) + (Ex21).*(Y-Ey(:,1))).*inv_detJ;
						ref_dfdy = 4*ref_y-1;
					case 4
						ref_x = ((Ey31).*(X-Ex(:,1)) - (Ex31).*(Y-Ey(:,1))).*inv_detJ;
						ref_y = ((-Ey21).*(X-Ex(:,1)) + (Ex21).*(Y-Ey(:,1))).*inv_detJ;
						ref_dfdx = -8*ref_x-4*ref_y+4;
						ref_dfdy = -4*ref_x;
					case 5
						ref_x = ((Ey31).*(X-Ex(:,1)) - (Ex31).*(Y-Ey(:,1))).*inv_detJ;
						ref_y = ((-Ey21).*(X-Ex(:,1)) + (Ex21).*(Y-Ey(:,1))).*inv_detJ;
						ref_dfdx = 4*ref_y;
						ref_dfdy = 4*ref_x;
					case 6
						ref_x = ((Ey31).*(X-Ex(:,1)) - (Ex31).*(Y-Ey(:,1))).*inv_detJ;
						ref_y = ((-Ey21).*(X-Ex(:,1)) + (Ex21).*(Y-Ey(:,1))).*inv_detJ;
						ref_dfdx = -4*ref_y;
						ref_dfdy = -8*ref_y-4*ref_x+4;
					otherwise
						error('MatFem:invalidBasisIdx', 'basis index %d is not valid', basisIdx)
				end
				
				switch ders(1)
					case 1
						result = ref_dfdx.*(Ey(:,3)-Ey(:,1)).*inv_detJ+ref_dfdy.*(Ey(:,1)-Ey(:,2)).*inv_detJ;
					case 0
						result = ref_dfdx.*(Ex(:,1)-Ex(:,3)).*inv_detJ+ref_dfdy.*(Ex(:,2)-Ex(:,1)).*inv_detJ;
					otherwise
						error('MatFem:invalidDerWith', 'cannot take order %d derivative with %s', der_order, der_with)
				end
			case 2
				switch basisIdx
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
						error('MatFem:invalidBasisIdx', 'basis index %d is not valid', basisIdx)
				end
				
				switch ders(1)
					case 2
						k1 = (Ey(:,3)-Ey(:,1)).^2;
						k2 = 2*(Ey(:,3)-Ey(:,1)).*(Ey(:,1)-Ey(:,2));
						k3 = (Ey(1)-Ey(2))^2;
					case 1
						k1 = (Ex(:,1)-Ex(:,3)).*(Ey(:,3)-Ey(:,1));
						k2 = (Ex(:,1)-Ex(:,3)).*(Ey(:,1)-Ey(:,2))+(Ex(:,1)-Ex(:,2))*(Ey(:,3)-Ey(:,1));
						k3 = (Ex(:,2)-Ex(:,1)).*(Ey(:,1)-Ey(:,2));
					case 0
						k1 = (Ex(:,1)-Ex(:,3)).^2;
						k2 = 2*(Ex(:,1)-Ex(:,3))*(Ex(:,2)-Ex(:,1));
						k3 = (Ex(:,2)-Ex(:,1))^2;
					otherwise
						error('MatFem:invalidDerWith', 'cannot take order %d derivative with %s', der_order, der_with)
				end
				
				result = k1.*ref_dfdxx+k2.*ref_dfdxy+k3.*ref_dfdyy;
				result = result.*(inv_detJ.^2);
				
			otherwise
				result = zeros(size(pts,1),1);
		end
	otherwise
		error('MatFem:unsptBasisType', 'basis type %d is currently not supported', basisType)
end
if size(result, 2) ~= size(X, 2)
	result = repmat(result, 1, size(X,2));
end
end


function result = evalTriElem_backup(X, Y, vE, basisIdx, basisType, ders)

% take care of line segment evaluation by adding an artificial point to
% form an element
if size(vE,1)==2
	pt1 = vE(1,:);
	pt2 = vE(2,:);
	dx = pt1(1) - pt2(1);
	dy = pt1(2) - pt2(2);
	if dy~=0
		artiPt = pt1+[dy,0];
	elseif dx~=0
		artiPt = pt1+[0,dx];
	else
		error('Two points at the same location');
	end
	vE = [vE; artiPt];
	if basisIdx==3
		basisIdx = 4;
	end
end

der_order = sum(ders);

Ex = vE(:,1);
Ey = vE(:,2);


Ex21 = Ex(2) - Ex(1);
Ex31 = Ex(3) - Ex(1);
Ey21 = Ey(2) - Ey(1);
Ey31 = Ey(3) - Ey(1);

detJ = Ex21*Ey31-Ex31*Ey21;
inv_detJ = 1/detJ;


switch basisType
	case 1
		switch der_order
			case 0
				switch basisIdx
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
						error('MatFem:invalidBasisIdx', 'basis index %d is not valid', basisIdx)
				end
			case 1
				switch ders(1)
					case 1
						switch basisIdx
							case 1
								result = (Ey(2)-Ey(3))*inv_detJ;
							case 2
								result = (Ey(3)-Ey(1))*inv_detJ;
							case 3
								result = -(Ey(2)-Ey(1))*inv_detJ;
							otherwise
								error('MatFem:invalidBasisIdx', 'basis index %d is not valid', basisIdx)
						end
					case 0
						switch basisIdx
							case 1
								result = (Ex(3)-Ex(2))*inv_detJ;
							case 2
								result = -(Ex(3)-Ex(1))*inv_detJ;
							case 3
								result = (Ex(2)-Ex(1))*inv_detJ;
							otherwise
								error('MatFem:invalidBasisIdx', 'basis index %d is not valid', basisIdx)
						end
					otherwise
						error('MatFem:invalidDerWith', 'cannot take order %d derivative with %s', der_order, der_with)
				end
			otherwise
				result = zeros(size(X(:),1),1);
		end
	case 2
		switch der_order
			case 0
				switch basisIdx
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
						error('MatFem:invalidBasisIdx', 'basis index %d is not valid', basisIdx)
				end
			case 1
				switch basisIdx
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
						error('MatFem:invalidBasisIdx', 'basis index %d is not valid', basisIdx)
				end
				
				switch ders(1)
					case 1
						result = ref_dfdx*(Ey(3)-Ey(1))*inv_detJ+ref_dfdy*(Ey(1)-Ey(2))*inv_detJ;
					case 0
						result = ref_dfdx*(Ex(1)-Ex(3))*inv_detJ+ref_dfdy*(Ex(2)-Ex(1))*inv_detJ;
					otherwise
						error('MatFem:invalidDerWith', 'cannot take order %d derivative with %s', der_order, der_with)
				end
			case 2
				switch basisIdx
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
						error('MatFem:invalidBasisIdx', 'basis index %d is not valid', basisIdx)
				end
				
				switch ders(1)
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
						error('MatFem:invalidDerWith', 'cannot take order %d derivative with %s', der_order, der_with)
				end
				
				result = k1*ref_dfdxx+k2*ref_dfdxy+k3*ref_dfdyy;
				result = result*(inv_detJ^2);
				
			otherwise
				result = zeros(size(pts,1),1);
		end
	otherwise
		error('MatFem:unsptBasisType', 'basis type %d is currently not supported', basisType)
end


end