function spc = mesh2spc(mesh, basisType, order)
if nargin < 3
	order = basisType;
end

switch basisType
	case 1
		Pb = mesh.P;
		Tb = mesh.T;
	case 2
		sides = [mesh.T(:,1) mesh.T(:,2);
			mesh.T(:,2) mesh.T(:,3);
			mesh.T(:,1) mesh.T(:,3)];
		
		sides = sort(sides, 2);
		[sides, ~, ic] = unique(sides,'rows');
		
		newP = (mesh.P(sides(:,1),:) + mesh.P(sides(:,2),:)) / 2;
		
		% use sparse matrix to index new points
		Pb = [mesh.P; newP];
		Tb = [mesh.T reshape(ic+mesh.nv, [], 3)];
end
spc = MatFem.FcnSpc(Pb, Tb, basisType, mesh.nlv, order);
end