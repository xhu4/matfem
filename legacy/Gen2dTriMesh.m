function [Pb, Tb] = Gen2dTriMesh(range, h_x, h_y, basis_type)
% GEN2DTRIMESH Generates 2D triangle mesh
%	
%	[Pb, Tb] = GEN2DTRIMESH([0, 1, 0, 2], 0.2, 0.4, 2) 
%
%	See also GENBOUNDARYEDGES, GENBOUNDARYNODES

left = range(1);
right = range(2);
bottom = range(3);
top = range(4);

inv_bt = 1/basis_type;

assert(mod(right-left,h_x)==0 && mod(top-bottom,h_y)==0);

Py = bottom:h_y*inv_bt:top;
Px = left:h_x*inv_bt:right;

n_nodes_x = numel(Px);
n_nodes_y = numel(Py);

Pb(2,:) = repmat(Py, size(Px));
if isOctave
	Pb(1,:) = repelems(Px, [1:numel(Px);n_nodes_y*ones(size(Px))]);
else
	Pb(1,:) = repelem(Px, n_nodes_y);
end

% idx_map = reshape(1:n_nodes_x*n_nodes_y,[n_nodes_y, n_nodes_x]);

n_square_x = (n_nodes_x-1)*inv_bt;
n_square_y = (n_nodes_y-1)*inv_bt;
n_ele = n_square_x*n_square_y*2;

Tb = zeros(sum(1:basis_type+1), n_ele);

% Consider the first square
% NW---N---NE
% |	 \ |    |
% |   \|    |
% W----C----E
% |	   | \  |
% |    |  \ |
% SW---S---SE

SW = 1;
NW = SW+basis_type;
SE = SW+n_nodes_y*basis_type;
NE = SE+basis_type;

Tb(1:3,1) = [SW; SE; NW];
Tb(1:3,2) = [NW; SE; NE];

switch basis_type
	case 1
	case 2
		S = 1+n_nodes_y;
		C = S+1;
		N = S+2;
		W = SW+1;
		E = SE+1;
		Tb(4:6,1) = [S, C, W];
		Tb(4:6,2) = [C, E, N];
	otherwise
		throw(MException('FEM:basisNotSupported',['basis_type %d is not '...
			'supported in Gen2dTriMesh.'], basis_type));
end

% copy the first square to the first column
for i = 2:n_square_y
	Tb(:,[2*i-1,2*i]) = Tb(:,[1,2])+(i-1)*basis_type;
end

% copy the first column to the whole mesh
n_ele_y = n_square_y*2;
for i = 2:n_square_x
	Tb(:,((i-1)*n_ele_y+1):i*n_ele_y) = Tb(:,1:n_square_y*2)+n_nodes_y*(i-1)*basis_type;
end

end
		
function retval = isOctave
  persistent cacheval;  % speeds up repeated calls

  if isempty (cacheval)
    cacheval = (exist ('OCTAVE_VERSION', 'builtin') > 0);
  end

  retval = cacheval;
end
