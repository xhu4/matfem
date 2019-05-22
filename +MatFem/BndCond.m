classdef BndCond
% BNDCOND controls boundaries of a mesh.
	properties (SetAccess = private)
		% a matrix containing information of nodes on boundary,
		%	bndNodes(i,1) = condition type of the ith boundary node
		%	bndNodes(i,2) = index of the ith boundary node
		bn
		
		% a matrix containing information of edges on boundary, in the form
		% bndEdges(i, :)=
		%	[condition type of edge i ...
		%	 indices of element containing edge i ...
		%	 indices of the start node of edge i ...
		%	 indices of the end node of edge i ...
		%	 indices of the center node of edge i if basisType=2]
		be
		
		% the function space that this boundary condition is based on
		spc
		
		nlv % # of local vertices in each edge
		
		quad
		wt
		len
		v
	end
	
	methods
		function obj = BndCond(bn, be, nlv, spc)
			obj.bn = bn;
			obj.be = be;
			obj.spc = spc;
			obj.nlv = nlv;
			if ~isempty(be)
% 				obj.v = vOfBnd(obj);
% 				[obj.quad, obj.wt, obj.len] = ...
% 					MatFem.quadPts(spc.order, obj.v, spc.ndim);
			end
		end
		
		
		function [i, X, Y] = nodes(obj, marker)
			i = obj.bn((obj.bn(:,1)==marker),2);
			[X, Y] = obj.spc.vertices(i);
		end
		
		
		function [b, A, M] = applyDir(obj, marker, f, b, A, M, shift)
		%APPLYDIR apply Dirchlet boundary condition to boundaries with
		% given marker.
		%
		% Args:
		%	marker:	marker of the Dirchlet boundary.
		%	f:		the values on the boundary. Can be scalar or function
		%			handle.
		%	b:		optional. Right hand side of the assembled system.
		%	A:		optional. Left hand side matrix of the assembled system.
		%	M:		optional. Mass matrix.
		%	shift:	optional. shift the dofs in row and column indeces to 
		%			apply Dirchlet condition. Default to 0. If is a scalar,
		%			shift = [shift shift].
		%
		% See also APPLYEDGE
			if nargin < 7
				shift = [0 0];
			elseif length(shift)==1
				shift = [shift shift];
			end
			
			if ischar(marker)
				switch marker
					case {'d', 'D'}
						marker=1;
					case {'n', 'N'}
						marker=2;
					case {'r', 'R'}
						marker=3;
				end
			end
			
			bn_ = obj.bn((obj.bn(:,1)==marker), 2:end);
			bn_i = bn_+shift(1);
			bn_j = bn_+shift(2);
			
			if ~isempty(b)
				pts = obj.spc.Pb(bn_,:);
				if isa(f, 'function_handle')
					b(bn_i,:) = f(pts(:,1),pts(:,2));
				else
					b(bn_i,:) = f;
				end
			end
			
			if nargin >= 5 && ~isempty(A) && nargout>1
% 				% Slow for sparse matrix
% 				A(bn_,:)=0;
% 				A(sub2ind(size(A), bn_, bn_)) = 1;
				[i, j, val] = find(A);
				idx = ismember(i, bn_i);
				i(idx) = [];
				j(idx) = [];
				val(idx) = [];
				i = [i;bn_i];
				j = [j;bn_j];
				val = [val;ones(size(bn_))];
				A = sparse(i,j,val);
			else
				A=[];
			end
			
			if nargin==6 && ~isempty(M) && nargout>2
				M(bn_i,:)=0;
			else
				M=[];
			end
		end
		
		function r = applyEdge(bnd, marker, ders, f)
		%APPLYEDGE similar to ASSEMBLE, but integrate on boundary with
		% given marker.
		%
		% Args:
		%	bnd: the BndCond object;
		%	marker: the marker of boundaries where the integration happens;
		%	ders: derivatives of test [and trial] functions, similar to
		%	that of function `assemble`.
		%	f: additional functions to be integrated. Default to 1.
		%
		% See also ASSEMBLE.
		
			if nargin < 4
				f = 1;
			end
			Tb = bnd.be(bnd.be(:,1)==marker,2:end);
			% FIXME: not correct if ders are not 0.
			bspc = MatFem.FcnSpc(bnd.spc.Pb, Tb, bnd.spc.basisType, ...
								 bnd.nlv, bnd.spc.order);
			r = MatFem.assemble(bspc, ders, f);
		end
		
	end
end
