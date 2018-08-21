classdef FcnSpc < handle
	% RectFcnSpc constructs a rectangular function space.
	
	properties (SetAccess = private)
		basisType	% basis type, 1 = linear, 2 = quadratic
		
		% node coordinate matrix, Pb(i, :) = [xi yi] where (xi, yi) is the
		% coordinate of the ith node.
		Pb
		
		% element matrix, Tb(i, :) = [i1 i2 i3 ...] where i's are
		% indices of nodes in ith element.
		Tb
		
		nlv			% # of local vertices in each element
		
		quad
		wt
		area
		v
		
		order
	end
	
	properties (Dependent)
		nlb			% # of local bases in each element
		nb			% # of total bases
		ne			% # of total elements
		ndim		% # of dimension
	end
	
	methods 
		function obj = FcnSpc(Pb, Tb, basisType, nlv, order)
			obj.basisType = basisType;
			obj.order = order;
			obj.nlv = nlv;
			
			obj.Pb = Pb;
			obj.Tb = Tb;
			
			obj.v = vOfElems(obj);
			[obj.quad, obj.wt, obj.area] = ...
				MatFem.quadPts(order, obj.v, obj.ndim);
		end
		
		function [vtc, idx] = elem(spc, i)
			idx = spc.Tb(i,:);
			vtc = spc.Pb(idx,:);
		end
		
		function r = get.nlb(spc)
			r = size(spc.Tb, 2);
		end
		
		function r = get.ne(spc)
			r = size(spc.Tb, 1);
		end
		
		function r = get.nb(spc)
			r = size(spc.Pb, 1);
		end
		
		function r = get.ndim(spc)
			r = size(spc.Pb, 2);
		end
		
		function [varargout] = vertices(spc,i)
			if nargin < 2
				varargout = num2cell(spc.Pb, 1);
			else
				varargout = num2cell(spc.Pb(i,:),1);
			end
		end
		
		function r = project(spc, u)
			if isa(u, 'function_handle')
				X = num2cell(spc.Pb, 1);
				r = u(X{:});
			elseif isnumeric(u)
				r = u;
			else
				error('cannot project u into space');
			end
		end
		
		function r = evalAll(spc, der)
			if isequal(der, [0 0]) && spc.nlv==3
				r = MatFem.phi(spc.order, spc.basisType);
				r = mat2cell(r', ones(1,spc.nlb), size(r,1));
			else
				r = cell(spc.nlb,1);
				[X, Y] = spc.quad{:};
				for i = 1:spc.nlb
					r{i} = MatFem.evalTriElem(X, Y, spc.v, i, spc.basisType, der);
				end
			end
		end
		
		function r = evalOnQuad(spc, u, der)
			if nargin < 3 || isempty(der)
				der = [0 0];
			end
			val = evalAll(spc, der);
			r = zeros(size(val{1}));
			for i = 1:spc.nlb
				r = r + val{i}.* u(spc.Tb(:,i));
			end
		end
		
		function plot(spc, varargin)
			patch('Faces', spc.Tb(:,1:spc.nlv), 'Vertices', spc.Pb, ...
				'FaceColor', 'White', varargin{:});
			hold on;
			plot(spc.quad{1}(:), spc.quad{2}(:), 'r.');
			hold off;
		end
		
		
		function quiver(spc, u, v)
			uu = spc.project(u);
			vv = spc.project(v);
			X = num2cell(spc.Pb, 1);
			quiver(X{:}, uu, vv, 'AutoScale', 'off');
		end
		
		function plotu(spc, u, varargin)
			uu = spc.project(u);
			if ~ismember('EdgeColor', varargin)
				edgecolor = {'EdgeColor', [.6 .6 .6]};
			else
				edgecolor = {};
			end
			patch('Faces', spc.Tb(:, 1:spc.nlv), 'Vertices', spc.Pb, ...
				'FaceVertexCData', uu, 'FaceColor', 'interp', ...
				edgecolor{:}, varargin{:});
			colorbar
		end
	end
	
end

function v = vOfElems(spc)
nd = spc.ndim;
nlv = spc.nlv;
v = spc.Pb(spc.Tb(:,1:nlv),:);
v = mat2cell(v, spc.ne*ones(nlv, 1), nd);
v = cell2mat(v');
end

