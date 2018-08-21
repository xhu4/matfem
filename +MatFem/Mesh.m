classdef Mesh < handle
	% MESH constructs a mesh object.
	% P: node vertices, each row is (x,y) coordinates of a vertex.
	% T: connectivity matrix, each row defines an element.
	
	properties (SetAccess = private)
		P
		T
	end
	
	properties (Dependent)
		nv
		ne
		nlv
	end
	
	methods
		function obj = Mesh(P, T)
			% RectMesh constructs a rectangular mesh.
			%	domain = [left, right, bottom, top] defines the boundary of
			%	the mesh,
			%	n = [nx, ny] defines the number of nodes in each direction.
			obj.P = P;
			obj.T = T;
		end
		
		function r = get.nv(obj)
			r = size(obj.P, 1);
		end
		
		function r = get.ne(obj)
			r = size(obj.T, 1);
		end
		
		function r = get.nlv(obj)
			r = size(obj.T, 2);
		end
		
		function [idx, vtc] = elem(obj, i)
			idx = obj.T(i, :);
			vtc = obj.P(idx, :);
		end
		
		function plot(obj, varargin)
			patch('Faces', obj.T, 'Vertices', obj.P, 'FaceColor', ...
				  'White', varargin{:});
		end
	end
end