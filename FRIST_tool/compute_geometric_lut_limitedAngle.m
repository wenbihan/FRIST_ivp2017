function [L,theta] = compute_geometric_lut_limitedAngle(w,stepAngle)

% compute_geometric_lut - compute geometric permutation
%
%   L = compute_geometric_lut(w);
%
%   w is the size of the patch.
%   L(:,i) is the ith permutation.
%   
%   For a (w,w) array M, the 1D signal corresdponding 
%       to this oriented flattening is M(L(:,i));
%   
%   Copyright (c) 2006 Gabriel Peyr?%   
%   Revised by Xiaobo Qu, July 26, 2011
%   Add the stepAngle to simplify the number of points on grids.


% compute the rotation angles
% [Y,X] = meshgrid(0:w-1, 0:w-1); 
if nargin<2
   stepAngle=1; 
end
[Y,X] = meshgrid(0:stepAngle:w-1, 0:stepAngle:w-1); 


X = X(:); Y = Y(:);
X(1) = []; Y(1) = [];
theta = atan2(Y(:),X(:));
theta = unique(theta);
theta = [-theta(end-1:-1:2); theta]';
% take mid points
theta = ( theta + [theta(2:end),theta(1)+pi] )/2;
m = length(theta); % number of directions

% perform rotation of the points
[Y,X] = meshgrid(0:w-1, 0:w-1); X = X(:); Y = Y(:);
X = repmat(X, [1 m]);
Y = repmat(Y, [1 m]);
Theta = repmat(theta, [w^2 1]);
D = -sin(Theta).*X + cos(Theta).*Y;

% the LUT correspond to the 1D ordering of the projections
[~,L] = sort(D,1);
%% added by Xiaobo Qu
%% compatiable for C language
% L=L-1;
%%