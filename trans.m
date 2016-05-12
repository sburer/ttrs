function [H,h,y,z] = trans(A1,a1,A2,a2,y1,z1)

% [H,Hsqrt,h] = trans(A1,a1,A2,a2) takes in two ellipsoids and translates
% the first one E1 = (A1,a1) affinely to the unit ball (I,[0;0]). The 
% second ellipsoid E2 = (A2,a2) is translated to (H,h) under the same
% affine translation.

if nargin == 4
    n = length(a1);
    y1 = zeros(n,1);
    z1 = zeros(n,1);
end

[V,D] = eig(A1);
temp = sqrt(D)*V';
H = sqrt(D)\V'*A2*V/sqrt(D);
h = temp*(a2-a1);
y = temp*(y1-a1);
z = temp*(z1-a1);