function [dy,dx,dX] = TRS1(A,a,H,H_half,h,l0,l)

% The function solves the following minimization problem via SDP relaxation
% min 1-(x-a)'*A*(x-a)
% S.t. ||H_half*(x-h)|| <= 1
%      l0 - l'*x   <= 0
% where H_half is a n-by-n symmetric matrix and H = H_half^2.
% Outputs: dx,dX--solution; dy--minimum value

myset = sdpsettings('verbose',0);

n = length(a);

X = sdpvar(n,n,'symm');
x = sdpvar(n,1);
Y = [1,x';x,X];

obj = -A(:)'*X(:) + 2*a'*A*x;

con = [ Y >= 0 ];
con = [ con ; H(:)'*X(:) - 2*h'*H*x + h'*H*h <= 1 ];
con = [ con ; norm( H_half*(X*l - l0*x - l'*x*h + l0*h) ) <= l'*x - l0 ];

solvesdp(con,obj,myset);
dx = double(x);
dX = double(X);
dy = 1 - A(:)'*dX(:) + 2*a'*A*dx - a'*A*a;

%% NEED A WAY TO RECOVER THE SOLUTION!!!!!!
% if max(max(abs(dx*dx' - dX))) > 1e-3
%     error('Solution to TRS1 is not well recovered.');
% else 
%     v1 = sqrt(diag(dX));
%     v2 = [v1(1);-v1(2)];
%     v3 = -v2;
%     v4 = -v1;
%     if (H(:)'*dX(:) - 2*h'*H*v1 + h'*H*h <= 1+1e-4) && (l'*v1 - l0 <=1e-4)
%         v = v1;
%     elseif (H(:)'*dX(:) - 2*h'*H*v2 + h'*H*h <= 1+1e-4) && (l'*v2 - l0 <=1e-4)
%         v = v2;
%     elseif (H(:)'*dX(:) - 2*h'*H*v3 + h'*H*h <= 1+1e-4) && (l'*v3 - l0 <=1e-4)
%         v = v3;
%     else
%         v = v4;
%     end
% end
