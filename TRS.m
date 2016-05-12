% val = TRS(Q,c)
% The function solves the standard TRS (as below) via SDP.
% min. x'*Q*x + 2*c'*x
% S.t. x'*x <= 1

function val = TRS(Q,c)
n = length(c);
myset = sdpsettings('verbose',0);

X = sdpvar(n,n,'symm');
x = sdpvar(n,1);
Y = [1,x';x,X];

obj = Q(:)'*X(:) + 2*c'*x;

con = [ Y >= 0 ];
con = [ con ; trace(X) <= 1 ];

solvesdp(con,obj,myset);
% dx = double(x)
% dX = double(X)
% eig(double(Y))
val = double(obj);