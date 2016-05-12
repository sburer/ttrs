% [lambda,v] = lambda_ball(H,h,y,z,v)
% This function gives a lower bound of lambda in the lifted RLT constraint.
% y on the boundary of the unit ball, z on the boundary of (H,h)
% alpha: a bound of Ty*Tz/Lyz2

function [lambda,v] = lambda_ball(H,h,y,z,v)

% Connect y,z and extend; find w in the unit ball
syms tt
st = solve( (y+tt*(z-y))'*(y+tt*(z-y)) == 1 );
t = double(st(1));

if abs(t) < 1e-7
    t = double(st(2));
end
w = y + t*(z-y);

% Bisiection on alpha
alpha0 = 0;
alpha1 = 1;
eps = 1e-7;
n = length(h);
%P = eye(n) - ((z-y)*(z-y)')/(norm(z-y)^2);

% %% The idea of decomposition
P = v*v';

while(alpha1 - alpha0 > eps)
    alpha = (alpha0 + alpha1)/2;
    Q = y*w' - alpha*P;
    c = alpha*P*y - 1/2*(y+w);
    rho = 1 - alpha*y'*P*y;
    u = TRS(Q,c) + rho;
    if u  < 0
        alpha1 = alpha;
    else
        alpha0 = alpha;
    end
end

% Check whether alpha == 0
if alpha < 1e-4
    error('Alpha is zero.');
end

% Bisection on lambda
s0 = 0;
Tw = @(x) 1 - w'*x;
Tz = @(x) 1 - (z-h)'*H*(x-h);
s1 = Tz(y)/Tw(y)/alpha;
if s1 < eps
    warning('Numerical error.');
    fprintf('\nTz(y) = %f',Tz(y));
    v = 'flag';
    lambda = 0;
    return;    
end
[V,D] = eig(H); 
Hsqrt = V*sqrt(D)*V';
l_0 = 1 + (z-h)'*H*h;

while ( s1 - s0 > eps )
    s = (s0 + s1)/2;
    l0 = l_0 - s;
    l = H*(z-h) - s*w;
    [r,dx,dX] = TRS1(eye(n),zeros(n,1),H,Hsqrt,h,l0,l);
    
    if r < 1e-4
        s1 = s;
    else
        s0 = s;
    end
end
lambda = alpha*s; 

% Make a perturbation if necessary
i = 0;
while (max(max(abs(dx*dx' - dX))) > 1e-3) && (i < 10)
    Rad = 1e-3*(rand(n,1)-.5);
    [~,dx,dX] = TRS1(eye(n),Rad,H,Hsqrt,h,l0,l);
    i = i+1;
end
if i == 10
    error('Solution to TRS1 is not well recovered.');
end
v = dx;