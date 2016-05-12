% lambda = lambda_generalhd(A1,a1,A2,a2,y,z)
% This function returns the lambda in the lifted RLT constraint given y in
% E1 and z in E2. E1 and E2 are given by (A1,a1) and (A2,a2), respectively.

% The function calls: trans.m and lambda_ball.m (which calls TRS1.m)

function lambda = lambda_generalhd(A1,a1,A2,a2,y,z,vv)

% Translate E1 to the unit ball 
[H1,h1,y1,z1] = trans(A1,a1,A2,a2,y,z);
[V1,D1] = eig(A1); 
vv1 = (V1*sqrt(D1))\vv;
vv1 = vv1/norm(vv1);

% Find where the optimal lambda is obtained
[~,v1] = lambda_ball(H1,h1,y1,z1,vv1);
if strcmp(v1,'flag')
    lambda = 0;
    return
end
v1 = V1*(sqrt(D1)\v1) + a1;

% Translate E2 to the unit ball 
[H2,h2,z2,y2] = trans(A2,a2,A1,a1,z,y);
[V2,D2] = eig(A2);
vv2 = (V2*sqrt(D2))\vv;
vv2 = vv2/norm(vv2);

% Find where the optimal lambda is obtained
[~,v2] = lambda_ball(H2,h2,z2,y2,vv2);
if strcmp(v2,'flag')
    lambda = 0;
    return
end
v2 = V2*(sqrt(D2)\v2) + a2;

% Calculate corresponding lambda's and compare
Ty = @(x) 1 - (y-a1)'*A1*(x-a1);
Tz = @(x) 1 - (z-a2)'*A2*(x-a2);
%n = length(a1);
%P = eye(n) - ((z-y)*(z-y)')/(norm(z-y)^2);
P = vv*vv';

Lyz2 = @(x) (x-y)'*P*(x-y);
lambda1 = Ty(v1)*Tz(v1)/Lyz2(v1);
lambda2 = Ty(v2)*Tz(v2)/Lyz2(v2);
lambda = min(lambda1,lambda2);