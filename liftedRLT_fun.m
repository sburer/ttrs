% Function to implement the lifted-RLT method. Input is the
% filename of the instance mat file from reference [8]. See
% https://github.com/sburer/soctrust
%
% rm_old = "rank measure" from the method of reference [8]
%
% rg = relative gap between final relaxation bound and primal solution
%
% rm_new = rank measure from the lifted-RLT method
%
% cuts = number of lifted-RLT cuts
%
% rv = final relaxation value

function [rm_old, rg, rm_new, cuts, rv] = liftedRLT_fun(filename)

% Load instance

load(filename);

% Set various parameters

tol = 1.0e-14;
maxiter = 10;
myset = sdpsettings('verbose', 0);

SOCRLT_or_liftedRLT_or_both = 2; % 1, 2, or 3. 2 is the method of the paper

% Initialize from input file
  
n = size(Q,1);
rm_old = log10(eigY(n+1)/eigY(n));

%% Derive the sqrt of H

[V,D] = eig(H); d = diag(D); Hsqrt = V*diag(sqrt(d))*V';

%% Setup first SDP

x = sdpvar(n,1);
X = sdpvar(n,n,'symm');
Y = [1,x';x,X];

con = [ Y >= 0 ];
con = [ con ; trace(X) <= r1^2 ];
con = [ con ; H(:)'*X(:) - 2*h'*H*x + h'*H*h <= r2^2 ];

con_SOCRLT    = [];
con_liftedRLT = [];

obj = Q(:)'*X(:) + c'*x;

%% Setup separation problem (objective TBD)

a = sdpvar(n,1);
A = sdpvar(n,n,'symm');
sepcon = [ [1,a';a,A] >= 0 ];
sepcon = [ sepcon ; trace(A) <= 1 ];

%% Setup and begin iterations...

iter = 1;
stop = 0;
while stop == 0 & iter <= maxiter

    %% Solve current relaxation

    if SOCRLT_or_liftedRLT_or_both == 1
        solvesdp([con;con_SOCRLT],obj,myset);
    elseif SOCRLT_or_liftedRLT_or_both == 2
        solvesdp([con;con_liftedRLT],obj,myset);
    else
        solvesdp([con;con_SOCRLT;con_liftedRLT],obj,myset);
    end

    % Save relaxation value

    rv = double(obj);

    %% Store portion of solution

    dx = double(x);
    dX = double(X);
    eigY = eig(double(Y));

    %% Setup objective for SOC-RLT separation

    tmp1 = r1*(dx-h);
    tmp2 = h*dx' - dX;

    sepobj1 = (r2^2)*(r1^2 - 2*r1*a'*dx + dx'*A*dx);

    sepobj21 = tmp1'*H*tmp1;
    sepobj22 = tmp1'*H*tmp2;
    sepobj23 = tmp2'*H*tmp2;

    sepobj = sepobj1 - (sepobj21 + 2*sepobj22*a + sepobj23(:)'*A(:));

    %% Solve separation subproblem

    solvesdp(sepcon,sepobj,myset);

    %% If no separation, quit

    if double(sepobj)/(norm(Q,'fro') + norm(c)) > tol 
        stop = 1;
    else

    %% Specify new SOCRLT contraint

    da = double(a);
    new_SOCRLT = [ norm( Hsqrt*(r1*(x-h) - X*da + h*da'*x )) <= r2*(r1 - da'*x) ];

    %% If we want liftedRLT constraints...

    if SOCRLT_or_liftedRLT_or_both > 1

        %% Add new_SOCRLT to latest relaxation and solve

        if SOCRLT_or_liftedRLT_or_both == 2
            solvesdp([con;con_liftedRLT;new_SOCRLT],obj,myset);
        else
            solvesdp([con;con_SOCRLT;con_liftedRLT;new_SOCRLT],obj,myset);
        end

        % Extract y and z. Also aa and bb

        y = r1*da;
        y = r1*y/norm(y);
        z = (r1*dx - dX*da)/(r1 - da'*dx);
        z = r2*z/norm(Hsqrt*z);
        aa = y/(r1*norm(y));
        bb = H*z/r2^2;
        % 1 - aa'*y; % Note: This quantity should be 0
        % 1 - bb'*z; % Note: This quantity Should be 0

        % Calculate projection matrix

        dx = double(x);
        dX = double(X);
        P = sdpvar(n);
        newcon = [P >= 0 ; trace(P) == n-1 ; P*(y-z) == 0];
        newobj = P(:)'*dX(:) - 2*y'*P*dx + y'*P*y;
        solvesdp(newcon,-newobj,myset);
        P = double(P);

        % Extract positive eigenvector

        [V,D] = eig(P);
        d = diag(D);
        [~,ind] = max(d);
        if numel(ind) > 1
            eig(P)
            error('Unexpected rank(P) > 1');
        end
        v = V(:,ind);

        % Calculate lambda
        
        lambda = lambda_generalhd(eye(n)/r1^2,zeros(n,1),H/r2^2,h,y,z,v);

        % Save new liftedRLT constraint

        con_liftedRLT = [con_liftedRLT; ...
           1 - aa'*x - bb'*x + aa'*X*bb - lambda*(v'*X*v - 2*y'*v*v'*x + (y'*v)^2) >= 0];

      end

      %% Add new_SOCRLT to con_SOCRLT

      con_SOCRLT = [ con_SOCRLT ; new_SOCRLT ];

    end

    % Increment iteration count
    
    iter = iter + 1;

end

% Calculate final values

dx = double(x);
dX = double(X);
eigY = eig(double(Y));
rg = ((dx'*Q*dx + c'*dx) - double(obj))/abs(dx'*Q*dx + c'*dx);
rm_new = log10(eigY(n+1)/eigY(n));
cuts = iter - 2;
