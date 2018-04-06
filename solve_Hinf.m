function [sol,uf] = solve_Hinf(sys,x,xG,Ax,Bx,Q,R,delta_x,ts)
% not complete amybe there are mistake
A = @(a)double(subs(Ax,x,a));

B1 = -Bx;    % w
B2 = Bx;     % u
B = [B1 B2];

% [X,L,G] = care(double(subs(Ax,x,xG)),B2,Q,R);
C1 = sqrt(Q);
C2 = eye(6);
C = [C1;C2];

D11 = zeros(6,3);
D12 = [zeros(3,3);sqrt(R)];
D21 = zeros(6,3);
D22 = zeros(6,3);
D = [D11 D12;D21 D22];

KD = @(a)solveK(A(a),B,C,D,R,B2);
% [K, Q, gopt] = solveK(A(a),B,C,D,R,B2);
% KD = @(a)lqr(double(subs(Ax,x,a)),Bx,Q,R);
f=[0;0;0];
uf = @(t,b)double(-KD(b)*delta_x(t,b));
handlesys = @(t,x)(sys.MOVEdynamics(t,x,uf(t,x),B1,f));
options = odeset('RelTol',1e-4);
sol = ode45(handlesys,[ts(1),ts(end)],xG,options);
% drawx = deval(sol,ts);
% for i=1:length(ts)
%     control(:,i) = uf(ts(i),drawx(:,i));
% end


end


function [K, Q, gopt] = solveK(A,B,C,D,R,B2)
global savegopt
r=[6,3];
G = 0;
gmin=0;
tol=1;

P = ltisys(A, B, C, D);
options = [0,0,0,1];
for i=1:10
    [gopt, ~, x1,x2, y1, y2] = hinflmi(P, r, G, tol,options);
    % output = HinfSDP(A(xG),B,C,D,r,gmin,tol);
    %  gopt = output.gamma;
    if gopt-gmin<tol % convergence
        break;
    end
    gmin = gopt;
    G = gopt+100;
end
P = x2/x1;
Q = y2/y1;
% P=output.Xinf;
K=R^-1*B2'*P;
savegopt = gopt;
end


