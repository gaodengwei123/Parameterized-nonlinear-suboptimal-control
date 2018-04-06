% build by dengwei gao 2018.1.16
% this is build to solve LOS rel motion spacecraft
function [sol,uf] = solve_H_PP(sys,x,xG,Ax,Bx,Q,R,delta_x,ts,order,gamma)
num = length(x);
clean_tol = 1e-4;
g = sym('g',[1 num],'real')';

% set coe
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

A0 = A(xG);
B0 = B;
S0 = solveP(A0,B0,C,D,R,B2);
%%
% Ag = subs(Ax-A0,x,g);
E1 = zeros(6);E1(2,5) = 1;              % 1/x1-A0(2,5);
E2 = zeros(6);E2(3,6) = 1;              % 1/x1/cos(x2)-A0(3,6)
E3 = zeros(6);E3(4,5) = 1;E3(5,4) = -1; % x5/x1-A0(4,5)
E4 = zeros(6);E4(4,6) = 1;E4(6,4) = -1; % x6/x1-A0(4,6)
E5 = zeros(6);E5(5,6) = -1;E5(6,5) = 1; % (x6*sin(x2))/(x1*cos(x2))-A0(6,5)
Ac = (A0-(B2*R^-1*B2')*S0)';%-gamma^-2*B1*B1'

S1 = TaylorLQR(S0,Ac,B0,R,1/g(1)-A0(2,5),E1,order);
S2 = TaylorLQR(S0,Ac,B0,R,1/g(1)/cos(g(2))-A0(3,6),E2,order);
S3 = TaylorLQR(S0,Ac,B0,R,g(5)/g(1)-A0(4,5),E3,order);
S4 = TaylorLQR(S0,Ac,B0,R,g(6)/g(1)-A0(4,6),E4,order);
S5 = TaylorLQR(S0,Ac,B0,R,(g(6)*sin(g(2)))/(g(1)*cos(g(2)))-A0(6,5),E5,order);
S = S0+S1+S2+S3+S4+S5;

% meke the S as ploynominal
% p_x = msspoly('x',num);
% for i=1:num
%     for j=1:num
%         [Ss, ~, inargs] = plfcnchk(char(S(i,j)));
% %         sym_inargs = sym(inargs);
% %         [~,a2]=ismember(sym_inargs,g);
%         Sn(i,j) = build_poly(@(b)Ss(b(1),b(2),b(5),b(6)),xG,5,p_x,clean_tol);
%     end
% end

Kk = R^-1*B2';
uf = @(t,b)double(-Kk*subs(S,g,b)*delta_x(t,b));

handlesys = @(t,x)(sys.dynamics(t,x,uf(t,x)));
sol = ode45(handlesys,[ts(1),ts(end)],xG);
 
end

function p=build_poly(fun,x0,order,p_x,clean_tol)
xu=TaylorVar.init(x0,order);
p=getmsspoly(fun(xu),p_x-x0);
p=clean(p,clean_tol);
end
 
function [P, Q, gopt] = solveP(A,B,C,D,R,B2)
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