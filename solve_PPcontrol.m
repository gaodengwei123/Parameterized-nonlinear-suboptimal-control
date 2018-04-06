% build by dengwei gao 2018.1.5
% this is build to solve LOS rel motion spacecraft
function [sol,uf,Sn] = solve_PPcontrol(sys,x,xG,Ax,Bx,Q,R,delta_x,ts,order)
num = length(x);
clean_tol = 1e-4;
g = sym('g',[1 num],'real')';
% g = msspoly('g',num);
A0 = double(subs(Ax,x,xG));
B0 = double(subs(Bx,x,xG));
[~,S0,~] = lqr(A0,B0,Q,R);
% Ag = subs(Ax-A0,x,g);
E1 = zeros(6);E1(2,5) = 1;              % 1/x1-A0(2,5);
E2 = zeros(6);E2(3,6) = 1;              % 1/x1/cos(x2)-A0(3,6)
E3 = zeros(6);E3(4,5) = 1;E3(5,4) = -1; % x5/x1-A0(4,5)
E4 = zeros(6);E4(4,6) = 1;E4(6,4) = -1; % x6/x1-A0(4,6)
E5 = zeros(6);E5(5,6) = -1;E5(6,5) = 1; % (x6*sin(x2))/(x1*cos(x2))-A0(6,5)
Ac = (A0-B0*R^-1*B0'*S0)';

S1 = TaylorLQR(S0,Ac,B0,R,1/g(1)-A0(2,5),E1,order);
S2 = TaylorLQR(S0,Ac,B0,R,1/g(1)/cos(g(2))-A0(3,6),E2,order);
S3 = TaylorLQR(S0,Ac,B0,R,g(5)/g(1)-A0(4,5),E3,order);
S4 = TaylorLQR(S0,Ac,B0,R,g(6)/g(1)-A0(4,6),E4,order);
S5 = TaylorLQR(S0,Ac,B0,R,(g(6)*sin(g(2)))/(g(1)*cos(g(2)))-A0(6,5),E5,order);
S = S0+S1+S2+S3+S4+S5;

% meke the S as ploynominal
p_x = msspoly('x',num);
for i=1:num
    for j=1:num
        [Ss, ~, inargs] = plfcnchk(char(S(i,j)));
%         sym_inargs = sym(inargs);
%         [~,a2]=ismember(sym_inargs,g);
        Sn(i,j) = build_poly(@(b)Ss(b(1),b(2),b(5),b(6)),xG,5,p_x,clean_tol);
    end
end
% % trans to sys (it's not necessary)
% Ss = msspoly2sym(p_x,x,Sn);
% ufs = @(t,b)double(-Kk*subs(Ss,x,b)*delta_x(t,b));
Kk = R^-1*Bx';
uf = @(t,b)double(-Kk*subs(Sn,p_x,b)*delta_x(t,b));

handlesys = @(t,x)(sys.dynamics(t,x,uf(t,x)));
sol = ode45(handlesys,[ts(1),ts(end)],xG);
 
end

function p=build_poly(fun,x0,order,p_x,clean_tol)
xu=TaylorVar.init(x0,order);
p=getmsspoly(fun(xu),p_x-x0);
p=clean(p,clean_tol);
end
 
 