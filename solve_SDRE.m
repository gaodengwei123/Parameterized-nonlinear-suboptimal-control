% build by dengwei gao 2018.1.5
function [sol,uf] = solve_SDRE(sys,x,xG,Ax,Bx,Q,R,delta_x,ts)

KD = @(a)lqr(double(subs(Ax,x,a)),Bx,Q,R);

uf = @(t,b)double(-KD(b)*delta_x(t,b));
handlesys = @(t,x)(sys.dynamics(t,x,uf(t,x)));

sol = ode45(handlesys,[ts(1),ts(end)],xG);

end

 



