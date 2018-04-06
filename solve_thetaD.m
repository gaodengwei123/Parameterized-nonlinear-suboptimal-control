% build by dengwei gao 2018.1.5
function [sol,uf] = solve_thetaD(sys,x,xG,Ax,Bx,Q,R,delta_x,ts,order)

A0 = double(subs(Ax,x,xG));
B0 = double(subs(Bx,x,xG));
Kk = R^-1*B0';

[K,T0] = lqr(A0,B0,Q,R);
% set parmeters
ki=0.9999;
li=0.001;

Ac = (A0-B0*Kk*T0)';

uf = @(t,b)usolve(t,x,b,ki,li,T0,Ax,A0,Ac,Kk,B0,delta_x,K,order);
handlesys = @(t,x)(sys.dynamics(t,x,uf(t,x)));
sol = ode45(handlesys,[ts(1),ts(end)],xG);
end


function uf = usolve(t,xval,x,ki,li,T0,Ax,A0,Ac,Kk,B0,delta_x,K,order)
Ax = double(subs(Ax-A0,xval,x));
switch order
    case 0
        uf = -K*delta_x(t,x);
    case 1
        D1 = ki*exp(-li*t)*(-T0*Ax-Ax'*T0);
        FF1 = -T0*Ax-Ax'*T0-D1;
        T1 = lyap(Ac,FF1);
        uf = (-K-Kk*T1)*delta_x(t,x);
    case 2
        D1 = ki*exp(-li*t)*(-T0*Ax-Ax'*T0);
        FF1 = -T0*Ax-Ax'*T0-D1;
        T1 = lyap(Ac,FF1);
        D2 = ki*exp(-li*t)*(-T1*Ax-Ax'*T1+T1*B0*Kk*T1);
        FF2 = -T1*Ax-Ax'*T1+T1*B0*Kk*T1-D2;
        T2 = lyap(Ac,FF2);
        uf = (-K-Kk*T1-Kk*T2)*delta_x(t,x);
    case 3
        D1 = ki*exp(-li*t)*(-T0*Ax-Ax'*T0);
        FF1 = -T0*Ax-Ax'*T0-D1;
        T1 = lyap(Ac,FF1);
        D2 = ki*exp(-li*t)*(-T1*Ax-Ax'*T1+T1*B0*Kk*T1);
        FF2 = -T1*Ax-Ax'*T1+T1*B0*Kk*T1-D2;
        T2 = lyap(Ac,FF2);
        D3 = ki*exp(-li*t)*(-T2*Ax-Ax'*T2 + T1*B0*Kk*T2 + T2*B0*Kk*T1);
        FF3 = -T2*Ax-Ax'*T2 + T1*B0*Kk*T2 + T2*B0*Kk*T1-D3;
        T3 = lyap(Ac,FF3);
        uf = (-K-Kk*T1-Kk*T2-Kk*T3)*delta_x(t,x);
    case 4
        D1 = ki*exp(-li*t)*(-T0*Ax-Ax'*T0);
        FF1 = -T0*Ax-Ax'*T0-D1;
        T1 = lyap(Ac,FF1);
        D2 = ki*exp(-li*t)*(-T1*Ax-Ax'*T1+T1*B0*Kk*T1);
        FF2 = -T1*Ax-Ax'*T1+T1*B0*Kk*T1-D2;
        T2 = lyap(Ac,FF2);
        D3 = ki*exp(-li*t)*(-T2*Ax-Ax'*T2 + T1*B0*Kk*T2 + T2*B0*Kk*T1);
        FF3 = -T2*Ax-Ax'*T2 + T1*B0*Kk*T2 + T2*B0*Kk*T1-D3;
        T3 = lyap(Ac,FF3);
        D4 = ki*exp(-li*t)*(-T3*Ax-Ax'*T3 + T1*B0*Kk*T3 + T2*B0*Kk*T2 + T3*B0*Kk*T1);
        FF4 = -T3*Ax-Ax'*T3 + T1*B0*Kk*T3 + T2*B0*Kk*T2 + T3*B0*Kk*T1-D4;
        T4 = lyap(Ac,FF4);
        uf = (-K-Kk*T1-Kk*T2-Kk*T3-Kk*T4)*delta_x(t,x);
end


end





