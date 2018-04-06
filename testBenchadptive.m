clear all
close all
clc
dbstop if error

sys = Benchmark();
% set A0 and parmeterized
xG = zeros(2,1);uG = zeros(2,1);
g = msspoly('g',3);
x = msspoly('x',2);
Ax = [1-x(1)^2   1;1 x(1)^2-1];
A0 = double(subs(Ax,x,xG));
B0 = eye(2);

figure
hold on
% xlabel('times(s)');
% ylabel('x_1');
% title('state response')

xlabel('times(s)');
ylabel('control u');
title('control response')
% SDRE control
Q = eye(2); R=2*eye(2);
KD = @(a)lqr(double(subs(Ax,x,a)),B0,Q,R);
uf =@(b)double(-KD(b)*b);
handlesys = @(t,x)(sys.dynamics(t,x,uf(x)+[1;1]));%
sol = ode45(handlesys,[0,10],[1;1]);
ts = 0:0.1:10;
drawx = deval(sol,ts);
% P1 = plot(ts,drawx(1,:),'k','lineWidth',2)
% plot(ts,0.7768*ones(1,length(ts)),'k--')
% plot control
for i=1:length(ts)
    uu(:,i) = uf(drawx(:,i));
end
P1 = plot(ts,uu(1,:),'k','lineWidth',2);
plot(ts,uu(2,:),'k','lineWidth',2)

%% algebraic control
[K,S0,~] = lqr(A0,B0,Q,R);
% Ag = subs(Ax-A0,x,g);
E = [-1 0;0 1]; % g(1)^2
% F = [1;1];      % g(3)

Ac = (A0-B0*R^-1*B0'*S0)';
% axis off

for i=1:4
    order = 2;
    S1 = TaylorLQR(S0,Ac,B0,R,g(1)^2,E,order);
%     S1 = TaylorLQR(S0,Ac,B0,R,g(1)^2,E,order,g(3),F,F);
%     S2 = TaylorLQR(S0,Ac,B,R,a,E,order,varargin) 
    S = S0+S1;%+S2;
    S = subs(S,g(1:2),x);
    Vs = x'*S*x;
    F = [1;1];
    PdotV = (diff(Vs,x)*F)';
%     PdotVs = diff(Vs,g(3));
    Gamma = 0.1;
    PV = @(a)Gamma*double(subs(PdotV,x,a));
    BL = pinv(B0);
    uf = -R^-1*B0'*S*x -BL*F*g(3);    
    uf = @(a)double(subs(uf,[x;g(3)],a)); 
    handlesys = @(t,a)(sys.adp_dynamics(t,a,uf(a),1,F,PV(a(1:2))));
    sol = ode45(handlesys,[0,10],[1;1;0]); 
    ts = 0:0.01:10; 
    drawx = deval(sol,ts);
    
%     P2=plot(ts,drawx(1,:),'r','lineWidth',2) 
%     plot(ts,zeros(1,length(ts)),'r--')
%     legend([P1,P2],'non-adptive','adptive')
    
    %% plot theta
%     plot(ts,drawx(3,:),'r','lineWidth',2)
%     ylabel('$$\hat{\theta}$$','interpreter','latex');
    for i=1:length(ts)
        uu(:,i) = uf(drawx(:,i));
    end
    P2 = plot(ts,uu(1,:),'r','lineWidth',2);
    plot(ts,uu(2,:),'r','lineWidth',2)
    
    
end 
 





