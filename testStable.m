clear
% close all
clc
dbstop if error
% figure
% hold on
% axHandle = gca;
% xlabel(axHandle,'x1');
% ylabel(axHandle,'x2');
% [x1,x2] = meshgrid(-5:0.5:5,-5:0.5:5);
% u = x1-x1.^3+x2;
% v = x1+x2.*x1.^2-x2;
% sumuv = sqrt(u.^2+v.^2);
% u = u./sumuv;
% v = v./sumuv;
% quiver(x1,x2,u,v,0.5)

sys = Benchmark();
% set A0 and parmeterized
xG = zeros(2,1);uG = zeros(2,1);
g = msspoly('g',2);
x = msspoly('x',2);
Ax = [1-x(1)^2   1;1 x(1)^2-1];
A0 = double(subs(Ax,x,xG));
B0 = eye(2);

% SDRE control
Q = eye(2); R=2*eye(2);
KD = @(a)lqr(double(subs(Ax,x,a)),B0,Q,R);
uf =@(b)double(-KD(b)*b);
handlesys = @(t,x)(sys.dynamics(t,x,uf(x)));
sol = ode45(handlesys,[0,5],[1;1]);
ts = 0:0.1:5;
drawx = deval(sol,ts);
% plot
% plot(ts,drawx(1,:),'y','lineWidth',2)
% plot(ts,drawx(2,:),'k','lineWidth',2)
% plot control
% for i=1:length(ts)
%     uu(:,i) = uf(drawx(:,i));
% end
% P1 = plot(ts,uu(1,:),'k','lineWidth',2);
% plot(ts,uu(2,:),'k','lineWidth',2)


%% algebraic control
[K,S0,~] = lqr(A0,B0,Q,R);
% Ag = subs(Ax-A0,x,g);
E = [-1 0;0 1];
Ac = (A0-B0*R^-1*B0'*S0)';
% axis off
xs = sym('x',[1 2],'real')';
for i=1:4
    order = 3;
    S1 = TaylorLQR(S0,Ac,B0,R,g(1)^2,E,order);
    S = S0+S1;

    Ss = msspoly2sym(g,xs,S);
    V = xs'*Ss*xs;
%     V = xs'*S*xs;
    Vs = plfcnchk(char(V));
    [x1,x2] = meshgrid(-5:0.5:5,-5:0.5:5);
    valueV = feval(Vs,x1,x2);
    figure
    hold on
    pv=surf(x1,x2,valueV,'edgecolor','none');
    shading interp;
    surf(x1,x2,0*ones(size(x1)),'edgecolor','none','facealpha',0.5,'FaceColor',[0 0 0]);
    axHandle = gca;
    xlabel(axHandle,'x1');
    ylabel(axHandle,'x2');
    zlabel(axHandle,'value V');
    
    uf = clean(-R^-1*B0'*S*g);
    us = msspoly2sym(g,xs,uf);
    dotxs = sys.dynamics(0,xs,us);
    dotV = [diff(V,xs(1)) diff(V,xs(2))]*dotxs;
    dotVs = plfcnchk(char(dotV));
    valuedotV = feval(dotVs,x1,x2);
    figure
    hold on
    pvdot=surf(x1,x2,valuedotV,'edgecolor','none');
    shading interp;
    surf(x1,x2,eps*ones(size(x1)),'edgecolor','none','facealpha',0.5,'FaceColor',[0 0 0]);
    axHandle = gca;
    xlabel(axHandle,'x1');
    ylabel(axHandle,'x2');
    zlabel(axHandle,'$\dot{V}$','interpreter','latex');
    
    uf = subs(uf,g,x);
    uf = @(a)double(subs(uf,x,a));
    
    handlesys = @(t,x)(sys.dynamics(t,x,uf(x)));
    sol = ode45(handlesys,[0,5],[1;1]);
    ts = 0:0.01:5; 
    drawx = deval(sol,ts);
    for j=1:length(ts)
        uu(:,j) = uf(drawx(:,j));
    end
    switch order
        %         case 0
        %             plot(ts,drawx(1,:),'r--')
        %             plot(ts,drawx(2,:),'b--')
        case 1
            P2 = plot(ts,uu(1,:),'b--');
            plot(ts,uu(2,:),'b--')
            % plot(ts,drawx(1,:),'b--')
        case 2
            P3 = plot(ts,uu(1,:),'m--');
            plot(ts,uu(2,:),'m--')
            % plot(ts,drawx(1,:),'m--')
        case 3
            P4 = plot(ts,uu(1,:),'g--');
            plot(ts,uu(2,:),'g--')
            % plot(ts,drawx(1,:),'g--')
        case 4
            P5 = plot(ts,uu(1,:),'r--');
            plot(ts,uu(2,:),'r--')
            % plot(ts,drawx(1,:),'r--')
    end
end
legend([P1,P2,P3,P4,P5],'SDRE','1-order','2-order','3-order','4-order')
% ingraph = axes('Position',[0.45 0.4 0.4 0.2]);
% set(ingraph,'xtick',[],'ytick',[])
% set(gca,'xtick',[],'ytick',[])





