clear
close all
clc
dbstop if error

sys = Subopti();
% set A0 and parmeterized
xG = zeros(2,1);uG = 0;
g = msspoly('g',2);
x = msspoly('x',2);

Ax = [0  1;  1  1];
Bx = [-1+x(2);-2-x(1)];
A0 = Ax;
B0 = double(subs(Bx,x,xG));

figure
hold on
xlabel('times(s)');
ylabel('x_1,x_2');
% title('state response')

% xlabel('times(s)');
% ylabel('control u');
% title('control response')
% SDRE control
Q = [1 0;0 2]; R=1;
KD = @(a)lqr(A0,double(subs(Bx,x,a)),Q,R);
uf =@(b)double(-KD(b)*b);
% uf=@(b)(b(1)+2*b(2));
handlesys = @(t,x)(sys.dynamics(t,x,uf(x)));
sol = ode45(handlesys,[0,6],[1;1]);
ts = 0:0.01:6;
drawx = deval(sol,ts);
P1=plot(ts,drawx(1,:),'k','lineWidth',2);
plot(ts,drawx(2,:),'k','lineWidth',2)
% plot control
% for i=1:length(ts)
%     uu(i) = uf(drawx(:,i));
% end
% P1 = plot(ts,uu,'k','lineWidth',2); 

%% algebraic control
[K,S0,~] = lqr(A0,B0,Q,R);
% Ag = subs(Ax-A0,x,g);
% E1 = [0 1;-1 0];    %x1*x2;
% E2 = [0 0;0 1];     %x1^2;
F1 = [1;0];           %x2
F2 = [0;-1];          %x1
Ac = (A0-B0*R^-1*B0'*S0)';
% axis off

for i=1:4
    order = i;
    S1 = TaylorLQR(S0,Ac,Bx,R,[],[],order,g(2),F1,B0);
    S2 = TaylorLQR(S0,Ac,Bx,R,[],[],order,g(1),F2,B0);
    S = S0+S1+S2;
    uf = -R^-1*Bx'*S*x;
    uf = subs(uf,g,x);
    uf = @(a)double(subs(uf,x,a));
    
    handlesys = @(t,x)(sys.dynamics(t,x,uf(x)));
    sol = ode45(handlesys,[0,6],[1;1]);
    ts = 0:0.01:6; 
    drawx = deval(sol,ts); 
%     for j=1:length(ts)
%         uu(j) = uf(drawx(:,j));
%     end
    switch order
        %         case 0
        %             plot(ts,drawx(1,:),'r--');
        %             plot(ts,drawx(2,:),'b--');
        case 1
%             P2 = plot(ts,uu,'b--'); 
            P2 = plot(ts,drawx(1,:),'b--');
            plot(ts,drawx(2,:),'b--');
        case 2
%             P3 = plot(ts,uu,'m--'); 
            P3 = plot(ts,drawx(1,:),'m--');
            plot(ts,drawx(2,:),'m--');
        case 3
%             P5 = plot(ts,uu,'r--'); 
            P5 = plot(ts,drawx(1,:),'r--');
            plot(ts,drawx(2,:),'r--');
        case 4
%             P4 = plot(ts,uu,'g--'); 
            P4 = plot(ts,drawx(1,:),'g--');
            plot(ts,drawx(2,:),'g--');
    end
end
legend([P1,P2,P3,P4,P5],'SDRE','1-order','2-order','3-order','4-order')
% ingraph = axes('Position',[0.45 0.4 0.4 0.2]);
% set(ingraph,'xtick',[],'ytick',[])
% set(gca,'xtick',[],'ytick',[])





