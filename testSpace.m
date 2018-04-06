clear
close all
clc
dbstop if error
deg = pi/180;
sys = LOSrelative();
sys = setInputLimits(sys,5,5);
nb=[1;0;0];
rouf=20;
drouf=0;
% initial state
rt0=[-6713501.42;2108702.64;3864206.20];
rc0=[-6713602.78;2108423.08;3864165.94];
drt0=[-3796.728;-4213.883;-4159.899];
drc0=[-3796.586;-4213.906;-4160.010];
% coe1=sys.coe_from_sv(rt0,drt0);
% coe2=sys.coe_from_sv(rc0,drc0);
% initial state
xr = sys.xi_initial(rt0,rc0,drt0,drc0);
xG = xr(:);uG = zeros(3,1);

%% desire state
% target quasi and angle-vel
wb = [0.2;0.2;0.2]*deg;
qt0 = [0.7071,-0.31,0.55,-0.32]';
handlesys_att = @(t,x)(sys.attitude(t,x));
sol = ode45(handlesys_att,[0,100],[qt0;wb]);
% solve desire trajectory
ts = 0:0.1:100;
qbtwb = @(t)deval(sol,t);
for i=1:length(ts)
    qef(i) = sys.desire_qe(qbtwb(ts(i)),nb);
    qbf(i) = sys.desire_qb(qbtwb(ts(i)),nb);
    dqef(i) = sys.desire_dqe(qbtwb(ts(i)),nb);
    dqbf(i) = sys.desire_dqb(qbtwb(ts(i)),nb);
end
fun1 = spline(ts,qef');
fun2 = spline(ts,qbf');
fun3 = spline(ts,dqef');
fun4 = spline(ts,dqbf');
qef = polyniminalTrajectory(fun1);
qbf = polyniminalTrajectory(fun2);
dqef = polyniminalTrajectory(fun3);
dqbf = polyniminalTrajectory(fun4);
% figure
% hold on
% plot(qef,'r')
% plot(qbf,'k')


%% control solve
% set A0 and parmeterized
x = sym('x',[1 6],'real')';

Axx = [1            0           0
    0            1/x(1)       0
    0            0           1/x(1)/cos(x(2))
    0            x(5)/x(1)   x(6)/x(1)
    -x(5)/x(1)    0           -x(6)/x(1)/cos(x(2))*sin(x(2))
    -x(6)/x(1)    x(6)/x(1)/cos(x(2))*sin(x(2)) 0];
Ax = [zeros(6,3),Axx];
Bx = sys.Bx;

% error state
% delta_x = @(t,b)[b(1)-rouf
%     b(2)-qef.eval(t)
%     b(3)-qbf.eval(t)
%     b(4)-drouf
%     b(5)-b(1)*dqef.eval(t)
%     b(6)-b(1)*dqbf.eval(t)*cos(b(2))];
delta_x = @(t,b)[b(1)-rouf
    b(2)-qef.eval(t)
    b(3)-qbf.eval(t)
    b(4)-drouf
    b(5)-b(1)*dqef.eval(t)
    b(6)-b(1)*dqbf.eval(t)*cos(qef.eval(t))];
% control paremters
Q = 0.01*diag([1 100000 10000 100 100000 100000]);
R = 1000*diag([1 10 0.1]);
% Q = diag([1 1000 50 50 20 50]);
% R = 50*diag([100 0.05 0.5]);
% start SDRE control ======================================================
% [sol_SDRE,uf] = solve_SDRE(sys,x,xG,Ax,Bx,Q,R,delta_x,ts);
% drawx = deval(sol_SDRE,ts);
% algebraic control  ======================================================
order = 1;
[sol_PP,ufPP] = solve_PPcontrol(sys,x,xG,Ax,Bx,Q,R,delta_x,ts,order);
drawxPP = deval(sol_PP,ts);
% theta-D control    ======================================================
% [sol_DPP,ufD] = solve_thetaD(sys,x,xG,Ax,Bx,Q,R,delta_x,ts,order);
% drawxD = deval(sol_DPP,ts);
%% plot
figure(1)
xlabel('times(s)');
title('state response')
subplot(3,1,1)
ylabel('\rho (m)');
hold on
plot(ts,drawx(1,:),'k','lineWidth',2)
plot(ts,drawxPP(1,:),'r--')
plot(ts,drawxD(1,:),'b:')
% plot(ts,repmat(rouf,1,length(ts)),'k','lineWidth',2)
subplot(3,1,2)
ylabel('q_\epsilon  (rad)');
hold on
plot(ts,drawx(2,:),'k','lineWidth',2)
plot(ts,drawxPP(2,:),'r--')
plot(ts,drawxD(2,:),'b:')
% qef = rebreaks(qef,ts);
% hh=plot(qef,'k');set(hh,'lineWidth',2)
subplot(3,1,3)
ylabel('q_\beta (rad)');
hold on
plot(ts,drawx(3,:),'k','lineWidth',2)
plot(ts,drawxPP(3,:),'r--')
plot(ts,drawxD(3,:),'b:')
% qef = rebreaks(qef,ts);
% hh=plot(qbf,'k');set(hh,'lineWidth',2)
figure(2)
grid on
hold on
Traj3D = @(x1,x2,x3)[cos(x2)*cos(x3),sin(x2),-cos(x2)*sin(x3);-sin(x2)*cos(x3),cos(x2),sin(x2)*sin(x3);sin(x3),0,cos(x3)]*[-x1;0;0];
for i=1:length(ts)
    Traj3D_SDRE(:,i) = Traj3D(drawx(1,i),drawx(2,i),drawx(3,i));
    Traj3D_PP(:,i) = Traj3D(drawxPP(1,i),drawxPP(2,i),drawxPP(3,i));
    Traj3D_D(:,i) = Traj3D(drawxD(1,i),drawxD(2,i),drawxD(3,i));
end
plot3(Traj3D_SDRE(1,:),Traj3D_SDRE(2,:),Traj3D_SDRE(3,:),'k','lineWidth',2)
plot3(Traj3D_PP(1,:),Traj3D_PP(2,:),Traj3D_PP(3,:),'r--')
plot3(Traj3D_D(1,:),Traj3D_D(2,:),Traj3D_D(3,:),'b:')
legend('SDRE','our algorithm','\theta-D')
% minpage
% fig=figure;
% hold on
% plot3(Traj3D_SDRE(1,500:end),Traj3D_SDRE(2,500:end),Traj3D_SDRE(3,500:end),'k','lineWidth',2)
% plot3(Traj3D_PP(1,500:end),Traj3D_PP(2,500:end),Traj3D_PP(3,500:end),'r--')
% plot3(Traj3D_D(1,500:end),Traj3D_D(2,500:end),Traj3D_D(3,500:end),'b:')
 

 
figure(3)
title('control response')
tic
for i =1:length(ts)
    t = ts(i);
    uu(:,i) = uf(t,drawx(:,i));
end
time1=toc
tic
for i =1:length(ts)
    t = ts(i);
    uuPP(:,i) = ufPP(t,drawxPP(:,i));
end
time2=toc
tic
for i =1:length(ts)
    t = ts(i);
    uuD(:,i) = ufD(t,drawxD(:,i));
end
time3=toc
subplot(3,1,1)
ylabel('ux (m/s^2)');
hold on
plot(ts,uu(1,:),'k','lineWidth',2);
plot(ts,uuPP(1,:),'r--');
plot(ts,uuD(1,:),'b:')
subplot(3,1,2)
ylabel('uy (m/s^2)');
hold on
plot(ts,uu(2,:),'k','lineWidth',2)
plot(ts,uuPP(2,:),'r--');
plot(ts,uuD(2,:),'b:')
subplot(3,1,3)
ylabel('uz (m/s^2)');
xlabel('times (s)');
hold on
plot(ts,uu(3,:),'k','lineWidth',2)
plot(ts,uuPP(3,:),'r--');
plot(ts,uuD(3,:),'b:')
% legend('SDRE','our algorithm','\theta-D')

  