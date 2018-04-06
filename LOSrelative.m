classdef LOSrelative<GaoSystem
    properties
        Bx = [zeros(3);[-1 0 0;0 -1 0;0 0 1]];
        miu = 3.9860044e14;
    end
    
    methods
        function obj = LOSrelative()
            mark = 0;
            obj = obj@GaoSystem(mark,6,3);
        end
        
        function xdot= dynamics(obj,t,x,u)
            xdot = [x(4)
                x(5)/x(1)
                x(6)/(x(1)*cos(x(2)))
                (x(5)^2+x(6)^2)/x(1)
                -x(4)*x(5)/x(1)-x(6)^2*sin(x(2))/(x(1)*cos(x(2)))
                -x(4)*x(6)/x(1)+x(5)*x(6)*sin(x(2))/(x(1)*cos(x(2)))]...
                +obj.Bx*u;
        end
        
        function xcdot= MOVEdynamics(obj,t,x,u,Bw,f,PV)
            xcdot = [dynamics(obj,t,x,u)+Bw*f;PV];
        end
        
        function dq = attitude(obj,t,q)
            % calculate the Attitude dynamic and kinematics
            % target interia: diag([4000 5000 1000])
            dq = [(-q(2)*q(5)-q(3)*q(6)-q(4)*q(7))/2
                (q(1)*q(5)-q(4)*q(6)+q(3)*q(7))/2
                (q(4)*q(5)+q(1)*q(6)-q(2)*q(7))/2
                (-q(3)*q(5)+q(2)*q(6)+q(1)*q(7))/2
                q(7)*q(6)
                -3/5*q(5)*q(7)
                -q(6)*q(5)];
        end
        
        function output = xi_initial(obj,rt,rc,drt,drc)
            %initial state in LOS
            rct=rt-rc;
            drct=drt-drc;
            x=rct(1);
            y=rct(2);
            z=rct(3);
            dx=drct(1);
            dy=drct(2);
            dz=drct(3);
            
            p=sqrt(x^2+y^2+z^2);
            qE=atan(y/sqrt(x^2+z^2));
            
            qB=atan2(-z,x);
            
            dp=(x*dx+y*dy+z*dz)/p;
            dqE=(dy*sqrt(x^2+z^2)-(x*y*dx+y*z*dz)/sqrt(x^2+z^2))/(p^2);
            dqB=(z*dx-x*dz)/(x^2+z^2);
            x1=p;
            x2=qE;
            x3=qB;
            x4=dp;
            x5=p*dqE;
            x6=p*cos(qE)*dqB;
            output = [x1; x2; x3; x4; x5; x6];
        end
        
        function qef = desire_qe(obj,q,nb)
            [x,y,z,~,~,~] = obj.desire_dqefdqbf(q,nb);
            qef = atan(y/sqrt(x^2+z^2)); % -pi/2~pi/2
        end
        function qbf = desire_qb(obj,q,nb)
            [x,~,z,~,~,~] = obj.desire_dqefdqbf(q,nb);% -pi~pi
            qbf = atan2(-z,x);
        end
        function dqef = desire_dqe(obj,q,nb)
            [x,y,z,dx,dy,dz] = obj.desire_dqefdqbf(q,nb);
            if((x^2+z^2) > eps)
                dqef=(dy*sqrt(x^2+z^2)-(x*y*dx+y*z*dz)/sqrt(x^2+z^2))/(x^2+y^2+z^2);
            else
                dqef=(dz*sin(qbf)-dx*cos(qbf))/(x^2+y^2+z^2);
            end
        end
        function dqbf = desire_dqb(obj,q,nb)
            [x,~,z,dx,~,dz] = obj.desire_dqefdqbf(q,nb);
            if((x^2+z^2) > eps)
                dqbf=(z*dx-x*dz)/(x^2+z^2);
            else
                dqbf = 0;
            end
        end
        
        function [x,y,z,dx,dy,dz] = desire_dqefdqbf(obj,q,nb)
            % desire position in LOS
            q1 = q(1);q2 = q(2);q3 = q(3);q4 = q(4); wb=q(5:7);wb=wb(:);
            %             global qbf_i
            %             global qbf0
            % Cross-product matrix
            mys_x = @(qwe1,qwe2,qwe3)([0 -qwe3 qwe2;qwe3 0 -qwe1;-qwe2 qwe1 0]);
            q = [q1;q2;q3];   % quaternion of target
            Cbi = (q4^2-q'*q)*eye(3)+2*q*q'-2*q4*mys_x(q1,q2,q3);% I->b
            wbi = Cbi'*wb;    % body angle vel project in I
            ni = Cbi'*(-nb);  % desire pointing in I
            x = ni(1);
            y = ni(2);
            z = ni(3);
            dd = mys_x(wbi(1),wbi(2),wbi(3))*ni;
            dx = dd(1);
            dy = dd(2);
            dz = dd(3);
        end
        
        function coe = coe_from_sv(obj,R,V)
            mu = obj.miu;
            eps=1e-10;
            r=norm(R);
            v=norm(V);
            vr=dot(R,V)/r;
            H=cross(R,V);
            h=norm(H);
            incl=acos(H(3)/h);
            N=cross([0 0 1],H);
            n=norm(N);
            % Right ascension of ascending node
            if n~=0
                RA=acos(N(1)/n);
                if N(2)<0
                    RA=2*pi-RA;
                end
            else
                RA=0;
            end
            % Eccentricity
            E=1/mu*((v^2-mu/r)*R-r*vr*V);
            e=norm(E);
            % Argument of perigee
            if n~=0
                if e>eps
                    w=acos(dot(N,E)/n/e);
                    if E(3)<0
                        w=2*pi-w;
                    end
                else
                    w=0;
                end
            else
                w=0;
            end
            % Initial true anomaly
            if e>eps
                TA=acos(dot(E,R)/e/r);
                if vr<0
                    TA=2*pi-TA;
                end
            else
                cp=cross(N,R);
                if cp(3)>=0
                    TA=acos(dot(N,R)/n/r);
                else
                    TA=2*pi-acos(dot(N,R)/n/r);
                end
            end
            %             Semimajor axis
            a=h^2/mu/(1-e^2);
            %             Orbit parameters
            coe=[h e RA incl w TA a];
        end
        
    end
end