classdef Benchmark<GaoSystem
    properties
        b = 0.1;
        g=9.8;
        xG = [pi;0];
        uG=0;
        Q;R;uppercontrol
    end
    
    methods
        function obj = Benchmark()
            mark = 0;
            obj = obj@GaoSystem(mark,2,2);
        end
        function [xdot,df]= dynamics(obj,t,x,u)
            xdot = [x(1)-x(1)^3+x(2)+u(1)
                    x(1)+x(1)^2*x(2)-x(2)+u(2)];
            if (nargout>1) 
                % Compute Gradients:
                df = sparse(2,5);
                df(1,2) = 1-3*x(1)^2;
                df(1,3) = 1;
                df(2,2) = 1+2*x(1)*x(2);
                df(2,3) = x(1)^2-1;
                df(1,4) = 1;
                df(2,5) = 1;
            end 
        end
        
        % systme output
        function y = output(obj,t,x,u)
            y = x;
        end
        
        function A = dynamicA(obj, t, x, u)
            [~,df] = obj.dynamics(t,x,u);
            nX = obj.getNumStates;
            A = full(df(:,1+(1:nX)));
        end
        
        function B = dynamicB(obj, t, x, u)
            [~,df] = obj.dynamics(t,x,u);
            nX = obj.getNumStates;
            nU = obj.getNumInput;
            B = full(df(:,nX+1+(1:nU)));
        end
        
        function xdot= adp_dynamics(obj,t,x,u,f,F,pV)
%             xdot = [x(1)-x(1)^3+x(2)+u(1)-x(3)
%                     x(1)+x(1)^2*x(2)-x(2)+u(2)-x(3)
%                     pV]; 
            xdot =[dynamics(obj,t,x,u)+f*F;pV];
        end
        
    end
end