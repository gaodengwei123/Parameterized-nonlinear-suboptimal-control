classdef Subopti<GaoSystem
    properties
        b = 0.1;
        g=9.8;
        xG = [pi;0];
        uG=0;
        Q;R;uppercontrol
    end
    
    methods
        function obj = Subopti()
            mark = 0;
            obj = obj@GaoSystem(mark,2,1);
        end
        function xdot= dynamics(obj,t,x,u)
            xdot = [x(2)+(x(2)-1)*u
                    x(1)+x(2)-(x(1)+2)*u];
        end
        
        % systme output
        function y = output(obj,t,x,u)
            y = x;
        end 
    end
end