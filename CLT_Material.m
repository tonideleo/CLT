classdef CLT_Material < handle
    properties
        Name
        E1
        E2
        G12
        v12       
        S
        Q
    end
    
    properties (Hidden)
        % General Stress Allowables
        s1t
        s1c
        s2t
        s2c
        t12
        
        % Tsai Wu
        F1(1,1)
        F11(1,1)
        F2(1,1)
        F22(1,1)
        F12(1,1)
        F66(1,1)
        
        Verbose(1,1) logical = false
    end
    
    
    methods
        function obj = CLT_Material(varargin)
            p               =   inputParser;
%             p.KeepUnmatched =   true;
            
            addOptional(p,'E1',1, @(x) x > 0);
            addOptional(p,'E2',1, @(x) x > 0);
            addOptional(p,'G12',0.15, @(x) x > 0);
            addOptional(p,'v12',@(x) (x > 0 && x < 0.5));
            addOptional(p,'Name','Mat', @(x) ischar(x));
            addOptional(p,'Verbose',false, @(x) islogical(x));
       
            parse(p,varargin{:});
            
            obj.E1          =   p.Results.E1;
            obj.E2          =   p.Results.E2;
            obj.G12         =   p.Results.G12;
            obj.v12         =   p.Results.v12;
            obj.Name        =   p.Results.Name;
            obj.Verbose     =   p.Results.Verbose;
            
            if nargin == 10
                obj.calculateS();
                obj.calculateQ();
            end
        end
        
        function calculateS(obj)
            obj.S   =   [1/obj.E1,-obj.v12/obj.E1,0;
                         -obj.v12/obj.E1,1/obj.E2,0;
                         0,0,1/obj.G12];
            
            
        end
        
        function calculateQ(obj)
            if isempty(obj.S)
                obj.calculateS();
            end
            
            sub     =   obj.S(1,1)*obj.S(2,2) - obj.S(1,2)^2;
            Q11     =   obj.S(2,2)/sub;
            Q22     =   obj.S(1,1)/sub;
            Q12     =   obj.S(1,2)/sub;
            Q66     =   1/obj.S(3,3);
            
            obj.Q   =   [Q11,Q12,0;
                         Q12,Q22,0;
                         0,0,Q66];
        end
        
        function addFailureTW(obj,varargin)
            p               =   inputParser;
%             p.KeepUnmatched =   true;
            
            addOptional(p,'s1t',1, @(x) x > 0);
            addOptional(p,'s1c',1, @(x) x > 0);
            addOptional(p,'s2t',1, @(x) x > 0);
            addOptional(p,'s2c',1, @(x) x > 0);
            addOptional(p,'t12',1, @(x) x > 0);
            addOptional(p,'F12',0,@(x) (x >= 0)||(strcmpi(x,'auto')));
       
            parse(p,varargin{2:end});
            
            obj.s1t         =   p.Results.s1t;
            obj.s1c         =   p.Results.s1c;
            obj.s2t         =   p.Results.s2t;
            obj.s2c         =   p.Results.s2c;
            obj.t12         =   p.Results.t12;
            obj.F12         =   p.Results.F12;
            
            obj.F1   	=   1/obj.s1t - 1/obj.s1c;
            obj.F11     =   1/(obj.s1t*obj.s1c);
            obj.F2      =   1/obj.s2t - 1/obj.s2c;
            obj.F22     =   1/(obj.s2t * obj.s2c);
            obj.F66     =   1/(obj.t12)^2;
            
            if strcmpi(obj.F12,'auto')
                obj.F12     =   -(1/2) * sqrt(1 / (obj.s1t*obj.s1c*obj.s2t*obj.s2c));
            end
            
            if obj.F11*obj.F22 - obj.F12^2 > 0; obj.F12 = 0; end
        end
        
        function TW = calculateTW(obj,S12)
            s1  =   S12(1);
            s2  =   S12(2);
            t12 =   S12(3);
            TW  =   obj.F1 * s1 + obj.F2 * s2 + obj.F11 * s1^2 + obj.F22 * s2^2 + ...
                    obj.F66 * (t12^2) + 2*obj.F12*s1*s2;
            if TW < 0; TW = 0; end
        end

        function [E1out] = calculateE1Parametric(obj,varargin)
            p               =   inputParser;
            
            addOptional(p,'Angles',(@(x) x >= 0 && x <= 90));

            parse(p,varargin{:});
            
            x           =   p.Results.Angles;
            P(length(x))=   CLT_Ply();
            E1out       =   zeros(size(x));
            
            counter     =   0;
            for angle = x
                counter     =   counter + 1;
                P(counter)  =   CLT_Ply(Angle = angle,...
                                        Material = obj,...
                                        Thickness = 1);
                a               =   inv(P(counter).Qbar);
                E1out(counter)  =   1/a(1,1);         
            end

        end
        
    end
    
end