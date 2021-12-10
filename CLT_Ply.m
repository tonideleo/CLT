classdef CLT_Ply < handle
    properties 
        S12
        Sxy
        e12
        exy
        kxy
        
        Qbar
        T
        Angle
        Thickness
        z_MidPlane
        Material
        
        FlagFail = false
        Qbar_reduced
    end
    
    
    
    methods 
        function obj = CLT_Ply(varargin)
            if nargin == 0
                return
            end
            p               =   inputParser;
%             p.KeepUnmatched =   true;
            
            addOptional(p,'Angle',@isdouble);
%             addOptional(p,'Material', @(x) isa(x,'CLT_Material'));
            addOptional(p,'Material',[]);
            addOptional(p,'Thickness',@(x) x > 0);
            parse(p,varargin{:});
            
            obj.Angle        	=   p.Results.Angle;
            obj.Material        =   p.Results.Material;
            obj.Thickness       =   p.Results.Thickness;
            
            obj.Qbar    =   obj.qbar(obj.Material.E1,...
                                     obj.Material.E2,...
                                     obj.Material.v12,...
                                     obj.Material.G12,...
                                     obj.Angle);
            obj.T       =   obj.calculateT(obj.Angle);
        end
        
        function calculateStress(obj)
                obj.Sxy     =   obj.Qbar * (obj.exy + obj.z_MidPlane.* obj.kxy);
                obj.S12     =   (inv(obj.T))\obj.Sxy;
                % Here there might be issues with obj.Material which does
                % not calculate S automatically! Need to double check!
                obj.Material.calculateS();
                obj.e12     =   obj.Material.S * obj.S12; %This might be wrong if using PFA
        end
       
        
        function failedPly(obj,Reduction)
            obj.FlagFail    =   true;
            obj.Qbar    =   obj.qbar(obj.Material.E1 * Reduction(1),...
                                     obj.Material.E2 * Reduction(2),...
                                     obj.Material.v12 ,...
                                     obj.Material.G12 * Reduction(3),...
                                     obj.Angle);
        end
       
    end
    
    methods (Static)
        
        function [qb] = qbar(E1,E2,v12,G12,theta)

        v21=v12*(E2/E1);
        q11=E1/(1-(v12*v21));
        q22=E2/(1-(v12*v21));
        q12=(v12*E2)/(1-(v12*v21));
        q66=G12;

        c=cosd(theta);
        s=sind(theta);

        qb11=(q11*c^4)+(2*(q12+2*q66)*s^2*c^2)+(q22*s^4);
        qb12=((q11+q22-4*q66)*s^2*c^2)+(q12*(s^4+c^4));
        qb22=(q11*s^4)+(2*(q12+2*q66)*s^2*c^2)+(q22*c^4);
        qb16=((q11-q12-2*q66)*s*c^3)+((q12-q22+2*q66)*s^3*c);
        qb26=((q11-q12-2*q66)*s^3*c)+((q12-q22+2*q66)*s*c^3);
        qb66=((q11+q22-2*q12-2*q66)*s^2*c^2)+(q66*(s^4+c^4));

        qb=[qb11 qb12 qb16
            qb12 qb22 qb26
            qb16 qb26 qb66];

        end
        
        function T = calculateT(theta)
            x   =   cosd(theta);
            y   =   sind(theta);
            
            T   =   [x^2, y^2, 2*x*y;
                     y^2, x^2, -2*x*y;
                     -x*y, x*y, x^2 - y^2];
            
        end
 
    end
    
    
end