classdef CLT_Laminate < handle
    properties
        Stack
        Material
        Thickness
        Symmetric
        RepeatLeft
        RepeatRight
        NumberOfPlies
        
        Plies = CLT_Ply()
        A (3,3) double = zeros(3,3)
        B (3,3) double = zeros(3,3)
        D (3,3) double = zeros(3,3)
        ABD (6,6) double = zeros(6,6)
        abd (6,6) double = zeros(6,6)
        a (3,3) double = zeros(3,3)
        b (3,3) double = zeros(3,3)
        d (3,3) double = zeros(3,3)
    end
    
    properties(Hidden)
        LengthDisplayOutput(1,1) uint16 = 70
        Verbose(1,1) logical = false
        
    end
    
    methods
        function obj = CLT_Laminate(varargin)
            
            obj.LicenseMessage();
            
            if nargin == 0; return; end
            
            p               =   inputParser;
%             p.KeepUnmatched =   true;
            
            addOptional(p,'Stack',[]);
            addOptional(p,'Material', @(x) isa(x,'CLT_Material'));
            addOptional(p,'Thickness',0.15);
            addOptional(p,'Symmetric',false,@(x) islogical(x));
            addOptional(p,'RepeatLeft',1, @(x) x > 0);
            addOptional(p,'RepeatRight',1, @(x) x > 0);
            addOptional(p,'Verbose',false,@(x) islogical(x));
            
            parse(p,varargin{:});
            
            obj.Stack        	=   p.Results.Stack;
            obj.Material        =   p.Results.Material;
            obj.Thickness     	=   p.Results.Thickness;
            obj.Symmetric    	=   p.Results.Symmetric;
            obj.RepeatLeft    	=   p.Results.RepeatLeft;
            obj.RepeatRight    	=   p.Results.RepeatRight;
            obj.Verbose       	=   p.Results.Verbose;
            
            if obj.Symmetric
                obj.NumberOfPlies   =   (length(obj.Stack) * obj.RepeatLeft) * 2 * obj.RepeatRight;
            else
                obj.NumberOfPlies   =   length(obj.Stack) * obj.RepeatLeft * obj.RepeatRight;
            end
            
            if length(obj.Material) == 1
                obj.Material(1 : length(obj.Stack))     =   obj.Material;
            end
            
            if length(obj.Thickness) == 1
                obj.Thickness   =   ones(1,length(obj.Stack)).*obj.Thickness;
            end
            
            obj.parseStackingSequence();
            obj.parsePly();
            obj.calculateABD()
            
        end
        
        function parseStackingSequence(obj)
            stack           =   obj.Stack;
            repeat_left     =   obj.RepeatLeft;
            repeat_right    =   obj.RepeatRight;
            thickness       =   obj.Thickness;
            material        =   obj.Material;
            
            if obj.Symmetric
            explicit_stack      =   repmat([repmat(stack,1,repeat_left),fliplr(repmat(stack,1,repeat_left))],1,repeat_right);
            explicit_thickness  =   repmat([repmat(thickness,1,repeat_left),fliplr(repmat(thickness,1,repeat_left))],1,repeat_right);
            explicit_material   =   repmat([repmat(material,1,repeat_left),fliplr(repmat(material,1,repeat_left))],1,repeat_right);
            else
            % If not symmetric, only repeat right is considered!!
            explicit_stack      =   repmat(stack,1,repeat_right);
            explicit_thickness  =   repmat(thickness,1,repeat_right);
            explicit_material   =   repmat(material,1,repeat_right);
            end

            obj.Stack       =   explicit_stack;
            obj.Thickness   =   explicit_thickness;
            obj.Material    =   explicit_material;

        end
        
        function parsePly(obj)
                obj.Plies(1:length(obj.Stack))   =   CLT_Ply();
            for i = 1 : length(obj.Stack)
                obj.Plies(i)    =   CLT_Ply(Angle       =   obj.Stack(i),...
                                            Material    =   obj.Material(i),...
                                            Thickness   =   obj.Thickness(i));
                
            end
        end
        
        function calculateABD(obj)
            % Laminate Local Rerefence System
            h       =   sum(obj.Thickness);
            z       =   zeros(length(obj.Stack),1);
            z(1)    =   -h/2;
            for i= 2 : length(obj.Stack) + 1
                z(i)    =   z(i-1) + obj.Thickness(i-1); 
            end
            
            for i = 1 : obj.NumberOfPlies
                obj.Plies(i).z_MidPlane = z(i) + obj.Thickness(i)./2;
            end
            
            for i = 1 : length(obj.Stack)
                Qbar    =   obj.Plies(i).Qbar;
                obj.A 	=   obj.A + Qbar*(z(i+1)-z(i));
                obj.B	=   obj.B +(1/2)*Qbar*(z(i+1)^2-z(i)^2);
                obj.D	=   obj.D +(1/3)*Qbar*(z(i+1)^3-z(i)^3);
            end
            
            obj.a   =   inv(obj.A);
            if all(abs(obj.B(:)) <= 1e-6)
                obj.b   =   nan(3,3);
            else
                obj.b   =   inv(obj.B);
            end
            obj.d   =   inv(obj.D);
            
            obj.ABD     =   [obj.A, obj.B; obj.B, obj.D];
            obj.abd     =   [obj.a, obj.b; obj.b, obj.d];
        end
        
        function ProgressiveFailureAnalysis(obj,varargin)
            
            if obj.Verbose
                obj.printLineText('Progressive Failure Analysis');
            end
            
            p               =   inputParser;
%             p.KeepUnmatched =   true;
            
            addOptional(p,'Nx',1, @isnumeric);
            addOptional(p,'Ny',0, @isnumeric);
            addOptional(p,'Nxy',0, @isnumeric);
            addOptional(p,'Mx',0, @isnumeric);
            addOptional(p,'My',0, @isnumeric);
            addOptional(p,'Mxy',0, @isnumeric);
            addOptional(p,'Verbose',false, @islogical);
            
            parse(p,varargin{2:end});
            
            Nx       	=   p.Results.Nx;
            Ny        	=   p.Results.Ny;
            Nxy      	=   p.Results.Nxy;
            Mx       	=   p.Results.Mx;
            My          =   p.Results.My;
            Mxy        	=   p.Results.Mxy;
            verbose     =   p.Results.Verbose;
            
            MCD         =   max(abs([Nx,Ny,Nxy,Mx,My,Mxy]));
            Nxold       =   Nx;
            F           =   [Nx,Ny,Nxy,Mx,My,Mxy]./MCD;
            if all(F == 0 | isnan(F))
                fprintf('Initial Loads are unusable. Aborting...\n');
                if obj.Verbose; obj.printLine(); end
                return
            end
            
            tmax        =   MCD;
            nsteps      =   10000;
            dt          =   linspace(0,tmax * 5,nsteps);
            t           =   0;
            counter     =   1;
            flag_FPF    =   false;

            while ~all([obj.Plies.FlagFail]) && t < tmax * 5
                counter     =   counter + 1;
                if counter > nsteps
                    fprintf(2,'The applied load is not enough to finish the PFA.\nIncrease initial value!\n');
                    return
                end
                NM          =   dt(counter).*F;
                eps         =   obj.ABD\NM';
                
                for i = 1 : obj.NumberOfPlies
                    obj.Plies(i).exy    =   eps(1:3);
                    obj.Plies(i).kxy    =   eps(4:6);

                    obj.Plies(i).calculateStress();
                    TW  =   obj.Material(i).calculateTW(obj.Plies(i).S12);
                    if TW > 1 && ~obj.Plies(i).FlagFail
                        obj.Plies(i).failedPly([0.4 0.4 0.15]);
                        if verbose
                            fprintf('Ply %i (Angle = %i deg) FAILED at %.2f %% of given loads!\n',...
                                i,obj.Plies(i).Angle,(NM(1)/Nxold)*100);
                        end
                        if ~flag_FPF % NEEDS TO BE IMPLEMENTED FOR MORE DETAILED FPF 
%                             fprintf('Following First Fly Failure, the laminate would fail at %f\n',i);
                            flag_FPF    =   true;
                        else
                            
                        end
                    end

                end 
            end

            if verbose
                fprintf('ALL plies failed!\n');
            end
            
            if obj.Verbose
                fprintf('Progressive Failure Analysis Successfully Terminated!\n');
                obj.printLine();
            end
            
        end
        
        function [varargout] = calculateMacroBehavior(obj)
            Ex  =   1 / (sum(obj.Thickness) * obj.abd(1,1));
            Ey  =   1 / (sum(obj.Thickness) * obj.abd(2,2));
            Gxy =   1 / (sum(obj.Thickness) * obj.abd(6,6));
            vxy =   -obj.abd(2,1) / obj.abd(1,1);
            
            varargout{1}    =   Ex;
            varargout{2}    =   Ey;
            varargout{3}    =   Gxy;
            varargout{4}    =   vxy;
            
        end
        
        function [varargout] = calculateBuckling(obj,varargin)
            
            obj.printLineText('AE 553 - Fall 2021 - Guest Lecture Code')
            fprintf('WARNING: Buckling Load for Symmetric Laminate Only!\n');
            fprintf('This is a show demonstration for AE 553!\n');
            fprintf('No parts of this method were verified and validated!\n');
            fprintf('Use at your own risk!\n');
            
            
            p               =   inputParser;
            p.KeepUnmatched =   true;
            
            addOptional(p,'Lx',1, @isnumeric);
            addOptional(p,'Ly',0, @isnumeric);
            addOptional(p,'k',0, @isnumeric);
            addOptional(p,'BC',{'S4'}, @ischar);
            addOptional(p,'DisplayPlot',false, @islogical);
            addOptional(p,'DisplayTable',false, @islogical);
       
            parse(p,varargin{2:end});
            
            Lx       	=   p.Results.Lx;
            Ly        	=   p.Results.Ly;
            k         	=   p.Results.k;
            BC          =   p.Results.BC;
            flag_plot 	=   p.Results.DisplayPlot;
            flag_table  =   p.Results.DisplayTable;
            
            if ~strcmpi(BC,{'s4'})
                fprintf(2,'Boundary Condition NOT IMPLEMENTED!\n');
                return
            end
            
            R       =   Lx/Ly;
            
            D11     =   obj.D(1,1);
            D22     =   obj.D(2,2);
            D12     =   obj.D(1,2);
            D66     =   obj.D(3,3);
            
            Nc  =   @(n,m) pi^2 * (D11*m^4 + 2*(D12 + 2*D66)*m^2*n^2*R^2 + D22*n^4*R^4) / ...
                    (Lx^2 * (m^2 + k * n^2 * R^2));
            
            counter =   0;
            Ncrit   =   zeros(1,9);
            N       =   zeros(1,9);
            M       =   zeros(1,9);
            for i = 1 : 3
                for j = 1 : 3
                    counter = counter + 1;
                    Ncrit(counter)  =   Nc(i,j);
                    N(counter)      =   i;
                    M(counter)      =   j;
                end
            end
            
            idx     =   find(Ncrit <= 0);
            
            if ~isempty(idx)
                fprintf('WARNING: Some critical buckling values were negative.\n');
                fprintf('         The plot (if requested) will not be shown!\n');
                flag_plot   =   false;
            end
            % Comment if you want to keep negative Eigenvalues
            Ncrit(idx)  =   [];
            N(idx)      =   [];
            M(idx)      =   [];
            
            Res     =   table(Ncrit',N',M');
            Res.Properties.VariableNames{1} = 'Critical Load';
            Res.Properties.VariableNames{2} = 'n';
            Res.Properties.VariableNames{3} = 'm';
            
            if flag_table
                ResultsTable = sortrows(Res,1)
            end
            
            if flag_plot
                figure
                warning('off')
                fsurf(Nc,[1 3 1 3])
                box on;
                grid on; grid minor;
                options = {'Interpreter','latex'};
                xlabel('n',options{:}); ylabel('m',options{:}); zlabel('Critical Load',options{:});
                set(gca,'TickLabelInterpreter','latex')
                set(gca,'FontSize',20)

                warning('on')
            end
            
        end
        
        function printLine(obj,varargin)
            if nargin > 1
                str     =   varargin{1};
            else
                str     =   '=';
            end
            
            fprintf(repelem(str,obj.LengthDisplayOutput));
            fprintf('\n');
        end        
        
        function printLineText(obj,text,varargin)
            if nargin > 2
                str     =   varargin{1};
            else
                str     =   '=';
            end
            
            len     =   length(text);
            toprint =   repelem(str,obj.LengthDisplayOutput/2 - 1 - ...
                        floor(len/2)); 
                    
            fprintf(toprint);
            fprintf([' ',text,' ']);
            
            if rem(len,2) == 1 %odd
                fprintf(toprint(1:end-1));
            else %even
                fprintf(toprint); 
            end
            
            fprintf('\n');
        end
        
        function LicenseMessage(obj)
            obj.printLineText('Classical Lamination Theory');
            fprintf('Copyright (C) 2021 Antonio Alessandro Deleo (adeleo@uw.edu)\n');
            fprintf('Department of Aeronautics & Astronautics\n');
            fprintf('University of Washington, Seattle, WA\n');
            fprintf('This program is free software: you can redistribute it and/or modify\n');
            fprintf('it under the terms of the GNU General Public License as published by\n');
            fprintf('the Free Software Foundation, either version 3 of the License, or\n');
            fprintf('any later version.\n');
            fprintf('This program is distributed in the hope that it will be useful,\n');
            fprintf('but WITHOUT ANY WARRANTY; without even the implied warranty of\n');
            fprintf('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n');
            fprintf('GNU General Public License for more details.\n');
            fprintf('You should have received a copy of the GNU General Public License\n');
            fprintf('along with this program.  If not, see <http://www.gnu.org/licenses/>.\n');
            obj.printLine();
            fprintf('\n');
            
        end

    end

 
end