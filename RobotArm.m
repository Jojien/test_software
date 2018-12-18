classdef RobotArm < handle
% function : Read physical parameters and calculate gravity and inertia of the robots.
% mass     : Mass of link;
% center   : Center of link;
% inertia  : Moments of inertia;
% gravityEqu  : Gravity equation of link;
% inertiaEqu  : Inertia equation of link;
% By Jiang Lei ,2018/12/17

    properties(SetAccess = private)
        mass     ;
        center   ;
        inertia  ;
        gravityEqu  ;
        inertiaEqu  ;
    end
    
    properties(Constant,Hidden)
        l1 = 0.34690;
        l2 = 0.200  ;
        l3 = 0.390  ;
        l4 = 0.12542;
        l5 = 0.4938 ;
        l6 = 0.009  ;
        g0 = [0, 0, -9.8];
        tilt        = pi/6   ;
        thetaZero   = 21.1635/180*pi;
        jnt2ZeroOff = 70/180*pi;
    end
    
    properties(Hidden)
        linkParms;
        filePath ;
        j3offset ;
    end
    
    methods
        function obj = RobotArm(arDir,arVersion)
            % arDir    : Folder where the parameter files are located.
            % arVersion: Robot number.
            obj.filePath = strcat(arDir,'\',arVersion);
            obj.j3offset = pi/2 - obj.tilt - obj.thetaZero;
            for i = 1:6
                CalcMassProperity (obj,strcat('link',num2str(i)),'all');
                obj.linkParms{i}.mass    = obj.mass  ;
                obj.linkParms{i}.center  = obj.center  ;
                obj.linkParms{i}.inertia = obj.inertia  ;
            end
        end
                
        function CalcMassProperity (obj,arLinkIdx,arObjName)
            % function : Read physical parameters of the robot.
            % arLinkIdx: The target joint you want to read.
            % arVersion: Robot number.
            massSign    = 'Mass =';
            centerSign  = 'center of mass: ( ';
            inertiaSign = 'Taken at the center of mass and aligned with the output coordinate system.';
            readMassRequired    = false;
            readcenterRequired  = false;
            readinertiaRequired = false;                                                                                            
            readFilePath = strcat(obj.filePath,'\',arLinkIdx,'.txt');
            if strcmpi(arObjName,'all')
                readMassRequired    = true;
                readcenterRequired  = true;
                readinertiaRequired = true;
            elseif strcmpi(arObjName,'mass')
                readMassRequired    = true;
            elseif strcmpi(arObjName,'center')
                readcenterRequired  = true;
            elseif strcmpi(arObjName,'inertia')
                readinertiaRequired = true;
            else
                warning('Please confirm the input parameters');
                return;
            end
            
            fid = fopen(readFilePath,'r');
            if fid<0
                error('Failed to open the file!');
            else
                tline = fgetl(fid);
                while ischar(tline)
                    if readMassRequired
                        if strncmpi(tline,massSign,6)                                       %Search for mass
                            obj.mass = sscanf(tline,'%*s%*s%f');
                            tline    = fgetl(fid);
                            if isempty(obj.mass)
                                error('Failed to read the data!');
                            end
                        end
                    end
                    
                    if readcenterRequired
                        if strncmpi(tline,centerSign,16)                                    %Search for center of mass
                            obj.center.x = sscanf(fgetl(fid),'%*s%*s%f');
                            obj.center.y = sscanf(fgetl(fid),'%*s%*s%f');
                            obj.center.z = sscanf(fgetl(fid),'%*s%*s%f');
                            tline = fgetl(fid);
                            if isempty(obj.center)
                                error('Failed to read the data!');
                            end
                        end
                    end
                    
                    if readinertiaRequired
                        if strncmpi(tline,inertiaSign,74)                                   %Search for center of moments of inertia
                            obj.inertia(1,1:3) = (sscanf(fgetl(fid),'%*s%*s%f%*s%*s%f%*s%*s%f'))';
                            obj.inertia(2,1:3) = (sscanf(fgetl(fid),'%*s%*s%f%*s%*s%f%*s%*s%f'))';
                            obj.inertia(3,1:3) = (sscanf(fgetl(fid),'%*s%*s%f%*s%*s%f%*s%*s%f'))';
                            tline = fgetl(fid);                 %#ok<*NASGU>
                            if isempty(obj.inertia)
                                error('Failed to read the data!');
                            end
                        end
                    end
                    tline = fgetl(fid);
                end
                if fclose(fid) ~= 0
                    error('Failed to close file!');
                end
            end
        end
        
        function gravityEqu = CalcGravityEqu (obj,arJoint)
            % function : Calcuate gravity of the robot.
            % arLinkIdx: The target joint you want to calculate.
            syms theta1 theta2 theta6 theta7 theta8 d5 a alph d theta
            DhTable = obj.GetDhTable();
            T = obj.CalcTransferMatrix();
            %Center of mass(4*1).
            for i = 1:6
                rUnit = 0.001;  
                r{i}  = [obj.linkParms{i}.center.x*rUnit, obj.linkParms{i}.center.y*rUnit, obj.linkParms{i}.center.z*rUnit, 1]'; %#ok<*AGROW>
            end
            %Convert to 0 coordinate system.
            T06 = T{2}*T{3}*T{4}*T{5}*T{6}*T{7};
            T05 = T{2}*T{3}*T{4}*T{5}*T{6};
            T04 = T{2}*T{3}*T{4}*T{5};
            T03 = T{2}*T{3}*T{4};
            T02 = T{2}*T{3};
            T{2}= T{2};      
            r06 = T06 *r{6};
            r05 = T05 *r{5};
            r04 = T04 *r{4};
            r03 = T03 *r{3};
            r02 = T02 *r{2};
            r01 = T{2}*r{1};
            
            potEn = subs(obj.linkParms{1}.mass*obj.g0(3)*r01(3) + obj.linkParms{2}.mass*obj.g0(3)*r02(3) +...    %Calculate the gravitational potential energy.
                         obj.linkParms{3}.mass*obj.g0(3)*r03(3) + obj.linkParms{4}.mass*obj.g0(3)*r04(3) + ...
                         obj.linkParms{5}.mass*obj.g0(3)*r05(3), {theta2}, {theta2 + obj.jnt2ZeroOff});
            number = arJoint(isstrprop(arJoint,'digit'));
            if DhTable(str2double(number)).type == 1
                argument =  strcat('theta',number);
                gravityEqu  = -obj.ChopSyms(vpa(simplify(expand(subs(diff(potEn,argument))))));    %Calculate partial derivative of joint.         
            end
            
            if DhTable(str2double(number)).type == 0
                argument =  strcat('d',number+2);
                gravityEqu  = -obj.ChopSyms(vpa(expand(simplify(subs(diff(potEn,argument))))));
            end
            
        end
       
        function inertiaEqu = CalcInertiaEqu (obj,arJoint)
            % function : Calcuate inertia of the robot.
            % arLinkIdx: The target joint you want to calculate.
            syms theta1 theta2 theta6 theta7 theta8 d5 a alph d theta
            DhTable = obj.GetDhTable();
            %Rotation matrix.
            T = obj.CalcTransferMatrix();
            %Center of mass(4*1).
            for i = 1:6
                rUnit = 0.001;  
                r{i}  = [obj.linkParms{i}.center.x*rUnit, obj.linkParms{i}.center.y*rUnit, obj.linkParms{i}.center.z*rUnit, 1]'; %#ok<*AGROW>
            end
                iUnit = 1e-6;
            %Calculate the gravitational potential energy.
            if strcmpi(arJoint,'joint1')
                energy =  obj.ConvertInertia(eye(4),              r{1},   obj.linkParms{1}.mass,  obj.linkParms{1}.inertia*iUnit) + ...
                          obj.ConvertInertia(T{3},                r{2},   obj.linkParms{2}.mass,  obj.linkParms{2}.inertia*iUnit) + ...
                          obj.ConvertInertia(T{3}*T{4},           r{3},   obj.linkParms{3}.mass,  obj.linkParms{3}.inertia*iUnit) + ...
                          obj.ConvertInertia(T{3}*T{4}*T{5},      r{4},   obj.linkParms{4}.mass,  obj.linkParms{4}.inertia*iUnit) + ...
                          obj.ConvertInertia(T{3}*T{4}*T{5}*T{6}, r{5},   obj.linkParms{5}.mass,  obj.linkParms{5}.inertia*iUnit) ;
                inertiaEqu = subs(energy(3,3),{theta2},{theta2 + 70/180*pi});
                inertiaEqu = obj.ChopSyms(simplify(vpa(inertiaEqu)));
            elseif strcmpi(arJoint,'joint2')
                auxFrame34 = [1,0,0,-obj.l4;0,1,0,(obj.l2^2-obj.l4^2)^0.5;0,0,1,0;0,0,0,1];
                energy = obj.ConvertInertia(eye(4),     r{2}, obj.linkParms{2}.mass, obj.linkParms{2}.inertia*iUnit)...
                    + obj.ConvertInertia(auxFrame34, r{4}, obj.linkParms{4}.mass, obj.linkParms{4}.inertia*iUnit);
                a = auxFrame34*T{6};
                inertiaEqu = energy(3,3) + obj.linkParms{3}.mass * T{4}(1,4)^2 + obj.linkParms{5}.mass *(a(1,4)^2+a(2,4)^2);
                inertiaEqu = subs(inertiaEqu,{theta2},{theta2 + 70/180*pi});
                inertiaEqu = obj.ChopSyms(simplify(expand(vpa(inertiaEqu))));
            else
                warning('Please confirm the input joint number');
                return;
            end
        end
        
        function T = CalcTransferMatrix(obj)
         % function : Calcuate rotation matrix..  
         DhTable = obj.GetDhTable;
         syms a alph d theta theta1 theta2 theta6 theta7 theta8 d5
         TT   = [cos(theta),            -sin(theta),              0,                  a        ;  %Transformation matrix of two coordinate systems.
                sin(theta)*cos(alph),    cos(theta)*cos(alph),   -sin(alph),     -sin(alph)*d;
                sin(theta)*sin(alph),    cos(theta)*sin(alph),    cos(alph),      cos(alph)*d;
                0,                        0,                        0,                  1       ];
         T{1} = subs(TT, {a, alph, d, theta}, {DhTable(1).a, DhTable(1).alpha, DhTable(1).d,      DhTable(1).theta         });  %TB0 0-to-base transformation
         T{2} = subs(TT, {a, alph, d, theta}, {DhTable(2).a, DhTable(2).alpha, DhTable(2).d,      DhTable(2).theta + theta1});  %T01 1-to-0    transformation
         T{3} = subs(TT, {a, alph, d, theta}, {DhTable(3).a, DhTable(3).alpha, DhTable(3).d,      DhTable(3).theta + theta2});  %T12 2-to-1    transformation
         T{4} = subs(TT, {a, alph, d, theta}, {DhTable(4).a, DhTable(4).alpha, DhTable(4).d,      DhTable(4).theta - theta2});  %T23 3-to-2    transformation
         T{5} = subs(TT, {a, alph, d, theta}, {DhTable(5).a, DhTable(5).alpha, DhTable(5).d,      DhTable(5).theta + theta2});  %T34 4-to-3    transformation
         T{6} = subs(TT, {a, alph, d, theta}, {DhTable(6).a, DhTable(6).alpha, DhTable(6).d + d5, DhTable(6).theta         });  %T45 5-to-4    transformation
         T{7} = subs(TT, {a, alph, d, theta}, {DhTable(7).a, DhTable(7).alpha, DhTable(7).d,      DhTable(7).theta + theta6});  %T56 6-to-5    transformation
         T{8} = subs(TT, {a, alph, d, theta}, {DhTable(8).a, DhTable(8).alpha, DhTable(8).d,      DhTable(8).theta + theta7});  %T67 7-to-6    transformation
         T{9} = subs(TT, {a, alph, d, theta}, {DhTable(9).a, DhTable(9).alpha, DhTable(9).d,      DhTable(9).theta + theta8});  %T78 8-to-7    transformation
               
        end
        
        function refVal = ConvertInertia(~,arT,arC,arM,arI)
            % function : Convert the inertia.
            interVal = (arT(1:3,1:4)*arC)'*arT(1:3,1:4)*arC;
            refVal   =  arM*([interVal,0,0;0,interVal,0;0,0,interVal]- arT(1:3,1:4)*arC*(arT(1:3,1:4)*arC)')...
                      + arT(1:3,1:3)*[arI(1,1),-arI(1,2),-arI(1,3);-arI(2,1),arI(2,2),-arI(2,3);-arI(3,1),-arI(3,2),arI(3,3)]*arT(1:3,1:3)';
        end
              
        function DhTable = GetDhTable(obj)
            % function : DH table.
            % type     : 1:Translational joints. 
            %            0:Rotational joints.
            % a,alpha,d,theta:Typical DH parameters.         
            DhTable = struct('type',  {1,      1,              0,     1,       1,             1,      1,        1,      1     },...   
                             'a',     {0,      0,              0,     obj.l2,  obj.l3,        obj.l4, 0,        0,      obj.l6},...
                             'alpha', {0,     -pi/2-obj.tilt, -pi/2,  0,       0,            -pi/2,   0,        pi/2,  -pi/2  },...
                             'd',     {-0.16,  obj.l1,         0,     0,       0,             0,      -obj.l5,  0,      0     },...
                             'theta', {0,     -pi/2,           0,    -pi/2,   -obj.j3offset,  0,      0,       -pi/2,   pi/2  });
        end
        
        function PlotGravity(obj,arJoint)
            % function : Plot the gravity of the robot.
            syms theta1 theta2 theta6 theta7 theta8 d5
            obj.gravityEqu = obj.CalcGravityEqu(arJoint);
            [fx,fy,fd] = meshgrid(linspace(deg2rad(-120),deg2rad(120),20),linspace(-40/180*pi-obj.jnt2ZeroOff,pi/2-obj.jnt2ZeroOff,20),linspace(0,0.4,20));
            fz         = double(subs(obj.gravityEqu, {theta1,theta2,d5}, {fx,fy,fd}));
            [fzMax,maxIndx] = max(fz(:));
            [fzMin,minIndx] = min(fz(:));
            fxmax = fx(maxIndx);
            fymax = fy(maxIndx);
            fxmin = fx(minIndx);
            fymin = fy(minIndx);
            %Gravity torque distribution:
            fz      = subs(obj.gravityEqu, {d5}, {0.4});
            [fx,fy] = meshgrid(linspace(deg2rad(-120),deg2rad(120),20),linspace(-40/180*pi-obj.jnt2ZeroOff,pi/2-obj.jnt2ZeroOff,20));
            fz      = double(subs(fz, {theta1,theta2}, {fx,fy}));
            figure,subplot(10,1,1:8),surf(fx*57.23 , fy*57.23 , fz);
            obj.SetLabel( '1st Joint Position /deg', '2nd Joint Position /deg', 'Gravity Torque tau/Nm', strcat(upper(arJoint(1)),arJoint(2:end),' Gravity Torque Distribution:'));
            if strcmpi(arJoint,'joint3')
                obj.SetLabel( '1st Joint Position /deg', '2nd Joint Position /deg', 'Gravity Torque tau/N', strcat(upper(arJoint(1)),arJoint(2:end),' Gravity Torque Distribution:'));
            end
            thandle = subplot(10,1,10);
            pos     = get(thandle,'position');
            delete(thandle);
            string = ['Max.gravity torque(Nm): ',num2str(roundn(fzMax,-2)),' at ','[',num2str(roundn(fxmax*57.23,-2)),',',num2str(roundn(fymax*57.23,-2) ),',',num2str(roundn(0.4,-2) ),'].        ',...
                'Min.gravity torque(Nm): ',num2str(roundn(fzMin,-2)),' at ','[',num2str(roundn(fxmin*57.23,-2) ),',',num2str(roundn(fymin*57.23,-2) ),',',num2str(roundn(0.4,-2) ),'].'];
            annotation('textbox',pos,'LineStyle','-','FitBoxToText','on','String',string,'FontSize',14);
        end
        
        function PlotInertia(obj,arJoint)
            % function : Plot inertia of the robot..
            syms theta1 theta2 theta6 theta7 theta8 d5
            obj.inertiaEqu = obj.CalcInertiaEqu(arJoint);
            [fy,fd] = meshgrid(linspace(-110/180*pi,pi/2-70/180*pi,20),linspace(0,0.4,20));
            fz      = double(subs(obj.inertiaEqu, {theta2,d5}, {fy,fd}));
            [fzMax,maxIndx] = max(fz(:));
            [fzMin,minIndx] = min(fz(:));
            fymax = fy(maxIndx);
            fdmax = fd(maxIndx);
            fymin = fy(minIndx);
            fdmin = fd(minIndx);
            %Inertia distribution:
            figure,subplot(10,1,1:8),surf(fy*57.23 , fd , fz);
                obj.SetLabel( '2nd Joint Position /deg', '3rd Joint Position /m', '1st Joint Load Inertia /kg.m^2', strcat(upper(arJoint(1)),arJoint(2:end),' Inertia Distribution:'));
           if strcmpi(arJoint,'joint2')
                obj.SetLabel( '2nd Joint Position /deg', '3rd Joint Position /m', '2st Joint Load Inertia /kg.m^2', strcat(upper(arJoint(1)),arJoint(2:end),' Inertia Distribution:')) ;   
           end
            set(get(gca, 'XLabel'), 'Rotation', 0)
            set(get(gca, 'YLabel'), 'Rotation', 0)
            thandle = subplot(10,1,10);
            pos     = get(thandle,'position');
            delete(thandle);
            string = ['Max.inertia(kg.m^2): ',num2str(roundn(fzMax,-2)),' at ','[',num2str(roundn(fymax*57.23,-2)),',',num2str(roundn(fdmax,-2) ),'].        ',...
                      'Min.inertia(kg.m^2): ',num2str(roundn(fzMin,-2)),' at ','[',num2str(roundn(fymin*57.23,-2)),',',num2str(roundn(fdmin,-2) ),'].'];
            annotation('textbox',pos,'LineStyle','-','FitBoxToText','on','String',string,'FontSize',14);
         end
        
        function SetLabel(~, arLableX, arLableY, arLableZ, arTitle)
             % function : Set label of figure.
             xlim = get(gca,'XLim');
             ylim = get(gca,'YLim');
             zlim = get(gca,'ZLim');
             xlabel (arLableX, 'FontSize', 14, 'position',[xlim(1)+0.5*(xlim(2)-xlim(1)),ylim(1)-0.3*(ylim(2)-ylim(1)),zlim(1)]);
             ylabel (arLableY, 'FontSize', 14, 'position',[xlim(1)-0.45*(xlim(2)-xlim(1)),ylim(1)+0.3*(ylim(2)-ylim(1)),zlim(1)]);
             zlabel (arLableZ, 'FontSize', 14);
             title  (arTitle,  'FontSize', 18);

%              zlim = get(gca,'ZLim');
%              set(gca, '',  0.5*xlim)%[xlim(1)+0.5*(xlim(2)-xlim(1)),ylim(1),zlim(1)]
%              set(gca, 'YTick',  0.5*ylim)%xlim(1),ylim(1)+0.5*(ylim(2)-ylim(1)),zlim(1)]
         end
         
        function chSyms = ChopSyms(~,arSyms)
            %function : Replace numbers near zero in expression with exact integer 0.
            syms theta1 theta2 theta6 theta7 theta8 d5
            [c,t]  = coeffs(arSyms);
             c     = (round(c*10000))/10000;
            chSyms = vpa(dot(c,t)); 
        end
    end
    
end