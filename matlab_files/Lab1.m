%% Computation of receiver's position
addpath('/Users/kevin/SkyDrive/KTH Work/Period 4 2014/GNSS/Labs/L1, L2 - Computation of satellite and receiver position')
addpath('/Users/kevin/SkyDrive/KTH Work/Period 4 2014/GNSS/Labs/L1, L2 - Computation of satellite and receiver position/matlab_files/')
cd('/Users/kevin/SkyDrive/KTH Work/Period 4 2014/GNSS/Labs/L1, L2 - Computation of satellite and receiver position')
clear all;
clc;
c = 299792458; % speed of light (m/s)
mu = 3.986005e14; % universal gravitational parameter (m/s)^3
omega_e_dot = 7.2921151467e-5; % earth rotation rate (rad/s)
F = -4.442807633e-10; % s/m^1/2
%% Import P1 numbers and satellite numbers
% observation file
p1_numbers = importObsP1numbers('0lov033b.04o', 1371, 1392); % Doesn't change
satelliteNumbers = importObsSatelliteNumbers('0lov033b.04o', 1370, 1370); % Doesn't change
[XA0,YA0,ZA0] = importApproxPosition('0lov033b.04o',8, 8); % Doesnt change
approxPos = [XA0,YA0,ZA0];
%% Import navigation file
% includes parameters
navfiles = importNavigationFiles('0lov033b.04n');
%% Match up satellite number with row in data
satNumMatch = navfiles(1:8:96,1); % Order of satellite numbers import
sortedSatelliteNumbers = sortrows([satelliteNumbers',p1_numbers],1);

%% Now it becomes satelite specific
for satelliteNumberOrder = 1:7
    [ LmatInside(satelliteNumberOrder,:), ...
        AmatInside(satelliteNumberOrder,:),...
        rho_A0_to_s(satelliteNumberOrder,:),...
        Xs(satelliteNumberOrder,:),Ys(satelliteNumberOrder,:),...
        Zs(satelliteNumberOrder,:),P1(satelliteNumberOrder,:),...
        dtsL1_with_dtr(satelliteNumberOrder,:)]...
        = satLandP( satelliteNumberOrder,sortedSatelliteNumbers(satelliteNumberOrder,2),navfiles,XA0,YA0,ZA0 );
end
numbcheck = [sortedSatelliteNumbers(1:7,1),rho_A0_to_s]
fprintf('%3f\n',numbcheck)
%%
for i = 1:20
    rho_A0_to_s = sqrt(...
        (Xs - XA0).^2 + ... % x^2
        (Ys - YA0).^2 + ... % y^2
        (Zs - ZA0).^2   ... % z^2
        );
    fprintf('%3f\n',rho_A0_to_s)
    Amat = [(Xs-XA0)./rho_A0_to_s,...
        (Ys-YA0)./rho_A0_to_s,...
        (Zs-ZA0)./rho_A0_to_s,...
        rho_A0_to_s./rho_A0_to_s];
    Lmat = P1 - rho_A0_to_s + c*dtsL1_with_dtr;
    changeX = (Amat'*Amat)\(Amat'*Lmat);
    newXYZ = [XA0,YA0,ZA0] + changeX(1:3)';
    newxyzcell = num2cell(newXYZ);
    [XA0,YA0,ZA0] = newxyzcell{:};
    [histXA0(i,:),histYA0(i,:),histZA0(i,:)] = newxyzcell{:};
    v(:,i) = -Amat*changeX + Lmat;
    %%
    
    %     if i>1
    %         condition = abs(v(:,end)'*v(:,end)-v(:,end-1)'*v(:,end-1));
    %         if condition < 1e-4
    %             fprintf('value of X is %7.3f\n',XA0)
    %             fprintf('value of Y is %7.3f\n',YA0)
    %             fprintf('value of Z is %7.3f\n',ZA0)
    %             fprintf('change of X is %3f\n',XA0-approxPos(1))
    %             fprintf('change of Y is %3f\n',YA0-approxPos(2))
    %             fprintf('change of Z is %3f\n',ZA0-approxPos(3))
    %             return
    %         end
    %     end
end

%% 17. Repeat steps 11 -16 until the solution has converged.
% The solution has converged if the following condition is
% fulfilled: (vTv)i −(vTv)i−1 <ε, where ε is a small number
% and depends on the numerical accuracy, ε = 1e-5 should
% suffice to preserve mm numerical precision of the computed
% coordinates; i is iteration number. The vector v is computed
% after step 13 by Equation (17).
condition = abs(v(:,end)'*v(:,end)-v(:,end-1)'*v(:,end-1));
fprintf('value of X is %7.3f\n',XA0)
fprintf('value of Y is %7.3f\n',YA0)
fprintf('value of Z is %7.3f\n',ZA0)
fprintf('change of X is %3f\n',XA0-approxPos(1))
fprintf('change of Y is %3f\n',YA0-approxPos(2))
fprintf('change of Z is %3f\n',ZA0-approxPos(3))
for i = 1:length(histXA0)
    fprintf('change of X is %3f, %3f, %3f \n',...
        [histXA0(i,:),histYA0(i,:),histZA0(i,:)])
end