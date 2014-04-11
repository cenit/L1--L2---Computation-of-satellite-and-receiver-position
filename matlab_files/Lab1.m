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
count = 1;
for i = 1:length(satNumMatch)
    %% Steps 1-11 done inside
    if cell2mat(satNumMatch(i))==sortedSatelliteNumbers(count,1)
    [ Lmat(count,:), ...
        Amat(count,:),...
        rho(count,:),...
        Xs(count,:),Ys(count,:),...
        Zs(count,:),P1(count,:),...
        dtsL1_with_dtr(count,:),...
        change_tsv(count,:),...
        ts(count,:),...
        tAtoS(count,:)]...
        = satLandP( i,sortedSatelliteNumbers(count,2),navfiles,XA0,YA0,ZA0);
        count = count + 1;
    else
        fprintf('No data for Satellite%3d\n',cell2mat(satNumMatch(i)))
    end
end

%% %% 17. Repeat steps 11 -16 until the solution has converged.
% The solution has converged if the following condition is
% fulfilled: (vTv)i −(vTv)i−1 <ε, where ε is a small number
% and depends on the numerical accuracy, ε = 1e-5 should
% suffice to preserve mm numerical precision of the computed
% coordinates; i is iteration number. The vector v is computed
% after step 13 by Equation (17).
for i = 1:10
    changeX = (Amat'*Amat)\(Amat'*Lmat);
    v(:,i) = -Amat*changeX + Lmat;
    newXYZ = [XA0,YA0,ZA0] + changeX(1:3)';
    newxyzcell = num2cell(newXYZ);
    [XA0,YA0,ZA0] = newxyzcell{:};
    v(:,i) = -Amat*changeX + Lmat;
    rho = sqrt(...
        (Xs - XA0 + omega_e_dot * YA0 * tAtoS).^2 + ... % x^2
        (Ys - YA0 - omega_e_dot * XA0 * tAtoS).^2 + ... % y^2
        (Zs - ZA0).^2   ... % z^2
        );
    Amat = [-(Xs - XA0)./rho,...
        -(Ys-YA0)./rho,...
        -(Zs-ZA0)./rho,...
        rho./rho];
    Lmat = P1 - rho + c*dtsL1_with_dtr;
    %% 17. Convergence condition
    if i>1
        condition = abs(v(:,end)'*v(:,end)-v(:,end-1)'*v(:,end-1));
        if condition < 1e-4
            fprintf('Convergence condition met = %d\n',condition)
            fprintf('X = %7.3f\n',XA0)
            fprintf('Y = %7.3f\n',YA0)
            fprintf('Z = %7.3f\n',ZA0)
            return
        end
    end
end
