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
%% Initialize variables
% rho = zeros(length(satelliteNumbers))
%% Import navigation file
% includes parameters
navfiles = importNavigationFiles('0lov033b.04n');
%% Now it becomes time specific
% Steps 1->11 loop
numberOfSat = 11;
for satelliteCount = 1:numberOfSat
    [ rho(satelliteCount,:),...
        Xs(satelliteCount,:),...
        Ys(satelliteCount,:),...
        Zs(satelliteCount,:),...
        dtsL1_with_dtr(satelliteCount,:)] = satLandP( satelliteCount,p1_numbers(satelliteCount),navfiles,XA0,YA0,ZA0 );
end
%% Repetition of 11-16
% for i = 1:2 % iterations for receiver coordinate estimates
    %% 13. Compute elements of vector L (19).
    Lmat = p1_numbers - rho + c*dtsL1_with_dtr;
    %% 14. Compute elements of matrix A (20); a_x_to_s , a_y_to_s , a_z_to_s by (12)
    Afake = [(Xs - XA0),(Ys - YA0),(Zs - ZA0),ones(length(Xs),1)];
    %%
    unknownX = (Amat'*Amat)\(Amat'*Lmat);
    %%
    newXaYaZa = [XA0;YA0;ZA0] + unknownX(1:3,i);
    newXYZcell = num2cell(newXaYaZa);
    [XA0,YA0,ZA0] = newXYZcell{:};
    v(:,i) = - Amat*unknownX(:,i) + Lmat;
% end
%% 15. Estimate unknown parameters by (18)
% unknownParameterX1 = (Amat'*Amat)\(Amat'*Lmat);
%% 16. Update receiver coordinates by (22)
%% 17. Repeat steps 11 -16 until the solution has converged.
% The solution has converged if the following condition is
% fulfilled: (vTv)i −(vTv)i−1 <ε, where ε is a small number
% and depends on the numerical accuracy, ε = 1e-5 should
% suffice to preserve mm numerical precision of the computed
% coordinates; i is iteration number. The vector v is computed
% after step 13 by Equation (17).
condition = abs(v(:,end)'*v(:,end)-v(:,end-1)'*v(:,end-1));