%% Computation of receiver's position
addpath('/Users/kevin/SkyDrive/KTH Work/Period 4 2014/GNSS/Labs/L1, L2 - Computation of satellite and receiver position')
clear all;
clc;
c = 299792458; % speed of light (m/s)
mu = 3.986005e14; % universal gravitational parameter (m/s)^3
omega_e_dot = 7.2921151467e-5; % earth rotation rate (rad/s)
F = -4.442807633e-10; % s/m^1/2
%% Import navigation file
navfiles = importNavigationFiles('0lov033b.04n');
i = 1:8:112;
sat = navfiles(1:8,:);
sat = cell2mat(sat);
%% Compute signal propagation time by (13)
P1 = 23640467.92143;
ta_nom = seconds_in_week(1,1,14,0); % 1 hour and 14 minutes
tAtoS = P1/c; % signal propagation time
%% Compute signal transmission time by (14)
ts_nom = ta_nom - tAtoS;
%% Compute satellite clock correction dtsL1
% by (24) and (25), neglect dtr
af0 = sat(1,2);
af1 = sat(1,3);
af2 = sat(1,4);
tgd = sat(7,3);
t_oc = seconds_in_week(1,2,0,0); %
tsv = af0 + af1*(ts_nom-t_oc)+af2*(ts_nom-t_oc)^2; % (25)
dtsL1 = tsv - tgd; % (24)
%% Compute ts using the correction from the step 3.
ts = ts_nom - dtsL1;
%% Compute eccentric anomaly (Table 2)
% ek = mk + ec*sin(ek)
dn = sat(2,3);
m0 = sat(2,4);
A = (sat(3,4))^2;
n0 = sqrt(mu/A^3); % Computed mean motion
n = n0 + dn;
ec = sat(3,2);
toe = sat(4,1);
tk = ts - toe;
tk = fixTk(tk); % if,then for table 2 of tk
mk = m0 + n*tk;
Ek = keplersEquation(mk,ec);
%% Compute dtr by (26) and ts by (15).
dtr = F*ec*sqrt(A)*sin(Ek); %(26)
ts_with_dtr = ts - dtr;
%% Compute satellite coordinates Xs, Ys, Zs, for time ts - Table 2
% Calculate rk
crs = sat(2,2);
crc = sat(5,2);
vk = atan2((sqrt(1-ec^2)*sin(Ek)/(1-ec*cos(Ek))),((cos(Ek)-ec)/(1-ec*cos(Ek))));
w = sat(5,3);
Phik = vk + w;
drk = crs*sin(2*Phik) + crc*cos(2*Phik);
rk = A*(1-ec*cos(Ek)) + drk;  % Corrected radius
% Calculate uk
cus = sat(3,3);
cuc = sat(3,1);
duk = cus*sin(2*Phik) + cuc*cos(2*Phik);
uk = Phik + duk;
% Calculate ik
i0 = sat(5,1);
cis = sat(4,4);
cic = sat(4,2);
dik = cis*sin(2*Phik) + cic*cos(2*Phik);
idot = sat(6,1);
ik = i0 + dik + idot*tk;
% Calculate omega's
omega0 = sat(4,3);
omegadot = sat(5,4);
omegak = omega0 + (omegadot-omega_e_dot)*tk - omega_e_dot*toe;
% Calculate xkp and ykp
xkp = rk*cos(uk);
ykp = rk*sin(uk);
% Calculate xk,yk,zk -> Xs, Ys, Zs for time ts
xk = xkp*cos(omegak) - ykp*cos(ik)*sin(omegak);
yk = xkp*sin(omegak) + ykp*cos(ik)*cos(omegak);
zk = ykp*sin(ik);
Xs = xk;
Ys = yk;
Zs = zk;
%% Compute satellite clock correction dtsL1 by (24) - (27)
dtsL1_with_dtr = dtsL1 + dtr; % (24)
%% 9. Compute tropospheric correction T_A_to_s (tA)
%% 10. Compute ionospheric correction I_A_to_s (tA)
%% 11. Compute approximate distance rho_A0_to_s (tA) by (11).
dts = 0; % terms with dts are negligible, so I set it to zero
XA0 = 3104219.4530; % 
YA0 = 998383.9820; % approximate receiver coordinates
ZA0 = 5463290.5080; %
rho_A0_to_s = sqrt(...
    (Xs - XA0 + omega_e_dot*YA0*dts)^2 + ... % x^2
    (Ys - YA0 + omega_e_dot*XA0*dts)^2 + ... % y^2
    (Zs - ZA0)^2   ... % z^2
);
% dtA = 0;
% rho_A_to_s = P1 + c*dtsL1_with_dtr - c*dtA; % (8) dtA =\= 0 
%% 12. Repeat steps 1 - 11 for all measured satellites.
%% 13. Compute elements of vector L (19).
%% 14. Compute elements of matrix A (20); a_x_to_s , a_y_to_s , a_z_to_s by (12)
%% 15. Estimate unknown parameters by (18)
%% 16. Update receiver coordinates by (22)
%% 17. Repeat steps 11 -16 until the solution has converged. 
% The solution has converged if the following condition is
% fulfilled: (vTv)i −(vTv)i−1 <ε, where ε is a small number
% and depends on the numerical accuracy, ε = 1e-5 should 
% suffice to preserve mm numerical precision of the computed
% coordinates; i is iteration number. The vector v is computed
% after step 13 by Equation (17).