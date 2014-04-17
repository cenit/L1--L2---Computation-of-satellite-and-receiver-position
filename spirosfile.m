    % SATELITE 24
    clear all
c=299792458;%light speed
m=3.986005e14;%Earth universal gravitational parameter
PAs=25384646.34743;                                  %Pseudorange
tA=((60+10)*60)+86400;                                  %Personal epoch in seconds +day seconds


%______Question 1,2__________________

DtAs= (PAs)/c;                                                   %signal propagation time (Eq.13)

ts=(tA)-((PAs)/c);                                               %nomimnal transmission time (Eq.14)


af0=1.093372702600e-6;                                % sv clock bias
af1=2.955857780760e-12;                              % sv clock drift
af2=0;                                                                  % sv clock drift rate
toe=9.360000000000e4;                                  % table 3
toc= 86400+(3600*2);                                          %TIME OF CLOCK



TGD=-9.313225746160e-10;                            %Timing group delay

%___________ Question 3_________________________

Dtsv= af0+af1*(ts-toc)+af2*(ts-toc)*(ts-toc);    % (3) Equation 25
dtL1s=Dtsv-TGD;                                                   % (3)satelite clock correction (Eq.24)

 %______Question (4)________________________                                       

  TRts= ts-dtL1s;                                                 %System transmission time ts(Eq.15)
  
  
  
  
%_______ECCENTRIC ANOMALLY (5)_______________

Dn=4.158030341510e-9;
A=(5.153573421480e3)^2;
Aabs=sqrt(A*A);                                                  %semimajor axis
n0=sqrt(m/(Aabs^3));                                         %computed mean motion
n=n0+Dn;                                                               %corrected mean motion
tk=ts-toe;                                                               %time for ephemeris reference epoch
M0=-1.234439725290e0;
ec=9.439246961850e-3;
Mk=M0+n*tk;                                                      %MEAN ANOMALLY
Ek(1) = Mk;
for i = 1:5 
    Ek(i+1) = Mk + ec*sin(Ek(i))
end
Ek = Ek(end);                                                     %(5) Eccentric anommaly
  
  
%__Computing Dtr and new system transmission time ts  (6)_____

F=-4.442807633e-10;                                         % F value  eq.(27)
epsilon=9.439246961850e-3;                            %e value table 3

Dtr=F*epsilon*(sqrt(Aabs))*(sin(Ek));               % Eq.26

Dtsvnew=Dtsv+Dtr;                                               % New values eq.(24,25)
dtL1snew=Dtsvnew-TGD;                                         

TRtsnew=ts-dtL1snew;                                        %New system transmission time

%_______Question 7: Satelite  coordinates_______
sinvk= (sqrt(1-epsilon)* sin(Ek))/(1-(epsilon*cos(Ek)));    %table 2 True anomally
cosvk=(cos(Ek)-epsilon)/(1-(epsilon*cos(Ek)));                  %table 2
vk=atan2(sinvk,cosvk);                                                         %table 2
crs=1.055312500000e2;                                                 %table 3
crc=2.370625000000e2;                                                  %table 3

w=-1.409026029280D+00;                                  % omega parameter     (Table 2)
Fk=w+vk;                                                               % argument of latitude

drk=crs*sin(2*Fk)+crc*cos(2*Fk);                   %radius correction
rk=Aabs*(1-epsilon*cos(Ek))*drk;                    % corrected radius

cus=7.957220077510e-6;                                %table 3
cuc=5.569308996200e-6;                                 % table 3
duk=cus*sin(2*Fk)+cuc*cos(2*Fk);                 %argument of latitude correction
uk=Fk+duk;
xki=rk*cos(uk);                                                     % Position in orbital plane 
yki=rk*sin(uk);

OMEGAe=7.2921151467e-5;
OMEGAdot=-7.928544541420e-09;
OMEGA0=2.953800868980e+00;   
OMEGAk= OMEGA0+(OMEGAdot-OMEGAe)*tk-OMEGAe*toe; %corrected long. ascend. node table 2


cis=1.303851604460e-07;
cic=-1.676380634310e-07;
i0=9.719990595480e-01 ;
dik=cis*sin(2*Fk)+cic*cos(2*Fk);                      %inclination correction
IDOT=-1.821504444400e-11;                            %table 3
ik=i0+dik+(IDOT)*tk   ;                                      % corrected inclination

%INPUT THE REST OF THE PARAMETERS

Xk=xki*cos(OMEGAk)-(yki*cos(ik)*sin(OMEGAk))
Yk=xki*sin(OMEGAk)+(yki*cos(ik)*cos(OMEGAk))
Zk=yki*sin(ik)
