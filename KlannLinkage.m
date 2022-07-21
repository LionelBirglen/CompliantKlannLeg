% Matlab Script for the kinetostatic analysis of a new version of Klann's
% robotic leg using compliant joints
%
% by Jerome Bastien and Lionel Birglen
% Polytechnique Montreal
% March 2022
%
% https://github.com/LionelBirglen/CompliantKlannLeg
%
% GNU GPL v2 license

close all;
clear all;


% Optimal geometric parameters, see paper for details
x=[39.6566   80.5157  -10.6073   18   48.4410   37.4648];
L2=x(1);
thetaAB=x(2)*pi/180;
alpha=x(3)*pi/180;
k0=x(4);
L6=x(5);
MEM=x(6)*pi/180;

% Other constant geometric parameters
L1=28;
L3=91;
L4=94;
L7=95;
L8=76;rho=2.1;
epaisseur=8;

% Point position lists
% All lists represent discretized data using a one degree step 
% of the actuator input angle, again see paper for nomenclature
listeM_prime=zeros(360,2); 
listeC=zeros(360,2); 
listeD=zeros(360,2); 
listeM=zeros(360,2); 
listeE=zeros(360,2);
listeB1=zeros(360,2); 
listeE1=zeros(360,2); 
listeD1=zeros(360,2); 
listeC1=zeros(360,2);
listeB=zeros(360,2);
listeM1=zeros(360,2);

% Joint angle list
liste_thetaA=zeros(360,1);
liste_thetaB=zeros(360,1); 
liste_thetaD=zeros(360,1); 
liste_thetaM=zeros(360,1); 
liste_thetaE=zeros(360,1);

% Kinematic coefficient list
liste_thetaA_theta0=zeros(360,1); 
liste_thetaB_theta0=zeros(360,1); 
liste_thetaD_theta0=zeros(360,1); 
liste_thetaM_theta0=zeros(360,1); 
liste_thetaE_theta0=zeros(360,1);

% Force, torque, and velocity list
liste_force=zeros(360,1); 
liste_couple_a_vide=zeros(360,1);
liste_vitesse=zeros(360,1);
liste_energie_potentielle=zeros(360,1);

% Fixed points of the linkage
O=[0, 0]; 
A=[-L4, 0]; 
B=[-L4 + L8 * cos(thetaAB), L8 * sin(thetaAB)];

% Other constants
V_envers=[0 0 1; 1 0 0; 0 1 0];
global mat_E
mat_E=[0 -1; 1 0];

% Simulation of a complete turn of the input crank, one degree dicrete step
for theta0=k0:1:359+k0
    % Point C computation and list update
    C=[L1*cos(theta0*pi/180),L1*sin(theta0*pi/180)];
    listeC(mod(theta0-1,360)+1,:)=C;
    
    % Point D
    d=sqrt((-L4-L1*cos(theta0*pi/180))^2+(-L1*sin(theta0*pi/180))^2);
    a=(L3^2-L2^2+d^2)/(2*d);
    h=sqrt(L3^2-a^2);
    x2=L1*cos(theta0*pi/180)+a*(-L4-L1*cos(theta0*pi/180))/d;
    y2=L1*sin(theta0*pi/180)+a*(-L1*sin(theta0*pi/180))/d;
    D=[x2+h*(-L1*sin(theta0*pi/180))/d,y2-h*(-L4-L1*cos(theta0*pi/180))/d];
    listeD(mod(theta0-1,360)+1,:)=D;
    
    % Point M
    XCD=D(1)-C(1);
    YCD=D(2)-C(2);
    M=[rho*(XCD*cos(alpha)-YCD*sin(alpha))+C(1),rho*(XCD*sin(alpha)+YCD*cos(alpha))+C(2)];
    listeM(mod(theta0-1,360)+1,:)=M;
     
    % Point E
    d=sqrt((-L4+L8*cos(thetaAB)-M(1))^2+(L8*sin(thetaAB)-M(2))^2);
    a=(L7^2-L6^2+d^2)/(2*d);
    h=sqrt(L7^2-a^2);
    x2=M(1)+a*(-L4+L8*cos(thetaAB)-M(1))/d;
    y2=M(2)+a*(L8*sin(thetaAB)-M(2))/d;
    E=[x2-h*(L8*sin(thetaAB)-M(2))/d,y2+h*(-L4+L8*cos(thetaAB)-M(1))/d];
    listeE(mod(theta0-1,360)+1,:)=E;
    
    % Point M'
    XEM=M(1)-E(1);
    YEM=M(2)-E(2);
    M_prime=[rho*(XEM*cos(MEM)-YEM*sin(MEM))+E(1),rho*(XEM*sin(MEM)+YEM*cos(MEM))+E(2)];
    listeM_prime(mod(theta0-1,360)+1,:)=M_prime;

    % Update of the joint angles
    liste_thetaA(mod(theta0-1,360)+1)=180 / pi * atan2(D(2), D(1)+L4);
    liste_thetaB(mod(theta0-1,360)+1)=180 / pi * (atan2(E(2) - L8 * sin(thetaAB), E(1) + L4 - L8 * cos(thetaAB)));
    xDA=A(1)-D(1); yDA=A(2)-D(2); xDC=C(1)-D(1); yDC=C(2)-D(2);
    liste_thetaD(mod(theta0-1,360)+1)=180 / pi * acos((xDA * xDC + yDA * yDC)/(sqrt(xDA^2 + yDA^2) * sqrt(xDC^2 + yDC^2)));
    xMD=D(1)-D(1); yMD=D(2)-E(2); xME=E(1)-M(1); yME=E(2)-M(2);
    liste_thetaM(mod(theta0-1,360)+1)=180 / pi * acos((xMD * xMD + yMD * yME)/(sqrt(xMD^2 + yMD^2) * sqrt(xME^2 + yME^2)));
    xEM=M(1)-E(1); yEM=M(2)-E(2); xEB=B(1)-E(1); yEB=B(2)-E(2);
    liste_thetaE(mod(theta0-1,360)+1)=180 / pi * acos((xEM * xEB + yEM * yEB)/(sqrt(xEM^2 + yEM^2) * sqrt(xEB^2 + yEB^2)));
    
    % Analytical computation of the velocities using planar screws    
    phi=atan2(D(2) - C(2), D(1) - C(1));

    mat_A=[prod_scal(vecteur(M(1), M(2), C(1), C(2), 1), vecteur(C(1), C(2), O(1), O(2), 0)), ...
        vecteur(M(1), M(2), C(1), C(2), 1); ...
        prod_scal(vecteur(M(1), M(2), D(1), D(2), 1), vecteur(D(1), D(2), O(1), O(2), 0)), ...
        vecteur(M(1), M(2), D(1), D(2), 1); ...
        prod_scal(vecteur(E(1), E(2), B(1), B(2), 1), vecteur(B(1), B(2), O(1), O(2), 0)), ...
        vecteur(E(1), E(2), B(1), B(2), 1)];

    mat_C=[L1*sin(phi-theta0*pi/180), cos(phi), sin(phi); 0, 1, 0; -L4, 0, 1];
    
    mat_d=[L1*sin(phi-theta0*pi/180); 0; 0];
    
    tau0D=[vecteur(M(1), M(2), D(1), D(2), 1)'; ...
        prod_scal(vecteur(M(1), M(2), D(1), D(2), 1), vecteur(D(1), D(2), O(1), O(2), 0))];
    
    mat_b=[prod_scal(vecteur(M(1), M(2), C(1), C(2), 1), vecteur(C(1), C(2), O(1), O(2), 0)); ...
        (mat_C\mat_d)'*V_envers*tau0D; 0];
    
    OM_prime=vecteur(M_prime(1), M_prime(2), O(1), O(2), 0);
    EOM_prime=mat_E*OM_prime';
    
    OE=vecteur(E(1), E(2), O(1), O(2), 0);
    OE_prime=mat_E*OE';
    
    mat_F=[1, 0, 0; -EOM_prime(1), 1, 0; -EOM_prime(2), 0, 1];
    mat_F_pour_E=[1, 0, 0; -OE_prime(1), 1, 0; -OE_prime(2), 0, 1];
    
    AF = mat_A*mat_F;
    AF_pour_E = mat_A*mat_F_pour_E;
    
    visseur=AF\mat_b;
    visseur_pour_E=AF_pour_E\mat_b;
    
    % Update of the foot velocity
    liste_vitesse(mod(theta0-1,360)+1)=sqrt(visseur(2)^2 + visseur(3)^2);
    
    % Computation of the kinematic coefficients for application of the 
    % virtual work principle with the compliant version of the leg
    liste_thetaB_theta0(mod(theta0-1,360)+1)=visseur_pour_E(3)/(E(1)-B(1));
    liste_thetaA_theta0(mod(theta0-1,360)+1)=L1*(sin(theta0*pi/180)-tan(phi)*cos(theta0*pi/180))/ ...
        (L2*(sin(liste_thetaA(mod(theta0-1,360)+1)*pi/180)-tan(phi)*cos(liste_thetaA(mod(theta0-1,360)+1)*pi/180)));
    thetaphi_theta0=L1*(sin(theta0*pi/180)-tan(liste_thetaA(mod(theta0-1,360)+1)*pi/180) * cos(theta0*pi/180)) / ...
        (L3 * (-sin(phi) + tan(liste_thetaA(mod(theta0-1,360)+1)*pi/180) * cos(phi)));
    liste_thetaD_theta0(mod(theta0-1,360)+1)=-liste_thetaA_theta0(mod(theta0-1,360)+1)+thetaphi_theta0;
    liste_thetaM_theta0(mod(theta0-1,360)+1)=visseur(1)+thetaphi_theta0*pi/180;
    liste_thetaE_theta0(mod(theta0-1,360)+1)=-visseur(1)+liste_thetaB_theta0(mod(theta0-1,360)+1);

    % Initial configuration save
    if theta0==k0
        C0=C;
        D0=D;
        M0=M;
        E0=E;
        Mp0=M_prime;
        k0=k0-17;
        A=[-L4,0];
        B=[-L4+L8*cos(thetaAB),L8*sin(thetaAB)];
    end
    
end

% Computation of the best angle to rotate the leg in order to keep a flat
% support phase of the trajectory
[~,indice1]=min(listeM_prime(:,1),[],'linear');
[~,indice2]=max(listeM_prime(:,2),[],'linear');
point1=[listeM_prime(indice1,:)];
point2=[listeM_prime(indice2,:)];
angle_opti=pi-atan((point1(2)-point2(2))/(point1(1)-point2(1)));

% Reorientation of the point trajectories considering the previous angle
listeM_prime_rot=Rotation(listeM_prime, angle_opti);
listeC_rot=Rotation(listeC, angle_opti);
listeD_rot=Rotation(listeD, angle_opti);
listeM_rot=Rotation(listeM, angle_opti);
listeE_rot=Rotation(listeE, angle_opti);
A_rot=Rotation(A, angle_opti);
B_rot=Rotation(B, angle_opti);
liste_pts=[listeC_rot; listeD_rot; listeM_rot; listeE_rot; A_rot; B_rot];
position_pas=min(listeM_prime_rot(:,2))-min(liste_pts(:,2));

% Animation of the leg motion
figure
i=1;
for theta0=0:1:359
    pause(0.01)
    clf
    patch([0 A_rot(1) B_rot(1) 0],[0 A_rot(2) B_rot(2) 0],'ko-','Markersize',5,'Linewidth',2,'FaceAlpha',0.25)
    hold on
    plot([A_rot(1) listeD_rot(i,1) listeC_rot(i,1) 0],[A_rot(2) listeD_rot(i,2) listeC_rot(i,2) 0],'ro-','Markersize',10,'Linewidth',2)
    patch([listeD_rot(i,1) listeM_rot(i,1) listeC_rot(i,1)],[listeD_rot(i,2) listeM_rot(i,2) listeC_rot(i,2)],'ro-','Markersize',10,'Linewidth',2,'FaceAlpha',0.5)
    plot([B_rot(1) listeE_rot(i,1) listeM_rot(i,1) listeM_prime_rot(i,1)],[B_rot(2) listeE_rot(i,2) listeM_rot(i,2) listeM_prime_rot(i,2)],'ro-','Markersize',10,'Linewidth',2)
    plot(listeM_prime_rot(1:i,1),listeM_prime_rot(1:i,2),'b','Linewidth',1)
    i=i+1;
    axis([-200 150 -250 50])
    title('Animation')
    grid on
end

% Original (before reorientation) optimal plots of the different point 
% trajectories
figure()
plot(listeM_prime(:,1),listeM_prime(:,2));
hold on
plot(listeD(:,1),listeD(:,2));
plot(listeM(:,1),listeM(:,2));
plot(listeE(:,1),listeE(:,2));
plot(listeC(:,1),listeC(:,2));
patch([0 A(1) B(1) 0],[0 A(2) B(2) 0],'ko-','Markersize',5,'Linewidth',2,'FaceAlpha',0.25)
plot([A(1) D0(1) C0(1) 0],[A(2) D0(2) C0(2) 0],'ro-','Markersize',10,'Linewidth',2)
patch([D0(1) M0(1) C0(1)],[D0(2) M0(2) C0(2)],'ro-','Markersize',10,'Linewidth',2,'FaceAlpha',0.5)
plot([B(1) E0(1) M0(1) Mp0(1)],[B(2) E0(2) M0(2) Mp0(2)],'ro-','Markersize',10,'Linewidth',2)
legend("M' Trajectory","D Trajectory","M Trajectory","E Trajectory","C Trajectory", "", "","","")
title('Point trajectories with optimal leg')
axis equal;
hold off

% Rotated trajectory plots without the leg
figure()
plot(listeM_prime_rot(:,1),listeM_prime_rot(:,2));
hold on
plot(listeD_rot(:,1),listeD_rot(:,2));
plot(listeM_rot(:,1),listeM_rot(:,2));
plot(listeE_rot(:,1),listeE_rot(:,2));
plot(listeC_rot(:,1),listeC_rot(:,2));
hold off
title('Reoriented point trajectories')

% Identification of the limits of the support/fly phases
figure()
plot(listeM_prime(:,1),listeM_prime(:,2),'b')
hold on
patch([0 A(1) B(1) 0],[0 A(2) B(2) 0],'ko-','Markersize',5,'Linewidth',2,'FaceAlpha',0.25)
plot([A(1) D0(1) C0(1) 0],[A(2) D0(2) C0(2) 0],'ro-','Markersize',10,'Linewidth',2)
patch([D0(1) M0(1) C0(1)],[D0(2) M0(2) C0(2)],'ro-','Markersize',10,'Linewidth',2,'FaceAlpha',0.5)
plot([B(1) E0(1) M0(1) Mp0(1)],[B(2) E0(2) M0(2) Mp0(2)],'ro-','Markersize',10,'Linewidth',2)
plot(point1(1),point1(2),'go-','Markersize',5,'Linewidth',2)
plot(point2(1),point2(2),'go-','Markersize',5,'Linewidth',2)
axis equal;
grid on
view(-90,90)
title('Fly/Support phase limits')

%Reoriented point trajectories with optimal leg
figure()
plot(listeM_prime_rot(:,1),listeM_prime_rot(:,2),'b')
hold on
plot(listeD_rot(:,1),listeD_rot(:,2));
plot(listeM_rot(:,1),listeM_rot(:,2));
plot(listeE_rot(:,1),listeE_rot(:,2));
plot(listeC_rot(:,1),listeC_rot(:,2));
liste=Rotation([[0, 0]; A; B; [0, 0]], angle_opti);
patch(liste(:,1),liste(:,2),'ko-','Markersize',5,'Linewidth',2,'FaceAlpha',0.25)
liste=Rotation([A; D0; C0; [0, 0]], angle_opti);
plot(liste(:,1),liste(:,2),'ro-','Markersize',10,'Linewidth',2)
liste=Rotation([D0; M0; C0], angle_opti);
patch(liste(:,1),liste(:,2),'ro-','Markersize',10,'Linewidth',2,'FaceAlpha',0.5)
liste=Rotation([B; E0; M0; Mp0], angle_opti);
plot(liste(:,1),liste(:,2),'ro-','Markersize',10,'Linewidth',2)
%plot([listeM_prime_rot(indice1,1)], [listeM_prime_rot(indice1,2)],'bo-','Markersize',5,'Linewidth',2)
%plot([listeM_prime_rot(indice2,1)], [listeM_prime_rot(indice2,2)],'go-','Markersize',5,'Linewidth',2)
legend("M' Trajectory","D Trajectory","M Trajectory","E Trajectory","C Trajectory", "", "","","")
title('Reoriented point trajectories with optimal leg')
axis equal;
grid on

% Display of the joint angles values
figure()
plot(linspace(0,359,360), liste_thetaA)
hold on 
plot(linspace(0,359,360), liste_thetaB)
plot(linspace(0,359,360), liste_thetaD)
plot(linspace(0,359,360), liste_thetaM)
plot(linspace(0,359,360), liste_thetaE)
legend('thetaA','thetaB','thetaD','thetaM','thetaE')
xlabel('theta0 (°)')
ylabel('Angles (°)')
grid on
title('Joint angle values in function of the input actuator angle')
hold off

% Foot trajectory
figure()
plot(listeM_prime_rot(:,1),listeM_prime_rot(:,2),'b')
title('Optimal foot trajectory')

% Computation of the required torque to drive the leg
for theta0=0:1:359
    % Compliance torque computation
    tD=CalculMoment(liste_thetaD(theta0+1)-liste_thetaD(k0));
    tE=CalculMoment(liste_thetaE(theta0+1)-liste_thetaE(k0));
    tM=CalculMoment(liste_thetaM(theta0+1)-liste_thetaM(k0));
    tA=CalculMoment(liste_thetaA(theta0+1)-liste_thetaA(k0));
    tB=CalculMoment(liste_thetaB(theta0+1)-liste_thetaB(k0));

    % Generated output foot force and actuator torque
    liste_force(theta0+1)=(-liste_thetaD_theta0(theta0+1)*tD-liste_thetaE_theta0(theta0+1)*tE-liste_thetaM_theta0(theta0+1)*tM-liste_thetaA_theta0(theta0+1)*tA-liste_thetaB_theta0(theta0+1)*tB)/liste_vitesse(theta0+1);
    liste_couple_a_vide(theta0+1)=-liste_thetaD_theta0(theta0+1)*tD-liste_thetaE_theta0(theta0+1)*tE-liste_thetaM_theta0(theta0+1)*tM-liste_thetaA_theta0(theta0+1)*tA-liste_thetaB_theta0(theta0+1)*tB;
    liste_energie_potentielle(theta0+1)=1/2*pi/180*(tD^2+tE^2+tM^2+tA^2+tB^2);
    
end


% Plot of the actuator torque
liste_couple_a_vide(1)=liste_couple_a_vide(2);
% Extract positive and negative part
yp = (liste_couple_a_vide + abs(liste_couple_a_vide))/2;
yn = (liste_couple_a_vide - abs(liste_couple_a_vide))/2;
% Plot the data using area function
figure()
plot(linspace(0,359,360), liste_couple_a_vide,'linewidth',2)
hold on 
area(linspace(0,359,360),yp,'FaceColor',[160/255 251/255 142/255])
area(linspace(0,359,360),yn,'FaceColor',[239/255 134/255 131/255])
plot([indice1], [liste_couple_a_vide(indice1)],'bo-','Markersize',15,'Linewidth',2)
plot([indice2], [liste_couple_a_vide(indice2)],'mo-','Markersize',15,'Linewidth',2)
legend("Motor torque",'Energy Release','Energy Storage','Support Phase start','Support Phase end','Location','SouthEast')
xlabel("theta0 (°)")
ylabel("Torque (N.m)")
hold off

%Print potential energy
figure()
plot(linspace(0,359,360), liste_energie_potentielle,'Linewidth',2)
hold on 
plot([indice1], [liste_energie_potentielle(indice1)],'bo-','Markersize',15,'Linewidth',2)
plot([indice2], [liste_energie_potentielle(indice2)],'mo-','Markersize',15,'Linewidth',2)
legend("Potential Energy",'Support Phase start','Support Phase end')
xlabel("theta0 (°)")
ylabel("Potential Energy (J)")
hold off


%%% SUBFUNCTIONS 

function Moment = CalculMoment(angle)
    % Computation of the torque due to joint stiffness
    Moment=-5.9479*10^-6*angle^2+1.2386*10^-3*angle;
end

function liste_pts_rot = Rotation(liste_pts, angle)
    %Computation of a rotated position list
    liste_pts_rot=zeros(size(liste_pts,[1]),2);
    for k=1:size(liste_pts,[1])
        liste_pts_rot(k,:)=[liste_pts(k,1)*cos(angle)-liste_pts(k,2)*sin(angle); liste_pts(k,1)*sin(angle)+liste_pts(k,2)*cos(angle)];
    end

end

function vect = vecteur(pt2x, pt2y, pt1x, pt1y, norme)
    %Unit vector between two points
    vect = [pt2x - pt1x, pt2y - pt1y];
    norme_vect = sqrt(vect(1)^2 + vect(2)^2);
    if norme == 1
        vect = [vect(1)/norme_vect vect(2)/norme_vect];
    end
end

function scalaire = prod_scal(vect1, vect2)
    %Dot product between two vectors
    global mat_E
    scalaire=vect1*mat_E*vect2';
end
