%Limpieza de pantalla
clear all
close all
clc
tic
%Declaración de variables simbólicas
syms th1(t) th2(t) l3(t) th4(t) th5(t) th6(t) t %Articulación
syms th1p(t) th2p(t) l3p(t) th4p(t) th5p(t) th6p(t)   %Velocidades de cada articulación
syms th1pp(t) th2pp(t) l3pp(t) th4pp(t) th5pp(t) th6pp(t) %Aceleraciones de cada articulación
syms m1 m2 m3 m4 m5 m6 Ixx1 Iyy1 Izz1 Ixx2 Iyy2 Izz2 Ixx3 Iyy3 Izz3 Ixx4 Iyy4 Izz4 Ixx5 Iyy5 Izz5 Ixx6 Iyy6 Izz6 %Masas y matrices de Inercia
syms l1 l2 l4 l5 l6 lc1 lc2 lc3 lc4 lc5 lc6 %l=longitud de eslabones y lc=distancia al centro de masa de cada eslabón
syms pi g a cero
%Creamos el vector de coordenadas articulares
Q= [th1; th2; l3; th4; th5; th6];
%disp('Coordenadas generalizadas');
%pretty (Q);
%Creamos el vector de velocidades articulares
Qp= [th1p; th2p; l3p; th4p; th5p; th6p];
%disp('Velocidades generalizadas');
%pretty (Qp);
%Creamos el vector de aceleraciones articulares
Qpp= [th1pp; th2pp; l3pp; th4pp; th5pp; th6pp];
%disp('Aceleraciones generalizadas');
%pretty (Qpp);
%Configuración del robot, 0 para junta rotacional, 1 para junta prismática
RP=[0 0 1 0 0 0];
%Número de grado de libertad del robot
GDL= size(RP,2);
GDL_str= num2str(GDL);
%Articulación 1 
%Posición de la articulación
P(:,:,1)= [0; 0; l1];
%Matriz de rotación 
Rz1= [cos(th1) -sin(th1)  0;
      sin(th1)  cos(th1)  0;
      0         0         1];
Rx1_90= [1 0 0;
         0 0 -1;
         0 1 0];
R(:,:,1)= Rz1 * Rx1_90;
%Articulación 2
%Posición de la articulación
P(:,:,2)= [0; 0; -l2];
%Matriz de rotación 
Rz2= [cos(th2) -sin(th2)  0;
      sin(th2)  cos(th2)  0;
      0         0         1];
Rx2_90n= [1 0 0;
          0 0 1;
          0 -1 0];
R(:,:,2)= Rz2 * Rx2_90n;
%Articulación 3
%Posición de la articulación
P(:,:,3)= [0; 0; l3];
%Matriz de rotación
Identidad3= [1 0 0;
             0 1 0;
             0 0 1];
R(:,:,3)= Identidad3;
%Articulación 4
%Posición de la articulación
P(:,:,4)= [0; 0; l4];
%Matriz de rotación
Rz4= [cos(th4) -sin(th4)  0;
      sin(th4)  cos(th4)  0;
      0         0         1];
Rx4_90= [1 0 0;
         0 0 -1;
         0 1 0];
R(:,:,4)= Rz4 * Rx4_90;
%Articulación 5
%Posición de la articulación
P(:,:,5)= [-l5*cos(th5); l5*sin(th5); 0];
%Matriz de rotación
Rz5= [cos(th5) -sin(th5)  0;
      sin(th5)  cos(th5)  0;
      0         0         1];
Rx5_90n= [1 0 0;
          0 0 1;
          0 -1 0];
R(:,:,5)= Rz5 * Rx5_90n;
%Articulación 6
%Posición de la articulación
P(:,:,6)= [0; 0; l6];
%Matriz de rotación
Identidad6= [1 0 0;
             0 1 0;
             0 0 1];
R(:,:,6)= Identidad6;
%Creamos un vector de ceros
Vector_Zeros= zeros(1, 3);
%Inicializamos las matrices de transformación Homogénea locales
A(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
%Inicializamos las matrices de transformación Homogénea globales
T(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
%Inicializamos las posiciones vistas desde el marco de referencia inercial
PO(:,:,GDL)= P(:,:,GDL); 
%Inicializamos las matrices de rotación vistas desde el marco de referencia inercial
RO(:,:,GDL)= R(:,:,GDL); 
for i = 1:GDL
    i_str= num2str(i);
   %disp(strcat('Matriz de Transformación local A', i_str));
    A(:,:,i)=simplify([R(:,:,i) P(:,:,i); Vector_Zeros 1]);
   %pretty (A(:,:,i));
   %Globales
    try
       T(:,:,i)= T(:,:,i-1)*A(:,:,i);
    catch
       T(:,:,i)= A(:,:,i);
    end
%     disp(strcat('Matriz de Transformación global T', i_str));
    T(:,:,i)= simplify(T(:,:,i));
%     pretty(T(:,:,i))
    RO(:,:,i)= T(1:3,1:3,i);
    PO(:,:,i)= T(1:3,4,i);
    %pretty(RO(:,:,i));
    %pretty(PO(:,:,i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CALCULAMOS LAS VELOCIDADES PARA CADA ESLABÓN%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Jacobiano eslabón 1
Jv_a1(:,1)=PO(:,:,1);
Jw_a1(:,1)=PO(:,:,1);
for k=1:1
    if RP(k)==0
        try
            Jv_a1(:,k)=cross(RO(:,3,k-1),PO(:,:,1)-PO(:,:,k-1));
            Jw_a1(:,k)=RO(:,3,k-1);
        catch
            Jv_a1(:,k)=cross([0;0;1],PO(:,:,1));
            Jw_a1(:,k)=[0;0;1];
        end
    else
        try
            Jv_a1(:,k)=RO(:,3,k-1);
        catch
            Jv_a1(:,k)=[0;0;1];
        end
        Jw_a1(:,k)=[0;0;0];
    end
end
Jv_a1=simplify(Jv_a1);
Jw_a1=simplify(Jw_a1);
Jac1=[Jv_a1;Jw_a1];
disp('Velocidad lineal eslabón 1')
Qp=Qp(t);
V1=simplify(Jv_a1*Qp(1:1));
pretty(V1)
disp('Velocidad angular eslabón 1')
W1=simplify(Jw_a1*Qp(1:1));
pretty(W1)
%Jacobiano eslabón 2
Jv_a2(:,2)=PO(:,:,2);
Jw_a2(:,2)=PO(:,:,2);
for k=1:2
    if RP(k)==0
        try
            Jv_a2(:,k)=cross(RO(:,3,k-1),PO(:,:,2)-PO(:,:,k-1));
            Jw_a2(:,k)=RO(:,3,k-1);
        catch
            Jv_a2(:,k)=cross([0;0;1],PO(:,:,2));
            Jw_a2(:,k)=[0;0;1];
        end
    else
        try
            Jv_a2(:,k)=RO(:,3,k-1);
        catch
            Jv_a2(:,k)=[0;0;1];
        end
        Jw_a2(:,k)=[0;0;0];
    end
end
Jv_a2=simplify(Jv_a2);
Jw_a2=simplify(Jw_a2);
Jac2=[Jv_a2;Jw_a2];
disp('Velocidad lineal eslabón 2')
V2=simplify(Jv_a2*Qp(1:2));
pretty(V2)
disp('Velocidad angular eslabón 2')
W2=simplify(Jw_a2*Qp(1:2));
pretty(W2)
%Jacobiano eslabón 3
Jv_a3(:,3)=PO(:,:,3);
Jw_a3(:,3)=PO(:,:,3);
for k=1:3
    if RP(k)==0
        try
            Jv_a3(:,k)=cross(RO(:,3,k-1),PO(:,:,3)-PO(:,:,k-1));
            Jw_a3(:,k)=RO(:,3,k-1);
        catch
            Jv_a3(:,k)=cross([0;0;1],PO(:,:,3));
            Jw_a3(:,k)=[0;0;1];
        end
    else
        try
            Jv_a3(:,k)=RO(:,3,k-1);
        catch
            Jv_a3(:,k)=[0;0;1];
        end
        Jw_a3(:,k)=[0;0;0];
    end
end
Jv_a3=simplify(Jv_a3);
Jw_a3=simplify(Jw_a3);
Jac3=[Jv_a3;Jw_a3];
disp('Velocidad lineal eslabón 3')
V3=simplify(Jv_a3*Qp(1:3));
pretty(V3)
disp('Velocidad angular eslabón 3')
W3=simplify(Jw_a3*Qp(1:3));
pretty(W3)
%Jacobiano eslabón 4
Jv_a4(:,4)=PO(:,:,4);
Jw_a4(:,4)=PO(:,:,4);
for k=1:4
    if RP(k)==0
        try
            Jv_a4(:,k)=cross(RO(:,3,k-1),PO(:,:,4)-PO(:,:,k-1));
            Jw_a4(:,k)=RO(:,3,k-1);
        catch
            Jv_a4(:,k)=cross([0;0;1],PO(:,:,4));
            Jw_a4(:,k)=[0;0;1];
        end
    else
        try
            Jv_a4(:,k)=RO(:,3,k-1);
        catch
            Jv_a4(:,k)=[0;0;1];
        end
        Jw_a4(:,k)=[0;0;0];
    end
end
Jv_a4=simplify(Jv_a4);
Jw_a4=simplify(Jw_a4);
Jac4=[Jv_a4;Jw_a4];
disp('Velocidad lineal eslabón 4')
V4=simplify(Jv_a4*Qp(1:4));
pretty(V4)
disp('Velocidad angular eslabón 4')
W4=simplify(Jw_a4*Qp(1:4));
pretty(W4)
%Jacobiano eslabón 5
Jv_a5(:,5)=PO(:,:,5);
Jw_a5(:,5)=PO(:,:,5);
for k=1:5
    if RP(k)==0
        try
            Jv_a5(:,k)=cross(RO(:,3,k-1),PO(:,:,5)-PO(:,:,k-1));
            Jw_a5(:,k)=RO(:,3,k-1);
        catch
            Jv_a5(:,k)=cross([0;0;1],PO(:,:,5));
            Jw_a5(:,k)=[0;0;1];
        end
    else
        try
            Jv_a5(:,k)=RO(:,3,k-1);
        catch
            Jv_a5(:,k)=[0;0;1];
        end
        Jw_a5(:,k)=[0;0;0];
    end
end
Jv_a5=simplify(Jv_a5);
Jw_a5=simplify(Jw_a5);
Jac5=[Jv_a5;Jw_a5];
disp('Velocidad lineal eslabón 5')
V5=simplify(Jv_a5*Qp(1:5));
pretty(V5)
disp('Velocidad angular eslabón 5')
W5=simplify(Jw_a5*Qp(1:5));
pretty(W5)
%Jacobiano eslabón 6
Jv_a6(:,6)=PO(:,:,6);
Jw_a6(:,6)=PO(:,:,6);
for k=1:6
    if RP(k)==0
        try
            Jv_a6(:,k)=cross(RO(:,3,k-1),PO(:,:,6)-PO(:,:,k-1));
            Jw_a6(:,k)=RO(:,3,k-1);
        catch
            Jv_a6(:,k)=cross([0;0;1],PO(:,:,6));
            Jw_a6(:,k)=[0;0;1];
        end
    else
        try
            Jv_a6(:,k)=RO(:,3,k-1);
        catch
            Jv_a6(:,k)=[0;0;1];
        end
        Jw_a6(:,k)=[0;0;0];
    end
end
Jv_a6=simplify(Jv_a6);
Jw_a6=simplify(Jw_a6);
Jac6=[Jv_a6;Jw_a6];
disp('Velocidad lineal eslabón 6')
V6=simplify(Jv_a6*Qp(1:6));
pretty(V6)
disp('Velocidad angular eslabón 6')
W6=simplify(Jw_a6*Qp(1:6));
pretty(W6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Energía Cinética
%%%%%%%%%%%%%%%%%%%%%%%%%%Omitimos la división de cada lc%%%%%%%%%%%%%%%
%Distancia del origen del eslabón a su centro de masa
%Vectores de posición respecto al centro de masa
P01=subs(P(:,:,1)/2, l1, lc1);
P12=subs(P(:,:,2)/2, l2, lc2);
P23=subs(P(:,:,3)/2, l3, lc3);
P34=subs(P(:,:,4)/2, l4, lc4);
P45=subs(P(:,:,5)/2, l5, lc5);
P56=subs(P(:,:,6)/2, l6, lc6);
%Creamos matrices de inercia para cada eslabón
I1=[Ixx1 0 0; 
    0 Iyy1 0; 
    0 0 Izz1];
I2=[Ixx2 0 0; 
    0 Iyy2 0; 
    0 0 Izz2];
I3=[Ixx3 0 0; 
    0 Iyy3 0; 
    0 0 Izz3];
I4=[Ixx4 0 0; 
    0 Iyy4 0; 
    0 0 Izz4];
I5=[Ixx5 0 0; 
    0 Iyy5 0; 
    0 0 Izz5];
I6=[Ixx6 0 0; 
    0 Iyy6 0; 
    0 0 Izz6];
%Calculamos la energía cinética para cada uno de los eslabones
%Eslabón 1
V1_Total= V1+cross(W1,P01);
K1= (1/2*m1*(V1_Total))'*((V1_Total)) + (1/2*W1)'*(I1*W1);
disp('Energía Cinética en el Eslabón 1');
K1= simplify(K1);
pretty(K1);
%Eslabón 2
V2_Total= V2+cross(W2,P12);
K2= (1/2*m2*(V2_Total))'*((V2_Total)) + (1/2*W2)'*(I2*W2);
disp('Energía Cinética en el Eslabón 2');
K2= simplify(K2);
pretty(K2);
%Eslabón 3
V3_Total= V3+cross(W3,P23);
K3= (1/2*m3*(V3_Total))'*((V3_Total)) + (1/2*W3)'*(I3*W3);
disp('Energía Cinética en el Eslabón 3');
K3= simplify(K3);
pretty(K3);
%Eslabón 4
V4_Total= V4+cross(W4,P34);
K4= (1/2*m4*(V4_Total))'*((V4_Total)) + (1/2*W4)'*(I4*W4);
disp('Energía Cinética en el Eslabón 4');
K4= simplify(K4);
pretty(K4);
%Eslabón 5
V5_Total= V5+cross(W5,P45);
K5= (1/2*m5*(V5_Total))'*((V5_Total)) + (1/2*W5)'*(I5*W5);
disp('Energía Cinética en el Eslabón 5');
K5= simplify(K5);
pretty(K5);
%Eslabón 6
V6_Total= V6+cross(W6,P56);
K6= (1/2*m6*(V6_Total))'*((V6_Total)) + (1/2*W6)'*(I6*W6);
disp('Energía Cinética en el Eslabón 6');
K6= simplify(K6);
pretty(K6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VECTORES A LOS CENTROS DE MASA (CoM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Vectores locales (Sin el /2)
P01_local = subs(P(:,:,1), l1, lc1);
P12_local = subs(P(:,:,2), l2, lc2);
P23_local = subs(P(:,:,3), l3, lc3);
P34_local = subs(P(:,:,4), l4, lc4);
P45_local = subs(P(:,:,5), l5, lc5);
P56_local = subs(P(:,:,6), l6, lc6);
% 2. Rotación de vectores locales al marco global
P01_g = P01_local; 
P12_g = RO(:,:,1) * P12_local;
P23_g = RO(:,:,2) * P23_local;
P34_g = RO(:,:,3) * P34_local;
P45_g = RO(:,:,4) * P45_local;
P56_g = RO(:,:,5) * P56_local;
% 3. Posiciones Globales Absolutas de cada CoM (Vital para la Energía Potencial)
Pos_cm1 = P01_g;
Pos_cm2 = PO(:,:,1) + P12_g;
Pos_cm3 = PO(:,:,2) + P23_g;
Pos_cm4 = PO(:,:,3) + P34_g;
Pos_cm5 = PO(:,:,4) + P45_g;
Pos_cm6 = PO(:,:,5) + P56_g;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENERGÍA CINÉTICA (K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrices de inercia locales
I1=[Ixx1 0 0; 0 Iyy1 0; 0 0 Izz1];
I2=[Ixx2 0 0; 0 Iyy2 0; 0 0 Izz2];
I3=[Ixx3 0 0; 0 Iyy3 0; 0 0 Izz3];
I4=[Ixx4 0 0; 0 Iyy4 0; 0 0 Izz4];
I5=[Ixx5 0 0; 0 Iyy5 0; 0 0 Izz5];
I6=[Ixx6 0 0; 0 Iyy6 0; 0 0 Izz6];
% Rotación de Inercias al marco global
I1_g = RO(:,:,1) * I1 * RO(:,:,1).';
I2_g = RO(:,:,2) * I2 * RO(:,:,2).';
I3_g = RO(:,:,3) * I3 * RO(:,:,3).';
I4_g = RO(:,:,4) * I4 * RO(:,:,4).';
I5_g = RO(:,:,5) * I5 * RO(:,:,5).';
I6_g = RO(:,:,6) * I6 * RO(:,:,6).';
% Velocidades Reales de los Centros de Masa (Parten de la junta ANTERIOR)
Vc1 = cross(W1, P01_g);
Vc2 = V1 + cross(W2, P12_g);
Vc3 = V2 + cross(W3, P23_g);
Vc4 = V3 + cross(W4, P34_g);
Vc5 = V4 + cross(W5, P45_g);
Vc6 = V5 + cross(W6, P56_g);
disp('Calculando Energías Cinéticas (Puede tardar por los 6 GDL)...');
K1 = simplify( 0.5*m1*(Vc1.' * Vc1) + 0.5*W1.' * I1_g * W1 );
K2 = simplify( 0.5*m2*(Vc2.' * Vc2) + 0.5*W2.' * I2_g * W2 );
K3 = simplify( 0.5*m3*(Vc3.' * Vc3) + 0.5*W3.' * I3_g * W3 );
K4 = simplify( 0.5*m4*(Vc4.' * Vc4) + 0.5*W4.' * I4_g * W4 );
K5 = simplify( 0.5*m5*(Vc5.' * Vc5) + 0.5*W5.' * I5_g * W5 );
K6 = simplify( 0.5*m6*(Vc6.' * Vc6) + 0.5*W6.' * I6_g * W6 );
K_Total = simplify(K1 + K2 + K3 + K4 + K5 + K6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENERGÍA POTENCIAL (U)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extraemos la altura (eje Y global, índice 2) de las posiciones absolutas
U1 = m1 * g * Pos_cm1(2);
U2 = m2 * g * Pos_cm2(2);
U3 = m3 * g * Pos_cm3(2);
U4 = m4 * g * Pos_cm4(2);
U5 = m5 * g * Pos_cm5(2);
U6 = m6 * g * Pos_cm6(2);
U_Total = simplify(U1 + U2 + U3 + U4 + U5 + U6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODELO DINÁMICO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Construyendo Lagrangiano y Hamiltoniano...');
Lagrangiano = simplify(K_Total - U_Total);
H = simplify(K_Total + U_Total);
disp('====================================================');
disp('¡Cálculos Finalizados!');
disp('Las variables Lagrangiano y H están en el Workspace.');
disp('====================================================');
toc
