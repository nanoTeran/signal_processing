%Design of a electromagnet

%% Iniciales
% Constantes siempre
C_lamina = 1500;
grosor_l = 1; % (mm)
densidad_cu = 8960; % kg/m3
Uo = 4*pi*10^(-7);
g = 9.8;
Ur = 4000;
L_n = (65+100)/1000;
L_a = 0.002;
C_kg_cu =100000;

% Constantes variables
Acap = 3;
m = 2; %masa kg


% V. Independientes
syms prof % Profundidad nucleo (mm)
syms AWG % Profundidad nucleo (mm)
% prof = 20;
% AWG = 36;

% V. Dependientes
F = m*g;
A = prof *80;
B = sqrt(2*F*Uo*(10^6)/A);
N_laminas = prof*grosor_l;
diametro = 2*((92^(((36-AWG)/39)))*0.0635)*0.001;
area_t = pi*((diametro/2)^2);
I = Acap*area_t*10^6;
NI = B*((L_n/Ur)+L_a)/Uo;
N = NI/I;
L = 2*N*(0.04+prof/1000);

m_cu =densidad_cu*area_t*L;

%% Costos
yC_n = C_lamina*N_laminas; % nucleo
yC_c = C_kg_cu*m_cu; % cobre
yC = yC_c+yC_n; % total

%% Limites de graficas
prof1 = 1;
prof2 = 100;
AWG1 = 14;
AWG2 = 36;
%% Gráficar costos VS profundidad
xC = prof; %Variable independiente

figure()
fplot(xC,yC_c,[prof1 prof2])
title('Costos cable');
xlabel('Profundidad (mm)');
ylabel('Costo cable ($)');

%% Gráficar costos VS profundidad
xC = prof; %Variable independiente

figure()
subplot(2,2,1)
fplot(xC,yC_n,[prof1 prof2])
title('Costos Núcleo');
xlabel('Profundidad (mm)');
ylabel('Costo núcleo ($)');

subplot(2,2,2)
fplot(xC,yC_c,[prof1 prof2])
title('Costos cable');
xlabel('Profundidad (mm)');
ylabel('Costo cable ($)');

subplot(2,2,3)
fplot(xC,yC,[prof1 prof2])
title('Costos Total');
xlabel('Profundidad (mm)');
ylabel('Costo Total ($)');



%% Graficas (prof, AWG)
x = prof;
y = AWG;

f1 = yC_n;
f2 = yC_c;
f3 = yC;

figure()
fsurf(f1, [AWG1 AWG2 prof1 prof2])
title('Costo núcleo(prof, AWG)')
ylabel('Profundidad (mm)');
xlabel('AWG');
zlabel('Costo núcleo ($)');

figure()
fsurf(f2, [AWG1 AWG2 prof1 prof2])
title('Costo cable (prof, AWG)')
ylabel('Profundidad (mm)');
xlabel('AWG');
zlabel('Costo cable ($)');

figure()
fsurf(f3, [AWG1 AWG2 prof1 prof2])
title('Costo total (prof, AWG)')
ylabel('Profundidad (mm)');
xlabel('AWG');
zlabel('Costo Total ($)');

%%
f4 = N;
figure()
subplot(1,2,1)
ylabel('y');
xlabel('x');
zlabel('z');

fsurf(f4, [prof1 prof2 AWG1 AWG2])

subplot(1,2,2)
fsurf(f4, [AWG1 AWG2 prof1 prof2])
title('f4(prof, AWG)');


