% clear all
clearvars
clear
clc
clearvars -GLOBAL
close all
format shorte

global C

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      %metres (32.1740 ft) per s²
C.am = 1.66053892e-27;

TotalElectrons = 1000;
lastTemperature= 300;
Tau = 0.2e-12;
length = 200e-9;
height = 100e-9;
Px= zeros(1,TotalElectrons);
NewPx= zeros(1,TotalElectrons);
Py= zeros(1,TotalElectrons);
NewPy= zeros(1,TotalElectrons);
Vx= zeros(1,TotalElectrons);
Vy= zeros(1,TotalElectrons);
temperatures = zeros(11,11);
lastLoop=0;
lastDensity=0;


xVolts = 0.1;
yVolts = 0;

xField = xVolts/length;
yField = yVolts/height;
fprintf('The x Field is %d\n', xField);
fprintf('The y Force is %d\n', yField);

xForce = xField*C.q_0;
yForce = yField*C.q_0;

fprintf('The x Force is %d\n', xForce);
fprintf('The y Force is %d\n', yForce);

xAcc = xForce/(C.m_0*.26);
yAcc = yForce/(C.m_0*.26);

fprintf('The x Acceleration is %d\n', xAcc);
fprintf('The y Acceleration is %d\n', yAcc);

currentConsentration = 1e15;
currentConsentration2 = currentConsentration/1e-4;


Vth = sqrt((C.kb*lastTemperature)/(C.m_0*0.26));
time = 150;
dt = 15e-15;
lambda = Vth*Tau;
fprintf('The thermal velocity is equal to %d\n',Vth)
fprintf('The mean free path is equal to %d\n',lambda)


figure('Name','Electron Paths')

Px = rand(1, TotalElectrons)*200e-9;       %Generate the random x location
Py = rand(1, TotalElectrons)*100e-9;       %generate the random y location

RandAng = rand(1, TotalElectrons)*2*pi;
Vx = Vth*sin(RandAng);
Vy = Vth*cos(RandAng);

hold on

for loops=1:time
    NewPx = Vx*dt+Px;
    NewPy = Vy*dt+Py;
    
    Vx = Vx + dt*xAcc;
    Vx = Vx + dt*yAcc;
    
    ix = NewPx>length;
    NewPx(ix) = NewPx(ix)-length;
    Px(ix) = Px(ix)-length;
    
    ix = NewPx<0;
    NewPx(ix) = NewPx(ix)+length;
    Px(ix) = Px(ix)+length;

    ix = NewPy<0;
    Vy(ix) = -Vy(ix);
    
    ix = NewPy>height;
    Vy(ix) = -Vy(ix);
    
    subplot (2,3,1)
    plot([Px(1) NewPx(1)], [Py(1) NewPy(1)], 'b')
    plot([Px(2) NewPx(2)], [Py(2) NewPy(2)], 'g')
    plot([Px(3) NewPx(3)], [Py(3) NewPy(3)], 'r')
    plot([Px(4) NewPx(4)], [Py(4) NewPy(4)], 'c')
    plot([Px(5) NewPx(5)], [Py(5) NewPy(5)], 'm')
    plot([Px(6) NewPx(6)], [Py(6) NewPy(6)], 'y')
    hold on
    title('Electron Paths');
    xlim([0 200e-9]);
    ylim ([0 100e-9]);
    pause(0.000001)
    
    Px=NewPx;
    Py=NewPy;

    % Drift Velocity to Temperature
    V(1, :) = sqrt(Vx(1, :).^2 + Vy(1, :).^2);
    V_mean = mean(V.^2);
    temperature = V_mean*(C.m_0*.26)/C.kb;

    % Plot Temperature Lines
    subplot(2,3,2)
    plot([lastLoop loops], [lastTemperature temperature], 'r');
    title('Temperature');
    xlabel('Time Step'); ylabel('Temperature (K)');
    xlim([1 time]);
    hold on;
    
    subplot(2,3,5)
    hist3([Px', Py'], 'CDataMode', 'auto', 'FaceColor', 'interp')
    colormap('default');
    colorbar;
    xlabel('Length'); ylabel('Height');
    title('Electron Density');
    view(2);

    
    %Map the temps
    for y = 1:10
        ymax = y*10;
        ymin = ymax-10;
        for x = 1:10
            xmax = x*20;
            xmin = xmax-20;
            side1 = Px > (xmin*1e-9);
            side2 = Px < (xmax*1e-9);
            side3 = Py > (ymin*1e-9);
            side4 = Py < (ymax*1e-9);
            between1 = bitand(side1, side2);
            between2 = bitand(side3, side4);
            hit = bitand(between1, between2);
            velocity = sqrt((Vx(hit).^2) + (Vy(hit).^2));
            v_mean = mean(velocity.*velocity);
            temperature_value = (((v_mean)*(C.m_0*.26))/(C.kb));
            temperatures(x, y) = temperature_value;
        end       
    end
    

    %Distribut the temps
    subplot(2,3,4)
    surf(transpose(temperatures));
    xlabel('Length');
    ylabel('Height');
    title('Temperature Distribution');
    colorbar;
    view(2);

     
     
     
     
     
     % Electron Drift Current Density
    drift = mean(Vx);
    density = C.q_0*currentConsentration*drift;
    
    % Plotting Current Density
    subplot(2,3,3)
    plot([lastLoop loops], [lastDensity density], 'r');
    title('Current Density');
    xlabel('Time Step'); ylabel('Current Density (A/cm2)');
    xlim([1 time]);
    hold on;

    pause(0.01);
    
    % Update electron coordinates
    ;
    lastLoop = loops;
    lastTemperature = temperature;
    lastDensity = density;
end

