clear 
clc
close all

global C

%Gotta include the constants just in case!
C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      %metres (32.1740 ft) per s


%Change to change area of calculation
nx = 200;       %length of boxed area
ny = 100;       %width of boxed area
lengthBox = 5;     %length of bottle neck box
heightBox = 2;     %height of bottle neck box



G = zeros(nx*ny, nx*ny);        %Create the G matrix to be used
B = zeros(nx*ny,1);             %Boundary condition vector
sigma = ones(nx,ny);           %Conductivity Matrix
% sigmaOutOfBox = 1;              %Free space conductivity
sigmaInBox = .01;               %Box conductivity
voltagePlot1 = zeros(nx,ny);    %Voltage matrix pre-initialization to save time



%Old Stuff, don't need this


% %Set the area for the calculations and plot a visualization of the region
% figure(1)
% %Bottom Box
% rectangle ('position', [(nx/2-lengthBox/2)  1 lengthBox heightBox],'FaceColor',[.5 .5 .5],'EdgeColor','k');
% %Top Box
% rectangle ('position', [(nx/2-lengthBox/2)  (ny-heightBox) lengthBox heightBox],'FaceColor',[.5 .5 .5],'EdgeColor','k');
% grid on;
% xlim([1 nx]);
% ylim ([1 ny]);
% title ('Area for Current Calculations') 
box1 = [80 120  1  30];
box2 = [80 120 70 100];

sigma(box1(1):box1(2), box1(3):box1(4)) = sigmaInBox;
sigma(box2(1):box2(2), box2(3):box2(4)) = sigmaInBox;

%Used to calculate sthe N value for the G matrix
fn = @(i, j) j + (i-1)*ny;

for w = 1:nx % x
    for s = 1:ny % y
        %Calculate the n values to be used for G matrix
        n = fn(w, s);
        nxm = fn(w-1, s);
        nxp = fn(w+1, s);
        nym = fn(w, s-1);
        nyp = fn(w, s+1);
        if w == 1
            G(n, n) = 1;
            B(n) = 0.1;
        elseif w == nx
            G(n, n) = 1;
            B(n) = 0;
        elseif s == 1
            s1 = (sigma(w,s)+sigma(w+1,s))/1.0;
            s2 = (sigma(w,s)+sigma(w-1,s))/1.0;
            s3 = (sigma(w,s)+sigma(w,s+1))/2.0;
            
            G(n, n) = -(s1+s2+s3);
            G(n, nxp) = s1;
            G(n, nxm) = s2;
            G(n, nyp) = s3;
        elseif s == ny
            s1 = (sigma(w,s)+sigma(w+1,s))/1.0;
            s2 = (sigma(w,s)+sigma(w-1,s))/1.0;
            s3 = (sigma(w,s)+sigma(w,s-1))/2.0;
            
            G(n, n) = -(s1 + s2 + sym);
            G(n, nxp) = s1;
            G(n, nxm) = s2;
            G(n, nym) = sym;
        else
            s2 = (sigma(w,s)+sigma(w-1,s))/2.0;
            s1 = (sigma(w,s)+sigma(w+1,s))/2.0;
            sym = (sigma(w,s)+sigma(w,s-1))/2.0;
            s3 = (sigma(w,s)+sigma(w,s+1))/2.0;
            
            G(n, n) = -(s2+s1+sym+s3);
            G(n, nxp) = s1;
            G(n, nxm) = s2;
            G(n, nyp) = s3;
            G(n, nym) = sym;
        end
    end
end

figure(2)
surf(sigma)
title('Visualization Plot of Sigma')

%Calculate the voltage matrix from the boundary and G
V=G\B;

%Remap V back to a familiar matrix
for w = 1:nx % x
    for s = 1:ny % y
        n = s+(w-1)*ny;
        voltagePlot1(w, s) = V(n);
    end 
end

%Plot the voltage matrix
figure(3)
surf(voltagePlot1)
title('Voltage Spreading Through Bottleneck');




%My old plots
%Calculate and plot the EX electrix field
figure(4)
[ex,ey]=gradient(voltagePlot1);
quiver(ex.', ey.');

title('Electric Field Ex');
% 
% %Calculate and plot the EY electrix field
% figure(5)
% surf(ey)
% title('Electric Field Ey');
% 
% %Calculate current density and plot
% J = sigma .* gradient(voltagePlot1);
% figure (6)
% surf(J)
% title('Current Density')
