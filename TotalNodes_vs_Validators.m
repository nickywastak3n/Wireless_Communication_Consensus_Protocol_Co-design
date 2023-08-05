clc;
close all;
clear;

%==============================
% Definitions
transmitPower_g = 2.5e-3; % Transmit power for gossip communication
transmitPower_b = 100e-3; % Transmit power for broadcast communication
noisePower = 1e-13; % Noise power
ro = 10; % Target SNR

R_0 = 1;
lambda = 0.125;
PL_R_0_db = (-20)*log10(lambda/(4*pi*R_0));
eta = 3; % Path-loss exponent
  
N_vector = (1:1:20).^2-1; % Number of nodes in system

alpha = 0.99;
beta = 1; % Acceptable consensus distortion
gamma = 0.9;

zeta = 0.9999; % Target dissemination success probability

a = 0.14; % Hermite-Pade approximation constant
g_inverse_positive = @(x) sqrt(-2/pi/a-log(1-x.^2)/2+sqrt((2/pi/a+log(1-x.^2)./2).^2-log(1-x.^2)./a)); % g inverse for x>0
g_inverse = @(x) (x>=0).*g_inverse_positive(x)+(x<0)*(-1).*g_inverse_positive(x); % g inverse for all x

%==============================
R2C_g_num = zeros(length(N_vector),1);
R2C_b_num = zeros(length(N_vector),1);

for n=1:length(N_vector)
    N = N_vector(n);
    F = ceil(N*0.1);

    R = 100/(sqrt(N+1)-1); % Distance between two neighboring nodes

    % Outage probability matrix
    outage_probability = @(x1,y1,x2,y2) 1-exp((-1)*10^(PL_R_0_db/10)*ro*noisePower/transmitPower_b*(sqrt((x1-x2)^2+(y1-y2)^2)*R/R_0)^eta);
    outage_probability_matrix = zeros(N+1,N+1);
    for i=1:N+1
        for k=1:N+1
            if i==k
                continue;
            end

            x1 = floor((i-1)/sqrt(N+1))+1;
            y1 = i-(x1-1)*sqrt(N+1);
            x2 = floor((k-1) /sqrt(N+1))+1;
            y2 = k-(x2-1)*sqrt(N+1);

            outage_probability_matrix(i,k) = outage_probability(x1,y1,x2,y2);
        end
    end

    A = (1/3)-(F/N);
    B = (F*(N-F)/((N-1)*N^2))*(g_inverse(2*alpha-1)^2);
    phi = 0.5; % Correction factor for normal approximation
    N_alpha = (phi*A+B*N+sqrt(2*phi*A*B*N-2*phi^2*B+B^2*N^2))/(A^2+2*B); 
    
    % psi_gossip_center = (13*N^2-4*N-8)*N/24/(N-1);
    % N_gamma_gossip_center = 1/(1/N+beta^2*N/2/g_inverse(gamma)^2/psi_gossip_center);
    psi_gossip_corner = (N+1)*((13*N-24*sqrt(N+1)+16)*N+12*(sqrt(N+1)-1))/6/(N-1);
    N_gamma_gossip_corner = 1/(1/N+beta^2*N/2/g_inverse(gamma)^2/psi_gossip_corner);
    
    psi_broadcast_corner = 0;
    p = 1;
    for v=1:N+1
        if v ~= p
            psi_broadcast_corner = psi_broadcast_corner + (1+outage_probability_matrix(p,v))/(1-outage_probability_matrix(p,v))^2;
        
            temp = 0;
            for j=1:N+1
                if j ~= p && j~= v
                    temp = temp + 1/(1-outage_probability_matrix(p,v))/(1-outage_probability_matrix(p,j));
                end
            end
            psi_broadcast_corner = psi_broadcast_corner + temp/(N-1);
        end
    end
    
    N_gamma_broadcast_corner = 1/(1/N+beta^2*N/2/g_inverse(gamma)^2/psi_broadcast_corner);
    
    R2C_g_num(n) = ceil(max(N_alpha, N_gamma_gossip_corner));
    R2C_b_num(n) = ceil(max(N_alpha, N_gamma_broadcast_corner));
end

%==============================
% Figure settings
figure;
title('# of validators vs Total # of nodes');
xlabel('Total # of nodes (N)');
ylabel('# of validators (N tilda)');
xlim([0 400]);
ylim([0 400]);
grid on;
hold on;

plot(N_vector, R2C_g_num, 'bo-', 'DisplayName', 'R2C with Gossip');
plot(N_vector, R2C_b_num, 'ro-', 'DisplayName', 'R2C with Broadcast');

legend('show');
hold off;
