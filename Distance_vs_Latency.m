clc;
close all;
clear;

%==============================
% Definitions
transmitPower_g = 2.5e-3; % Transmit power for gossip communication
transmitPower_b = 100e-3; % Transmit power for broadcast communication
noisePower = 1e-13; % Noise power
ro = 10; % Target SNR

R_vector = 10:0.5:35; % Distance between two neighboring nodes
R_0 = 1;
lambda = 0.125;
PL_R_0_db = (-20)*log10(lambda/(4*pi*R_0));
eta = 3; % Path-loss exponent

N = 80; % Number of nodes in system
F = 25; % Number of faulty nodes in system

alpha = 0.99; % Target resiliency vector
zeta = 0.9999; % Target dissemination success probability

% Calculation of N_tilda for R2C analysis
a = 0.14; % Hermite-Pade approximation constant
g_inverse_positive = @(x) sqrt(-2/pi/a-log(1-x^2)/2+sqrt((2/pi/a+log(1-x^2)/2)^2-log(1-x^2)/a)); % g inverse for x>0
g_inverse = @(x) (x>=0)*g_inverse_positive(x)+(x<0)*(-1)*g_inverse_positive(x); % g inverse for all x

A = (1/3)-(F/N);
B = (F*(N-F)/((N-1)*N^2))*(g_inverse(2*alpha-1)^2);
phi = 0.5; % Correction factor for normal approximation

N_alpha = (phi*A+B*N+sqrt(2*phi*A*B*N-2*phi^2*B+B^2*N^2))/(A^2+2*B); 
N_tilda = ceil(N_alpha); % Minimum number of validator nodes for target resiliency alpha

E2E_R2C_g = zeros(1, length(R_vector));
E2E_R2C_b = zeros(1, length(R_vector));
E2E_R2C_g_exp = zeros(1, length(R_vector));
E2E_R2C_b_exp = zeros(1, length(R_vector));

proposer = randi(N+1); % Random proposer node
outage_probability = @(x1,y1,x2,y2,distance) 1-exp((-1)*10^(PL_R_0_db/10)*ro*noisePower/transmitPower_b*(sqrt((x1-x2)^2+(y1-y2)^2)*distance/R_0)^eta);
outage_probability_matrix = zeros(N+1,N+1);

validator_pool = 1:N+1;
validator_pool(proposer) = [];   
validator_set = validator_pool(randperm(N, N_tilda));

%==============================
for r=1:length(R_vector)
    % Outage probability matrix
    R = R_vector(r);

    for i=1:N+1
        for k=1:N+1
            if i==k
                continue;
            end

            x1 = floor((i-1)/sqrt(N+1))+1;
            y1 = i-(x1-1)*sqrt(N+1);
            x2 = floor((k-1) /sqrt(N+1))+1;
            y2 = k-(x2-1)*sqrt(N+1);

            outage_probability_matrix(i,k) = outage_probability(x1,y1,x2,y2,R);
        end
    end
    
    %==============================
    % R2C gossip analysis
    E2E_R2C_g(r) = (3/2*sqrt(N+1)-1)*N_tilda+sqrt(N+1)-1; % E2E latency with proposer in center
    % R2C broadcast analysis
    outage_probability_max_matrix = max(outage_probability_matrix);
    E2E_RC_b = sum(ceil(log(1-zeta^(1/N))./log(outage_probability_max_matrix)));
    E2E_R2C_b(r) = N_tilda/N*E2E_RC_b + (N-N_tilda)/N*ceil((log(1-zeta^(1/N)))/log(outage_probability_max_matrix(proposer)));
    
    %==============================
    % R2C experiment
    outage_probability_g = outage_probability_matrix(1,2); % Pick any neighboring outage probability
    z_ik = @(i,k) random('Geometric', 1-outage_probability_matrix(i,k)) + 1; % Time for dissemination from node i to k

    dissemination_time_gossip_matrix = zeros(N+1,1);
    dissemination_time_broadcast_matrix = zeros(N+1,1);

    for source=1:N+1 
        dissemination_time_gossip_matrix(source) = dissemination_time_gossip(N,source,outage_probability_g);
    
        Z_i = zeros(N,1);
        for k=1:N+1
            Z_i(k) = z_ik(source,k);
        end

        dissemination_time_broadcast_matrix(source) = max(Z_i);
    end    
      
    % R2C gossip experiment
    E2E_R2C_g_exp(r)= dissemination_time_gossip_matrix(proposer) + sum(dissemination_time_gossip_matrix(validator_set));
    % R2C broadcast experiment
    E2E_R2C_b_exp(r) = dissemination_time_broadcast_matrix(proposer) + sum(dissemination_time_broadcast_matrix(validator_set));
end

%==============================
% Figure settings
figure;
title('E2E Latency vs Distance between neighboring nodes');
xlabel('Distance between nodes (R)');
ylabel('E2E Latency (Time slots)');
grid on;
hold on;

%plot(R_vector, E2E_R2C_g, '--', 'DisplayName', 'R2C(G) Analy.');
%plot(R_vector, E2E_R2C_b, '--', 'DisplayName', 'R2C(B) Analy.');
plot(R_vector, E2E_R2C_g_exp, '--', 'DisplayName', 'R2C(G) Exp.');
plot(R_vector, E2E_R2C_b_exp, '--', 'DisplayName', 'R2C(B) Exp.');

legend('show');
hold off;


