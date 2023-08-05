clc;
close all;
clear;

%==============================
% Definitions
transmitPower_g = 2.5e-3; % Transmit power for gossip communication
transmitPower_b = 100e-3; % Transmit power for broadcast communication
noisePower = 1e-13; % Noise power
ro = 10; % Target SNR

R = 10; % Distance between two neighboring nodes
R_0 = 1;
lambda = 0.125;
PL_R_0_db = (-20)*log10(lambda/(4*pi*R_0));
eta = 3; % Path-loss exponent
  
N = 80; % Number of nodes in system
F = 1:25; % Number of faulty nodes in system

alpha = 0.99;
beta = 1; % Acceptable consensus distortion
gamma = 0.9;

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

zeta = 0.9999; % Target dissemination success probability

%==============================
% Calculation of N_tilda for R2C analysis
a = 0.14; % Hermite-Pade approximation constant
g_inverse_positive = @(x) sqrt(-2/pi/a-log(1-x.^2)/2+sqrt((2/pi/a+log(1-x.^2)./2).^2-log(1-x.^2)./a)); % g inverse for x>0
g_inverse = @(x) (x>=0).*g_inverse_positive(x)+(x<0)*(-1).*g_inverse_positive(x); % g inverse for all x

A = (1/3)-(F/N);
B = (F.*(N-F)/((N-1)*N^2))*(g_inverse(2*alpha-1).^2);
phi = 0.5; % Correction factor for normal approximation

N_alpha = (phi*A+B*N+sqrt(2*phi*A.*B*N-2*phi^2*B+B.^2*N^2))./(A.^2+2*B); 

psi_gossip_center = (13*N^2-4*N-8)*N/24/(N-1);
psi_gossip_corner = (N+1)*((13*N-24*sqrt(N+1)+16)*N+12*(sqrt(N+1)-1))/6/(N-1);

N_gamma_gossip_center = 1/(1/N+beta^2*N/2/g_inverse(gamma)^2/psi_gossip_center);
N_tilda_gossip_center = ceil(max(N_alpha,N_gamma_gossip_center));
N_gamma_gossip_corner = 1/(1/N+beta^2*N/2/g_inverse(gamma)^2/psi_gossip_corner);
N_tilda_gossip_corner = ceil(max(N_alpha,N_gamma_gossip_corner));

psi_broadcast_center = 0;
psi_broadcast_corner = 0;

for v=1:N+1
    if v ~= 40
        psi_broadcast_center = psi_broadcast_center + (1+outage_probability_matrix(40,v))/(1-outage_probability_matrix(40,v))^2;
        
        temp = 0;
        for j=1:N+1
            if j ~= 40 && j~= v
                temp = temp + 1/(1-outage_probability_matrix(40,v))/(1-outage_probability_matrix(40,j));
            end
        end
        psi_broadcast_center = psi_broadcast_center + temp/(N-1);
    end
end

for v=1:N+1
    if v ~= 1
        psi_broadcast_corner = psi_broadcast_corner + (1+outage_probability_matrix(1,v))/(1-outage_probability_matrix(1,v))^2;
        
        temp = 0;
        for j=1:N+1
            if j ~= 1 && j~= v
                temp = temp + 1/(1-outage_probability_matrix(1,v))/(1-outage_probability_matrix(1,j));
            end
        end
        psi_broadcast_corner = psi_broadcast_corner + temp/(N-1);
    end
end

N_gamma_broadcast_center = 1/(1/N+beta^2*N/2/g_inverse(gamma)^2/psi_broadcast_center);
N_tilda_broadcast_center = ceil(max(N_alpha,N_gamma_broadcast_center));
N_gamma_broadcast_corner = 1/(1/N+beta^2*N/2/g_inverse(gamma)^2/psi_broadcast_corner);
N_tilda_broadcast_corner = ceil(max(N_alpha,N_gamma_broadcast_corner));

%==============================
% RC gossip analysis
E2E_RC_g = ((3*sqrt(N+1)-2)*(N+1)-sqrt(N+1))/2;

%==============================
% RC broadcast analysis
outage_probability_max_matrix = max(outage_probability_matrix);
E2E_RC_b = sum(ceil(log(1-zeta^(1/N))./log(outage_probability_max_matrix)));

%==============================
% R2C gossip analysis
E2E_R2C_g_cen = (3/2*sqrt(N+1)-1)*N_tilda_gossip_center+sqrt(N+1)-1; % E2E latency with proposer in center
E2E_R2C_g_cor = (3/2*sqrt(N+1)-(sqrt(N+1)-1)/N-1)*N_tilda_gossip_corner+2*(sqrt(N+1)-1); % E2E latency with proposer in corner

%==============================
% R2C broadcast analysis
E2E_R2C_b_cen = N_tilda_broadcast_center/N.*E2E_RC_b + (N-N_tilda_broadcast_center)/N*ceil((log(1-zeta^(1/N)))/log(outage_probability_max_matrix(40)));
E2E_R2C_b_cor = N_tilda_broadcast_corner/N.*E2E_RC_b + (N-N_tilda_broadcast_corner)/N*ceil((log(1-zeta^(1/N)))/log(outage_probability_max_matrix(1)));

%======================================================================================================================
% RC experiment
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

% RC gossip experiment
E2E_RC_g_exp = sum(dissemination_time_gossip_matrix);

% RC broadcast experiment
E2E_RC_b_exp = sum(dissemination_time_broadcast_matrix);

%==============================
%proposer = randi(N+1); % Random proposer node

% R2C experiment
proposer = 41;
E2E_R2C_g_cen_exp = zeros(1,length(N_tilda_gossip_center));
E2E_R2C_b_cen_exp = zeros(1,length(N_tilda_broadcast_center));

for index=1:length(N_tilda_gossip_center)
    n_tilda = N_tilda_gossip_center(index);
    
    validator_pool = 1:N+1;
    validator_pool(proposer) = [];   
    validator_set = validator_pool(randperm(N, n_tilda));
    
    % R2C gossip experiment
    E2E_R2C_g_cen_exp(index) = dissemination_time_gossip_matrix(proposer) + sum(dissemination_time_gossip_matrix(validator_set));
    

    n_tilda = N_tilda_broadcast_center(index);
    
    validator_pool = 1:N+1;
    validator_pool(proposer) = [];   
    validator_set = validator_pool(randperm(N, n_tilda));

    % R2C broadcast experiment
    E2E_R2C_b_cen_exp(index) = dissemination_time_broadcast_matrix(proposer) + sum(dissemination_time_broadcast_matrix(validator_set));
end

proposer = 1;
E2E_R2C_g_cor_exp = zeros(1,length(N_tilda_gossip_corner));
E2E_R2C_b_cor_exp = zeros(1,length(N_tilda_broadcast_corner));

for index=1:length(N_tilda_gossip_corner)
    n_tilda = N_tilda_gossip_corner(index);
    
    validator_pool = 1:N+1;
    validator_pool(proposer) = [];   
    validator_set = validator_pool(randperm(N, n_tilda));
    
    % R2C gossip experiment
    E2E_R2C_g_cor_exp(index) = dissemination_time_gossip_matrix(proposer) + sum(dissemination_time_gossip_matrix(validator_set));
    
    n_tilda = N_tilda_broadcast_corner(index);
    
    validator_pool = 1:N+1;
    validator_pool(proposer) = [];   
    validator_set = validator_pool(randperm(N, n_tilda));

    % R2C broadcast experiment
    E2E_R2C_b_cor_exp(index) = dissemination_time_broadcast_matrix(proposer) + sum(dissemination_time_broadcast_matrix(validator_set));
end

%==============================
% Figure settings
figure;
title('E2E Latency vs # of Faulty Nodes');
xlabel('# of Faulty Nodes (F)');
ylabel('E2E Latency (time slots)');
xlim([0 25]);
ylim([0 1200]);
grid on;
hold on;

plot(F, E2E_RC_g*ones(1,size(F,2)), 'bo', 'DisplayName', 'RC(G) Analy.');
plot(F, E2E_R2C_g_cen,'b--', 'DisplayName', 'R2C(G) Analy.');
plot(F, E2E_R2C_g_cor,'b--', 'HandleVisibility','off');
plot(F, E2E_RC_b*ones(1,size(F,2)), 'ro', 'DisplayName', 'RC(B) Analy.');
plot(F, E2E_R2C_b_cen,'r--', 'DisplayName', 'R2C(B) Analy.');
plot(F, E2E_R2C_b_cor,'r--', 'HandleVisibility','off');

plot(F, E2E_RC_g_exp*ones(1,size(F,2)),"co", 'DisplayName', 'RC(G) Exp.');
plot(F, E2E_R2C_g_cen_exp,"c--", 'DisplayName', 'R2C(G) Exp.');
plot(F, E2E_R2C_g_cor_exp,"c--", 'HandleVisibility','off');
plot(F, E2E_RC_b_exp*ones(1,size(F,2)),"mo", 'DisplayName', 'RC(B) Exp.');
plot(F, E2E_R2C_b_cen_exp,"m--", 'DisplayName', 'R2C(B) Exp.');
plot(F, E2E_R2C_b_cor_exp,"m--", 'HandleVisibility','off');

legend('show');
hold off;
