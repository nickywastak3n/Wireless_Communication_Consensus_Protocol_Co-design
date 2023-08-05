% Calculates single dissemination time for gossip protocol
function Z = dissemination_time_gossip(N, source, outage_probability)
    system = zeros(sqrt(N+1), sqrt(N+1));

    x_source = floor((source-1)/sqrt(N+1))+1;
    y_source = source-(x_source-1)*sqrt(N+1);
    system(x_source, y_source) = 1;

    Z = 0; % Number of loops, equals total dissemination time
    while true
        next_system = system;

        for i=1:sqrt(N+1)
            for j=1:sqrt(N+1)
                if system(i,j) == 1 % If node (i,j) has message, send message to neighboring nodes.
                    indices_of_neighbor = [i+1, j; i-1, j; i, j+1; i, j-1];
                    for index=1:length(indices_of_neighbor)
                        if indices_of_neighbor(index,1) >= 1 && indices_of_neighbor(index,1) <= sqrt(N+1) && indices_of_neighbor(index,2) >= 1 && indices_of_neighbor(index,2) <= sqrt(N+1)
                            if next_system(indices_of_neighbor(index,1),indices_of_neighbor(index,2)) ~= 1 % If neighboring node does not have message
                                next_system(indices_of_neighbor(index,1),indices_of_neighbor(index,2)) = random('Binomial',1,1-outage_probability); % Give message with probability
                            end
                        end
                    end        
                end    
            end
        end
        
        Z = Z+1; % Increase dissemination time by one time slot
        system = next_system;
        if sum(sum(system)) == N+1 % If message is received by all nodes
            break;
        end
    end
    