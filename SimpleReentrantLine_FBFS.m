function [size, ci] = SimpleReentrantLine_FBFS(lambda, mu1, mu2, mu3, N, I_N)
    % Simple re-entrant network
    % FBFS
    
    %            ^
    %            |
    % -> b1] ->(c1) -> b2] -> (c2)-
    %            ^                 |
    %            b3                |
    %            ^                 |
    %            |<-----------------
    
    
    %lambda - arriving rate (Poisson process mean)
    %mu1 - service time rate for server #1 class 1
    %mu2 - service time rate for server #2
    %mu3 - service time rate for server #1 class 3
    %N - number of time periods
    %I_N - number of simulations
    %ci - confidential intervals
    %size - sum of total queues sizes
    
    
    % %buffer #1 capacity
    % b1 = 10;
    % %buffer #2 capacity
    % b2 = 10;
    % %buffer #3 capacity
    % b3 = 10;
    
    % initial buffers sizes
    z0 = ones(3,1);
    
    
    
    
    
    
    uniform_rate = lambda + mu1 + mu2 + mu3;
    
    p_arriving = lambda / uniform_rate;
    p_compl1 = mu1 / uniform_rate;
    p_compl2 = mu2 / uniform_rate;
    p_compl3 = mu3 / uniform_rate;
    
    % we will store all history; z is queue lengths
    
    iter = zeros(I_N,1);
    
    
    
    parfor i=1:I_N
        % FBFS
        z = zeros(3, N+1);
        z(:, 1) = z0;
        for t = 2:N+1
            w = rand(1);
            if w < p_arriving
                %    if z(1, t - 1) + 1<=b1
                z(1, t) = z(1, t - 1) + 1;
                %       else
                %             z(1, t) = z(1, t - 1);
                %         end
                z(2, t) = z(2, t - 1);
                z(3, t) = z(3, t - 1);
                
            elseif w < (p_compl1 + p_arriving)
                
                if  z(1, t-1)>0
                    z(1, t) = z(1, t-1) - 1;
                    %        if z(2, t - 1) + 1<=b2
                    z(2, t) = z(2, t - 1) + 1;
                    %         else
                    %             z(2, t) = z(2, t - 1);
                    %          end
                    z(3, t) = z(3, t-1);
                else
                    z(1, t) = z(1, t-1);
                    z(2, t) = z(2, t-1);
                    z(3, t) = z(3, t-1);
                end
                
            elseif    w < (p_compl1 + p_compl2 + p_arriving)
                if z(2, t-1)>0
                    z(2, t) = z(2, t-1) - 1;
                    %         if z(3, t - 1) + 1<=b3
                    z(3, t) = z(3, t - 1) + 1;
                    %          else
                    %               z(3, t) = z(3, t - 1);
                    %          end
                    z(1, t) = z(1, t-1);
                else
                    z(3, t) = z(3, t-1);
                    z(2, t) = z(2, t-1);
                    z(1, t) = z(1, t-1);
                end
                
            elseif z(1, t-1)==0 && z(3, t-1)>0
                z(1, t) = z(1, t-1);
                z(2, t) = z(2, t-1);
                z(3, t) = z(3, t-1) - 1;
            else
                z(1, t) = z(1, t-1);
                z(2, t) = z(2, t-1);
                z(3, t) = z(3, t-1);
            end
            
            
            
            
        end
        
        iter(i) = sum(z(1, :)+ z(2, :)+z(3, :))/t;
    end
    
    
    ci = [mean(iter) - 1.96*var(iter)/sqrt(I_N), mean(iter) + 1.96*var(iter)/sqrt(I_N) ];
    
    size = mean(iter);
    
    
    
    
    
  
end
