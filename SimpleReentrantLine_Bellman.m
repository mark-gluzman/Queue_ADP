function [J, u_2D] = SimpleReentrantLine_Bellman(lambda, mu1, mu2, mu3, buffer1, buffer2, buffer3, c, factor)
    % Simple re-entrant network
    % DP method
    
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
    %Buffer 1,2,3 - buffers sizes
    
    
    
    % %buffer #1 capacity
    % b1 = 10;
    % %buffer #2 capacity
    % b2 = 10;
    % %buffer #3 capacity
    % b3 = 10;
    
    % Ax<=b
    
    
    g1 = zeros((buffer3+1)*100+(buffer2+1)*10+(buffer1+1),1);
    g3 = zeros((buffer3+1)*100+(buffer2+1)*10+(buffer1+1),1);
    
    
    P1 = zeros((buffer3+1)*100+(buffer2+1)*10+(buffer1+1), (buffer3+1)*100+(buffer2+1)*10+(buffer1+1));
    P3 = zeros((buffer3+1)*100+(buffer2+1)*10+(buffer1+1), (buffer3+1)*100+(buffer2+1)*10+(buffer1+1));
    
    for b1=1:buffer1+1
        for b2=1:buffer2+1
            for b3=1:buffer3+1
                
                
                % 1+1
                if b1~=buffer1+1
                    P1(b3*100+b2*10+b1, b3*100+b2*10+b1+1) = P1(b3*100+b2*10+b1, b3*100+b2*10+b1+1)+ lambda*(1-mu1)*(1-mu2);
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        lambda*(1-mu1)*(1-mu2)*c*[b1+1; b2; b3];
                else
                    P1(b3*100+b2*10+b1, b3*100+b2*10+b1) = P1(b3*100+b2*10+b1, b3*100+b2*10+b1)+...
                        lambda*(1-mu1)*(1-mu2);
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        lambda*(1-mu1)*(1-mu2)*c*[b1 ;b2; b3];
                end
                
                % 1+1, 2+1, 1-1, 3+1, 2-1
                if b3~=buffer3+1
                    P1(b3*100+b2*10+b1, (b3+1)*100+b2*10+b1) =P1(b3*100+b2*10+b1, (b3+1)*100+b2*10+b1)+ lambda*mu1*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        lambda*mu1*mu2*c*[b1; b2; b3+1];
                else
                    P1(b3*100+b2*10+b1, b3*100+b2*10+b1) =P1(b3*100+b2*10+b1, b3*100+b2*10+b1)+ lambda*mu1*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        lambda*mu1*mu2*c*[b1; b2; b3];
                end
                % 1+1, 2-1, 3+1
                if b1~=buffer1+1 && b2~=1 && b3~=buffer3+1
                    P1(b3*100+b2*10+b1, (b3+1)*100+(b2-1)*10+b1+1) =P1(b3*100+b2*10+b1, (b3+1)*100+(b2-1)*10+b1+1)+...
                        lambda*(1-mu1)*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        lambda*(1-mu1)*mu2*c*[b1+1; b2-1; b3+1];
                elseif b1~=buffer1+1 && b2~=1 && b3==buffer3+1
                    P1(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1+1) =P1(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1+1)+...
                        lambda*(1-mu1)*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        lambda*(1-mu1)*mu2*c*[b1+1; b2-1; b3];
                elseif b1~=buffer1+1 && b2==1 && b3==buffer3+1
                    P1(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1+1) =P1(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1+1)+...
                        lambda*(1-mu1)*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        lambda*(1-mu1)*mu2*c*[b1+1; b2; b3+1];
                elseif b1==buffer1+1 && b2==1 && b3==buffer3+1
                    P1(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1) =P1(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1)+...
                        lambda*(1-mu1)*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        lambda*(1-mu1)*mu2*c*[b1; b2; b3];
                elseif b1==buffer1+1 && b2~=1 && b3==buffer3+1
                    P1(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1) =P1(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1)+...
                        lambda*(1-mu1)*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        lambda*(1-mu1)*mu2*c*[b1; b2-1; b3];
                elseif b1==buffer1+1 && b2~=1 && b3~=buffer3+1
                    P1(b3*100+b2*10+b1, (b3+1)*100+(b2-1)*10+b1) =P1(b3*100+b2*10+b1, (b3+1)*100+(b2-1)*10+b1)+...
                        lambda*(1-mu1)*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        lambda*(1-mu1)*mu2*c*[b1; b2-1; b3+1];
                elseif b1~=buffer1+1 && b2==1 && b3~=buffer3+1
                    P1(b3*100+b2*10+b1, (b3+1)*100+(b2)*10+b1+1) =P1(b3*100+b2*10+b1, (b3+1)*100+(b2)*10+b1+1)+...
                        lambda*(1-mu1)*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        lambda*(1-mu1)*mu2*c*[b1+1; b2; b3+1];
                elseif b1==buffer1+1 && b2==1 && b3~=buffer3+1
                    P1(b3*100+b2*10+b1, (b3+1)*100+(b2)*10+b1) =P1(b3*100+b2*10+b1, (b3+1)*100+(b2)*10+b1)+...
                        lambda*(1-mu1)*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        lambda*(1-mu1)*mu2*c*[b1; b2; b3+1];
                else
                    warning('something wrong');
                end
                % 1+1, 1-1, 2+1
                if b2~=buffer2+1
                    P1(b3*100+b2*10+b1, b3*100+(b2+1)*10+b1) = P1(b3*100+b2*10+b1, b3*100+(b2+1)*10+b1)+...
                        lambda*mu1*(1-mu2);
                    
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        lambda*(mu1)*(1-mu2)*c*[b1; b2+1; b3];
                else
                    P1(b3*100+b2*10+b1, b3*100+(b2)*10+b1) = P1(b3*100+b2*10+b1, b3*100+(b2)*10+b1)+...
                        lambda*mu1*(1-mu2);
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        lambda*(mu1)*(1-mu2)*c*[b1; b2; b3];
                end
                
                %
                P1(b3*100+b2*10+b1, b3*100+b2*10+b1) =P1(b3*100+b2*10+b1, b3*100+b2*10+b1)+ (1-lambda)*(1-mu1)*(1-mu2);
                g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                    (1-lambda)*(1-mu1)*(1-mu2)*c*[b1; b2; b3];
                % 2-1, 3+1, 1-1, 2+1
                if b3~=buffer3+1 && b1~=1
                    P1(b3*100+b2*10+b1, (b3+1)*100+b2*10+b1-1) =P1(b3*100+b2*10+b1, (b3+1)*100+b2*10+b1-1)+...
                        (1-lambda)*mu1*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        (1-lambda)*(mu1)*(mu2)*c*[b1-1; b2; b3+1];
                    
                elseif b3~=buffer3+1 && b1==1
                    P1(b3*100+b2*10+b1, (b3+1)*100+b2*10+b1) =P1(b3*100+b2*10+b1, (b3+1)*100+b2*10+b1)+...
                        (1-lambda)*mu1*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        (1-lambda)*(mu1)*(mu2)*c*[b1; b2; b3+1];
                elseif b3==buffer3+1 && b1~=1
                    P1(b3*100+b2*10+b1, (b3)*100+b2*10+b1-1) =P1(b3*100+b2*10+b1, (b3)*100+b2*10+b1-1)+...
                        (1-lambda)*mu1*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        (1-lambda)*(mu1)*(mu2)*c*[b1-1; b2; b3];
                elseif b3==buffer3+1 && b1==1
                    P1(b3*100+b2*10+b1, (b3)*100+b2*10+b1) =P1(b3*100+b2*10+b1, (b3)*100+b2*10+b1)+...
                        (1-lambda)*mu1*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        (1-lambda)*(mu1)*(mu2)*c*[b1; b2; b3];
                else
                    warning('something wrong');
                end
                % 2-1, 3+1
                if b2~=1 && b3~=buffer3+1
                    P1(b3*100+b2*10+b1, (b3+1)*100+(b2-1)*10+b1) =P1(b3*100+b2*10+b1, (b3+1)*100+(b2-1)*10+b1)+...
                        (1-lambda)*(1-mu1)*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        (1-lambda)*(1-mu1)*mu2*c*[b1; b2-1; b3+1];
                elseif b2~=1 && b3==buffer3+1
                    P1(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1) =P1(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1)+...
                        (1-lambda)*(1-mu1)*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        (1-lambda)*(1-mu1)*mu2*c*[b1; b2-1; b3];
                elseif b2==1 && b3~=buffer3+1
                    P1(b3*100+b2*10+b1, (b3+1)*100+(b2)*10+b1) =P1(b3*100+b2*10+b1, (b3+1)*100+(b2)*10+b1)+...
                        (1-lambda)*(1-mu1)*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        (1-lambda)*(1-mu1)*mu2*c*[b1; b2; b3+1];
                elseif b2==1 && b3==buffer3+1
                    P1(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1) =P1(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1)+...
                        (1-lambda)*(1-mu1)*mu2;
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        (1-lambda)*(1-mu1)*mu2*c*[b1; b2; b3];
                else
                    warning('something wrong');
                end
                % 1-1, 2+1
                if b1~=1 && b2~=buffer2+1
                    P1(b3*100+b2*10+b1, b3*100+(b2+1)*10+b1-1) =P1(b3*100+b2*10+b1, b3*100+(b2+1)*10+b1-1)+...
                        (1-lambda)*mu1*(1-mu2);
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        (1-lambda)*mu1*(1-mu2)*c*[b1-1; b2+1; b3];
                elseif b1==1 && b2~=buffer2+1
                    P1(b3*100+b2*10+b1, b3*100+(b2+1)*10+b1) =P1(b3*100+b2*10+b1, b3*100+(b2+1)*10+b1)+...
                        (1-lambda)*mu1*(1-mu2);
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        (1-lambda)*mu1*(1-mu2)*c*[b1; b2+1; b3];
                elseif b1~=1 && b2==buffer2+1
                    P1(b3*100+b2*10+b1, b3*100+(b2)*10+b1-1) =P1(b3*100+b2*10+b1, b3*100+(b2)*10+b1-1)+...
                        (1-lambda)*mu1*(1-mu2);
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        (1-lambda)*mu1*(1-mu2)*c*[b1-1; b2; b3];
                elseif b1==1 && b2==buffer2+1
                    P1(b3*100+b2*10+b1, b3*100+(b2)*10+b1) =P1(b3*100+b2*10+b1, b3*100+(b2)*10+b1)+...
                        (1-lambda)*mu1*(1-mu2);
                    g1(b3*100+b2*10+b1) = g1(b3*100+b2*10+b1) +...
                        (1-lambda)*mu1*(1-mu2)*c*[b1; b2; b3];
                else
                    warning('something wrong');
                end
                
            
            end
        end
    end
    
    
    
    for b1=1:buffer1+1
        for b2=1:buffer2+1
            for b3=1:buffer3+1
                
                
                if b1~=buffer1+1
                    % 1+1
                    P3(b3*100+b2*10+b1, (b3)*100+b2*10+b1+1) = P3(b3*100+b2*10+b1, (b3)*100+b2*10+b1+1)+ lambda*(1-mu3)*(1-mu2);
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*(1-mu3)*(1-mu2)*c*[b1+1; b2; b3];
                else
                    P3(b3*100+b2*10+b1, b3*100+b2*10+b1) = P3(b3*100+b2*10+b1, b3*100+b2*10+b1)+...
                        lambda*(1-mu3)*(1-mu2);
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*(1-mu3)*(1-mu2)*c*[b1; b2; b3];
                end
                
                % 1+1, 3-1, 3+1, 2-1
                if b2~=1 && b1~=buffer1+1
                    P3(b3*100+(b2)*10+b1, (b3)*100+(b2-1)*10+b1+1) =P3(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1+1)...
                        + lambda*mu3*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*mu3*mu2*c*[b1+1; b2-1; b3];
                elseif b2~=1 && b1==buffer1+1
                    P3(b3*100+(b2)*10+b1, (b3)*100+(b2-1)*10+b1) =P3(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1)...
                        + lambda*mu3*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*mu3*mu2*c*[b1; b2-1; b3];
                elseif b2==1 && b1~=buffer1+1
                    P3(b3*100+(b2)*10+b1, (b3)*100+(b2)*10+b1+1) =P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1+1)...
                        + lambda*mu3*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*mu3*mu2*c*[b1+1; b2; b3];
                elseif b2==1 && b1==buffer1+1
                    P3(b3*100+(b2)*10+b1, (b3)*100+(b2)*10+b1) =P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1)...
                        + lambda*mu3*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*mu3*mu2*c*[b1; b2; b3];
                else
                    warning('something wrong');
                end
                
                
                % 1+1, 2-1, 3+1
                if b1~=buffer1+1 && b2~=1 && b3~=buffer3+1
                    P3(b3*100+b2*10+b1, (b3+1)*100+(b2-1)*10+b1+1) =P3(b3*100+b2*10+b1, (b3+1)*100+(b2-1)*10+b1+1)+...
                        lambda*(1-mu3)*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*(1-mu3)*mu2*c*[b1+1; b2-1; b3+1];
                elseif b1~=buffer1+1 && b2~=1 && b3==buffer3+1
                    P3(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1+1) =P3(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1+1)+...
                        lambda*(1-mu3)*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*(1-mu3)*mu2*c*[b1+1; b2-1; b3];
                elseif b1~=buffer1+1 && b2==1 && b3~=buffer3+1
                    P3(b3*100+b2*10+b1, (b3+1)*100+(b2)*10+b1+1) =P3(b3*100+b2*10+b1, (b3+1)*100+(b2)*10+b1+1)+...
                        lambda*(1-mu3)*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*(1-mu3)*mu2*c*[b1+1; b2; b3+1];
                elseif b1==buffer1+1 && b2~=1 && b3~=buffer3+1
                    P3(b3*100+b2*10+b1, (b3+1)*100+(b2-1)*10+b1) =P3(b3*100+b2*10+b1, (b3+1)*100+(b2-1)*10+b1)+...
                        lambda*(1-mu3)*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*(1-mu3)*mu2*c*[b1; b2-1; b3+1];
                elseif b1~=buffer1+1 && b2==1 && b3==buffer3+1
                    P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1+1) =P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1+1)+...
                        lambda*(1-mu3)*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*(1-mu3)*mu2*c*[b1+1; b2; b3];
                elseif b1==buffer1+1 && b2~=1 && b3==buffer3+1
                    P3(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1) =P3(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1)+...
                        lambda*(1-mu3)*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*(1-mu3)*mu2*c*[b1; b2-1; b3];
                elseif b1==buffer1+1 && b2==1 && b3~=buffer3+1
                    P3(b3*100+b2*10+b1, (b3+1)*100+(b2)*10+b1) =P3(b3*100+b2*10+b1, (b3+1)*100+(b2)*10+b1)+...
                        lambda*(1-mu3)*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*(1-mu3)*mu2*c*[b1; b2; b3+1];
                elseif b1==buffer1+1 && b2==1 && b3==buffer3+1
                    P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1) =P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1)+...
                        lambda*(1-mu3)*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*(1-mu3)*mu2*c*[b1; b2; b3];
                else
                    warning('something wrong');
                end
                
                
                % 1+1, 3-1
                if b1~=buffer1+1 && b3~=1
                    P3(b3*100+b2*10+b1, (b3-1)*100+(b2)*10+b1+1) = P3(b3*100+b2*10+b1, (b3-1)*100+(b2)*10+b1+1)+...
                        lambda*mu3*(1-mu2);
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*mu3*(1-mu2)*c*[b1+1; b2; b3-1];
                elseif b1~=buffer1+1 && b3==1
                    P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1+1) = P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1+1)+...
                        lambda*mu3*(1-mu2);
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*mu3*(1-mu2)*c*[b1+1; b2; b3];
                elseif b1==buffer1+1 && b3~=1
                    P3(b3*100+b2*10+b1, (b3-1)*100+(b2)*10+b1) = P3(b3*100+b2*10+b1, (b3-1)*100+(b2)*10+b1)+...
                        lambda*mu3*(1-mu2);
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*mu3*(1-mu2)*c*[b1; b2; b3-1];
                elseif b1==buffer1+1 && b3==1
                    P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1) = P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1)+...
                        lambda*mu3*(1-mu2);
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        lambda*mu3*(1-mu2)*c*[b1; b2; b3];
                else
                    warning('something wrong');
                end
                
                %
                P3(b3*100+b2*10+b1, b3*100+b2*10+b1) =P3(b3*100+b2*10+b1, b3*100+b2*10+b1)+ (1-lambda)*(1-mu3)*(1-mu2);
                g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                    (1-lambda)*(1-mu3)*(1-mu2)*c*[b1; b2; b3];
                % 2-1, 3+1, 3-1
                if b2~=1
                    P3(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1) =P3(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1)+...
                        (1-lambda)*mu3*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        (1-lambda)*mu3*mu2*c*[b1; b2-1; b3];
                elseif b2==1
                    P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1) =P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1)+...
                        (1-lambda)*mu3*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        (1-lambda)*mu3*mu2*c*[b1; b2; b3];
                else
                    warning('something wrong');
                end
                
                % 2-1, 3+1
                if b2~=1 && b3~=buffer3+1
                    P3(b3*100+b2*10+b1, (b3+1)*100+(b2-1)*10+b1) =P3(b3*100+b2*10+b1, (b3+1)*100+(b2-1)*10+b1)+...
                        (1-lambda)*(1-mu3)*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        (1-lambda)*(1-mu3)*mu2*c*[b1; b2-1; b3+1];
                elseif b2~=1 && b3==buffer3+1
                    P3(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1) =P3(b3*100+b2*10+b1, (b3)*100+(b2-1)*10+b1)+...
                        (1-lambda)*(1-mu3)*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        (1-lambda)*(1-mu3)*mu2*c*[b1; b2-1; b3];
                elseif b2==1 && b3~=buffer3+1
                    P3(b3*100+b2*10+b1, (b3+1)*100+(b2)*10+b1) =P3(b3*100+b2*10+b1, (b3+1)*100+(b2)*10+b1)+...
                        (1-lambda)*(1-mu3)*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        (1-lambda)*(1-mu3)*mu2*c*[b1; b2; b3+1];
                elseif b2==1 && b3==buffer3+1
                    P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1) =P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1)+...
                        (1-lambda)*(1-mu3)*mu2;
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        (1-lambda)*(1-mu3)*mu2*c*[b1; b2; b3];
                else
                    warning('something wrong');
                end
                
                % 3-1
                if b3~=1
                    P3(b3*100+b2*10+b1, (b3-1)*100+(b2)*10+b1) =P3(b3*100+b2*10+b1, (b3-1)*100+(b2)*10+b1)+...
                        (1-lambda)*mu3*(1-mu2);
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        (1-lambda)*mu3*(1-mu2)*c*[b1; b2; b3-1];
                    
                elseif b3==1
                    P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1) =P3(b3*100+b2*10+b1, (b3)*100+(b2)*10+b1)+...
                        (1-lambda)*mu3*(1-mu2);
                    g3(b3*100+b2*10+b1) = g3(b3*100+b2*10+b1) +...
                        (1-lambda)*mu3*(1-mu2)*c*[b1; b2; b3];
                else
                    warning('something wrong');
                end
            end
        end
    end
    
    
%     P1( ~any(P1,2), : ) = [];  %rows
%     P1( :, ~any(P1,1) ) = [];  %columns
%     P3( ~any(P3,2), : ) = [];  %rows
%     P3( :, ~any(P3,1) ) = [];  %columns
    
    b = [g1;g3];
%     b(~any(b)) = [];
    A = [eye(length(P1))-factor*P1;eye(length(P3))-factor*P3];
    



    lb = zeros(length(P1),1);
    f = ones(length(P1),1);
    
    J = linprog(-f,A,b,[],[],lb,[]);
%     del = find(J<0.0000001);
%     P1(del,del)=[];
%     P3(del,del)=[];
%     g1(del)=[];
%     g3(del)=[];
    u=zeros(length(g1),1);
    
     for b1=1:buffer1+1
        for b2=1:buffer2+1
            for b3=1:buffer3+1
               if g1(b3*100+b2*10+b1)+P1(b3*100+b2*10+b1,:)*J < g3(b3*100+b2*10+b1)+P3(b3*100+b2*10+b1,:)*J
                   u(b3*100+b2*10+b1)=1;
               else
                   u(b3*100+b2*10+b1)=3;
               end
            end
        end
     end
     
  
     u_3 = zeros((buffer1+1),(buffer2+1),(buffer3+1));
     u_2D = zeros((buffer1+1)*(buffer2+1)*(buffer3+1), 7);
      
     t=1;
     for b1=1:buffer1+1
        for b2=1:buffer2+1
            for b3=1:buffer3+1
              u_3(b1, b2, b3) =  u(b3*100+b2*10+b1);
              if b1>1
                  fbfs = 1;
              else
                  fbfs = 3;
              end
              if b3>1
                  lbfs = 3;
              else
                  lbfs = 1;
              end
              if b3>=b1
                  lqfs = 3;
              else
                  lqfs = 1;
              end
              
              u_2D(t, :) = [b1-1, b2-1, b3-1, u(b3*100+b2*10+b1), fbfs, lbfs, lqfs];
              t=t+1;
            end
        end
      end
     
      
      
end
