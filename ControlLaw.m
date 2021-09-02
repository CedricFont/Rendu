function input = ControlLaw(M, i, K, h, r, param, flag)

    addpath('Config\casadi-windows-matlabR2016a-v3.5.5')
    import casadi.*

    if strcmp(flag,'I-PID')
        
        % Controller variables
        Kp = K(1);
        Ki2 = K(2);
        Kd = K(3);
        Ki1 = K(4);
        
        switching = param(1);
        t_int = param(2);

        if i < switching

            % Integrate from 'real' begining, but last one, which should
            % not appear in present input
            input = Ki1*sum(r - M(1:i))*h;

        else

            input = Kp*(r - M(i)) + ...
                Ki2*sum(r - M(i - t_int:i))*h + ...
                Kd*(M(i) - M(i-1))/h;

        end
        
    end
    
    if strcmp(flag,'PID')
        
        % Controller variables
        Kp = K(1);
        Ki2 = K(2);
        Kd = K(3);
        t_int = param(1);
        
        if i > t_int
            
            input = Kp*(r - M(i)) + ...
                Ki2*sum(r - M(i - t_int:i))*h - ...
                Kd*(M(i) - M(i-1))/h;
            
        else
            
            input = 0;
            
        end
        
    end
    
    if strcmp(flag,'PI')
        
        % Controller variables
        Kp = K(1);
        Ki2 = K(2);
        t_int = param(2);
        
        if i > t_int
            
            input = Kp*(r - M(i)) + ...
                Ki2*sum(r - M(i - t_int:i))*h;
            
        else
            
            input = 0;
            
        end
        
    end
    
    if strcmp(flag,'P')
        
       Kp = K(1);
       
       input = Kp*(r - M);
        
    end

end