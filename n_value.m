function [n] = n_value(T, HTStapeName)

    if HTStapeName == 'Fujikura'
        n_meas = [7, 11, 12, 13, 15, 20, 34, 66]'; 
        T_meas = [77, 65, 50, 40, 30, 20, 10, 4.2]'; 	
        n = interp1(T_meas, n_meas, T); 
        n(T < min(T_meas)) = 66; 
        n(T > max(T_meas)) = 7;  
    end
    
end 