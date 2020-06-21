classdef GaN_str
    properties
        c = 5.1850e-10;
        a = 3.189e-10;
        
        del_cr = 0.010; % eV
        del_so = 0.017; % eV
        
        del1;
        del2;
        del3;
        
        A1 = -7.21;
        A2 = -0.44;
        A3 = 6.68;
        A4 = -3.46;
        A5 = -3.4;
        A6 = -4.9;
        
        % Deformation eV
        D1 = -3.6 ;
        D2 = 1.7 ;
        D3 = 5.2 ;
        D4 = -2.7 ;
        
        Eg = 3.437;
        
        a1 = -7.1;
        a2 = -9.9;
        
        me_t = 0.20;
        me_z = 0.20;
        
        eps = 9.7;
        
        Psp = -0.029;
        e31 = -0.35;
        e33 = 1.27;
        d31 = 1.05e-12;
        
        C11 = 390; % GPa
        C12 = 145;
        C13 = 106;
        C33 = 398;
    end
    
    methods
        function GaN = GaN_str()
            global eps0
            GaN.del1 = GaN.del_cr;
            GaN.del2 = GaN.del_so / 3;
            GaN.del3 = GaN.del_so / 3;
            GaN.eps = GaN.eps*eps0;
        end
    end
end