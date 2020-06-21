classdef InN_str
    properties
        c = 5.703e-10;
        a = 3.545e-10;
        
        A1 = -8.21;
        A2 = -0.68;
        A3 = 7.57;
        A4 = -5.23;
        A5 = -5.11;
        A6 = -5.96;
        
        % Deformation eV
        D1 = -3.6 ;
        D2 = 1.7 ;
        D3 = 5.2 ;
        D4 = -2.7 ;
        
        del_cr = 0.024 ; % eV
        del_so = 0.005 ; % eV
        
        del1;
        del2;
        del3;
        
        Eg = 0.69 ;
        
        a1 = -4.2;
        a2 = -4.2;
        
        me_t = 0.07;
        me_z = 0.07;
        
        eps = 9.3;
        
        Psp = -0.032;    
        e31 = -0.57;
        e33 = 0.97;
        d31 = 3.7e-12;

        C11 = 223; % GPa
        C12 = 115;
        C13 = 92;
        C33 = 224;
    end
    
    methods
        function InN = InN_str()
            global eps0
            InN.del1 = InN.del_cr;
            InN.del2 = InN.del_so / 3;
            InN.del3 = InN.del_so / 3;
            InN.eps = InN.eps*eps0;
        end
    end
end