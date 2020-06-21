classdef InGaN_str
    properties
        c;
        a;
        a0;
        
        del_cr; % eV
        del_so; % eV
        
        del1;
        del2;
        del3;
        
        A1;
        A2;
        A3;
        A4;
        A5;
        A6;
        
        % Deformation eV
        D1;
        D2;
        D3;
        D4;
        
        Eg;
        
        a1;
        a2;
        
        me_t;
        me_z;
        
        eps;

        Psp;   
        e31;
        e33;
        d31;

        C11; % GPa
        C12;
        C13;
        C33;
    end
    
    methods
        function InGaN = InGaN_str(GaN,InN,x,Cp)
            a_const = @(x,aInN,aGaN) x*aInN+(1-x)*aGaN;
            E_bow = @(x,EInN,EGaN) x*EInN+(1-x)*EGaN-Cp*x*(1-x);
            
            InGaN.c = a_const(x,InN.c,GaN.c);
            InGaN.a = a_const(x,InN.a,GaN.a); % InGaN a_lc @ 300 K
            InGaN.a0 = GaN.a; % GaN a_lc @ 300 K
            
            InGaN.A1 = a_const(x,InN.A1,GaN.A1);
            InGaN.A2 = a_const(x,InN.A2,GaN.A2);
            InGaN.A3 = a_const(x,InN.A3,GaN.A3);
            InGaN.A4 = a_const(x,InN.A4,GaN.A4);
            InGaN.A5 = a_const(x,InN.A5,GaN.A5);
            InGaN.A6 = a_const(x,InN.A6,GaN.A6);
            
            InGaN.D1 = a_const(x,InN.D1,GaN.D1);
            InGaN.D2 = a_const(x,InN.D2,GaN.D2);
            InGaN.D3 = a_const(x,InN.D3,GaN.D3);
            InGaN.D4 = a_const(x,InN.D4,GaN.D4);
            
            InGaN.del_cr = a_const(x,InN.del_cr,GaN.del_cr);
            InGaN.del_so = a_const(x,InN.del_so,GaN.del_so);
            
            InGaN.del1 = InGaN.del_cr;
            InGaN.del2 = InGaN.del_so / 3;
            InGaN.del3 = InGaN.del_so / 3;
            
            InGaN.C13 = a_const(x,InN.C13,GaN.C13);
            InGaN.C33 = a_const(x,InN.C33,GaN.C33);
            
            InGaN.a1 = a_const(x,InN.a1,GaN.a1);
            InGaN.a2 = a_const(x,InN.a2,GaN.a2);
            
            InGaN.me_t = a_const(x,InN.me_t,GaN.me_t);
            InGaN.me_z = a_const(x,InN.me_z,GaN.me_z);
            
            InGaN.Eg = E_bow(x,InN.Eg,GaN.Eg);
            InGaN.eps = a_const(x,InN.eps,GaN.eps);
                      
            InGaN.d31 = a_const(x,InN.d31,GaN.d31);
            InGaN.C11 = a_const(x,InN.C11,GaN.C11);
            InGaN.C12 = a_const(x,InN.C12,GaN.C12);
            InGaN.C13 = a_const(x,InN.C13,GaN.C13);
            InGaN.C33 = a_const(x,InN.C33,GaN.C33);
            
            InGaN.Psp = a_const(x,InN.Psp,GaN.Psp);
            InGaN.e31 = a_const(x,InN.e31,GaN.e31);
            InGaN.e33 = a_const(x,InN.e33,GaN.e33);
        
        end

    end
end