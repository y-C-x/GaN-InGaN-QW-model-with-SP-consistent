classdef sp_pe
    properties
        Ppe_GaN = 0;
        Ppe_InGaN;
        Ppz_GaN;
        Ppz_InGaN;
        E_b;
        E_w;
        error;
    end
    
    methods
        function Pz = sp_pe(GaN,InGaN,str)
            Pz.Ppe_InGaN = 2*(InGaN.a-InGaN.a0)/InGaN.a0*(InGaN.e31-InGaN.e33*InGaN.C13/InGaN.C33);
            Pz.Ppz_GaN = Pz.Ppe_GaN + GaN.Psp;
            Pz.Ppz_InGaN = Pz.Ppe_InGaN + InGaN.Psp;
            
            Sum1 = 2*str.Lb*Pz.Ppz_GaN/GaN.eps+str.Lw*Pz.Ppz_InGaN/InGaN.eps;
            Sum2 = 2*str.Lb/GaN.eps + str.Lw/InGaN.eps;
            
            Pz.E_b = -(Sum1-Pz.Ppz_GaN*Sum2)/(GaN.eps*Sum2);
            Pz.E_w = -(Sum1-Pz.Ppz_InGaN*Sum2)/(InGaN.eps*Sum2);
            
            Pz.error = Pz.E_b * str.Lb * 2 + Pz.E_w*str.Lw;
        end
    end
end