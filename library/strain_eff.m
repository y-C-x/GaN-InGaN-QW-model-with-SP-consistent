%% Class strain_eff
% This class includes the information about the strain effect in the
% structure.
classdef strain_eff
    properties
        % eps constant
        eps_xx;
        eps_yy;
        eps_zz;
        
        % well
        lemda_eps;
        theta_eps;
        
        Pce; % hydrostatic energy shift 
        
        % barrier
        lemda_eps_b;
        theta_eps_b;
        
        Pce_b;
    end
    methods
        function strain = strain_eff(InGaN,GaN)
            
            % eps constant
            strain.eps_xx = (InGaN.a0-InGaN.a)/InGaN.a;
            strain.eps_yy = strain.eps_xx;
            strain.eps_zz = -2*InGaN.C13/InGaN.C33*strain.eps_xx;
            
            strain.lemda_eps = InGaN.D1*strain.eps_zz+InGaN.D2*(strain.eps_xx+strain.eps_yy);
            strain.theta_eps = InGaN.D3*strain.eps_zz+InGaN.D4*(strain.eps_xx+strain.eps_yy);
            
            strain.Pce = 2*strain.eps_xx*(InGaN.a2-InGaN.a1*InGaN.C13/InGaN.C33);
            
            % barrier
            strain.lemda_eps_b = GaN.D1*strain.eps_zz+GaN.D2*(strain.eps_xx+strain.eps_yy);
            strain.theta_eps_b = GaN.D3*strain.eps_zz+GaN.D4*(strain.eps_xx+strain.eps_yy);
            
            strain.Pce_b = 2*strain.eps_xx*(GaN.a2-GaN.a1*GaN.C13/GaN.C33);
        end
    end 
end