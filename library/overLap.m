classdef overLap
    properties
        % CH1 - HH1
        Ic1_hh1;
        % CH1 - LH1
        Ic1_lh1;
        % CH2 - HH2
        Ic2_hh2;
        % CH2 - LH2
        Ic2_lh2;
    end
    
    methods
        function oL = overLap(ES)
            % CH1 - HH1
            oL.Ic1_hh1 = Inm2(ES.psi_e1,ES.psi_hh1);
            % CH1 - LH1
            oL.Ic1_lh1 = Inm2(ES.psi_e1,ES.psi_lh1);
            % CH1 - HH2
            oL.Ic2_hh2 = Inm2(ES.psi_e2,ES.psi_hh2);
            % CH1 - LH2
            oL.Ic2_lh2 = Inm2(ES.psi_e2,ES.psi_lh2);
        end
    end
end