classdef MME
    properties
        C1_HH1;
        C1_LH1;
        C2_HH2;
        C2_LH2;
    end
    
    methods
        function Mb = MME(Fermi,ES,InGaN,overLap)
            Mb.C1_HH1 = Mb2(Fermi.me1,ES.CH1,ES.HH1,InGaN.del2*sqrt(2))*overLap.Ic1_hh1;
            Mb.C1_LH1 = Mb2(Fermi.me1,ES.CH1,ES.LH1,InGaN.del2*sqrt(2))*overLap.Ic1_lh1;
            Mb.C2_HH2 = Mb2(Fermi.me2,ES.CH2,ES.HH2,InGaN.del2*sqrt(2))*overLap.Ic2_hh2;
            Mb.C2_LH2 = Mb2(Fermi.me2,ES.CH2,ES.LH2,InGaN.del2*sqrt(2))*overLap.Ic2_lh2;
        end 
    end
end