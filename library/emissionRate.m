classdef emissionRate
    properties
        C0;
        C1; % [1/J-m^2-s]
        oL; % over lap
        mmE; % momentum matrix element
        rho_CH1; % density of state for CH1-HH1 transition
        rho_CL1; % density of state for CH1-LH1 transition
        rho_CH2; % density of state for CH1-HH1 transition
        rho_CL2; % density of state for CH1-LH1 transition
        
        rhoC_ch1_hh1; %[1/J-m^3]
        rhoC_ch1_lh1;
        rhoC_ch2_hh2;
        rhoC_ch2_lh2;
            
        nr; % refractive index (n1-n2)/(n1+n2)
        wavLen; % wave length
        hw; % photon energy
        w; % frequency (rad/s)
        
        % reduced mass
        mr_ch1_hh1;
        mr_ch1_lh1;
        mr_ch2_hh2;
        mr_ch2_lh2;
        
        % fermi function
        fc_c1_hh1;
        fc_c1_lh1;
        fc_c2_hh2;
        fc_c2_lh2;
        
        fv_c1_hh1;
        fv_c1_lh1;
        fv_c2_hh2;
        fv_c2_lh2;
        
        g_linewidth_CH1; % C1 - HH1;
        g_linewidth_CL1; % C1 - LH1;
        g_linewidth_CH2; % C2 - HH2;
        g_linewidth_CL2; % C2 - LH2;
        
        g_linewidth; % [1/m]
        
        g_linewidth_CH1_broad; % C1 - HH1;
        g_linewidth_CL1_broad; % C1 - LH1;
        g_linewidth_CH2_broad; % C2 - HH2;
        g_linewidth_CL2_broad; % C2 - LH2;
        
        g_linewidth_broad; % [1/m]
        
        r_sp_CH1; % emission rate
        r_sp_CL1; % emission rate
        r_sp_CH2; % emission rate
        r_sp_CL2; % emission rate
        
        r_sp; % [1/eV-s-m^3]
        
        r_sp_CH1_broad; % emission rate
        r_sp_CL1_broad; % emission rate
        r_sp_CH2_broad; % emission rate
        r_sp_CL2_broad; % emission rate
        
        r_sp_broad; % [1/eV-s-m^3]
    end
    
    
    methods
        function rate = emissionRate(overLap,MMElement,waveLength)
            global h q c0
            rate.oL = overLap;
            rate.mmE = MMElement;
            rate.wavLen = waveLength * 1e-9;
            rate.hw = (h*c0./rate.wavLen /q); % unit of Joule to eV
            rate.w = c0./rate.wavLen;
        end
        
        function rate = updateC(rate)
            global h_b q c0 eps0 m0
            rate.C0 = pi*q^2/rate.nr/c0/eps0/m0^2.*(rate.w).^-1; % SI Unit [m^2/kg]
            rate.C1 = rate.nr^2.*rate.w.^2/pi^2/h_b/c0^2; % SI unit [1/J-m^2-s]
        end
        
        function rate = getRho(rate,Fermi,ES)
            global m0 h_b
            Lw = evalin('base','Lw');
            rate.mr_ch1_hh1 = (1/Fermi.me1+1/Fermi.mhh1)^-1;
            rate.mr_ch1_lh1 = (1/Fermi.me1+1/Fermi.mlh1)^-1;
            rate.mr_ch2_hh2 = (1/Fermi.me2+1/Fermi.mhh2)^-1;
            rate.mr_ch2_lh2 = (1/Fermi.me2+1/Fermi.mlh2)^-1;
            
            rhoC = @(mr) mr*m0 / pi / h_b^2 / Lw;
            
            rate.rhoC_ch1_hh1 = rhoC(rate.mr_ch1_hh1); %[1/J-m^3]
            rate.rhoC_ch1_lh1 = rhoC(rate.mr_ch1_lh1);
            rate.rhoC_ch2_hh2 = rhoC(rate.mr_ch2_hh2);
            rate.rhoC_ch2_lh2 = rhoC(rate.mr_ch2_lh2);
            
            H_ch1_hh1 = heaviside(rate.hw-ES.Ec1_hh1);
            H_ch1_lh1 = heaviside(rate.hw-ES.Ec1_lh1);
            H_ch2_hh2 = heaviside(rate.hw-ES.Ec2_hh2);
            H_ch2_lh2 = heaviside(rate.hw-ES.Ec2_hh2);
            
            rate.rho_CH1 = rate.rhoC_ch1_hh1 * H_ch1_hh1;
            rate.rho_CL1 = rate.rhoC_ch1_lh1 * H_ch1_lh1;
            rate.rho_CH2 = rate.rhoC_ch2_hh2 * H_ch2_hh2;
            rate.rho_CL2 = rate.rhoC_ch2_lh2 * H_ch2_lh2;
        end
        
        function rate = getFermi(rate,ES,Fermi)
            global kT;
            Et_c1_hh1 = rate.hw - ES.Ec1_hh1;
            Et_c1_lh1 = rate.hw - ES.Ec1_lh1;
            Et_c2_hh2 = rate.hw - ES.Ec2_hh2;
            Et_c2_lh2 = rate.hw - ES.Ec2_lh2;
            
            Fc = Fermi.Fc;
            Fv = Fermi.Fv;
            
            fc = @(Et,mr,me,En) (1+exp(((mr/me)*Et+(En-Fc))/kT)).^-1;
            fv = @(Et,mr,mh,Ev) (1+exp(-((mr/mh)*Et+(Fv-Ev))/kT)).^-1;
            
            rate.fc_c1_hh1 = fc(Et_c1_hh1,rate.mr_ch1_hh1,Fermi.me1,ES.CH1);
            rate.fc_c1_lh1 = fc(Et_c1_lh1,rate.mr_ch1_lh1,Fermi.me1,ES.CH1);
            rate.fc_c2_hh2 = fc(Et_c2_hh2,rate.mr_ch2_hh2,Fermi.me2,ES.CH2);
            rate.fc_c2_lh2 = fc(Et_c2_lh2,rate.mr_ch2_lh2,Fermi.me2,ES.CH2);
            
            rate.fv_c1_hh1 = fv(Et_c1_hh1,rate.mr_ch1_hh1,Fermi.mhh1,ES.HH1);
            rate.fv_c1_lh1 = fv(Et_c1_lh1,rate.mr_ch1_lh1,Fermi.mlh1,ES.LH1);
            rate.fv_c2_hh2 = fv(Et_c2_hh2,rate.mr_ch2_hh2,Fermi.mhh2,ES.HH2);
            rate.fv_c2_lh2 = fv(Et_c2_lh2,rate.mr_ch2_lh2,Fermi.mlh2,ES.LH2);
            
        end
        
        function rate = getSponRate(rate)
            global q;
            rate.r_sp_CH1 = rate.oL.Ic1_hh1*rate.C0.*rate.C1.*1.5*rate.mmE.C1_HH1.*rate.rho_CH1.*rate.fc_c1_hh1.*(1-rate.fv_c1_hh1)*q*q; % [1/eV-s-m^3]
            rate.r_sp_CL1 = rate.oL.Ic1_lh1*rate.C0.*rate.C1.*0.5*rate.mmE.C1_LH1.*rate.rho_CL1.*rate.fc_c1_lh1.*(1-rate.fv_c1_lh1)*q*q;
            rate.r_sp_CH2 = rate.oL.Ic2_hh2*rate.C0.*rate.C1.*1.5*rate.mmE.C2_HH2.*rate.rho_CH2.*rate.fc_c2_hh2.*(1-rate.fv_c2_hh2)*q*q;
            rate.r_sp_CL2 = rate.oL.Ic2_lh2*rate.C0.*rate.C1.*0.5*rate.mmE.C2_LH2.*rate.rho_CL2.*rate.fc_c2_lh2.*(1-rate.fv_c2_lh2)*q*q;
            rate.r_sp = rate.r_sp_CH1+rate.r_sp_CH2+rate.r_sp_CL1+rate.r_sp_CL2;
        end
        
        
        function rate = getGain(rate)
            global q;
            rate.g_linewidth_CH1 = rate.C0*rate.oL.Ic1_hh1*1.5*rate.mmE.C1_HH1.*rate.rho_CH1.*(rate.fc_c1_hh1-rate.fv_c1_hh1)*q; %[1/m]
            rate.g_linewidth_CL1 = rate.C0*rate.oL.Ic1_lh1*0.5*rate.mmE.C1_LH1.*rate.rho_CL1.*(rate.fc_c1_lh1-rate.fv_c1_lh1)*q; %[1/m]
            rate.g_linewidth_CH2 = rate.C0*rate.oL.Ic2_hh2*1.5*rate.mmE.C2_HH2.*rate.rho_CH2.*(rate.fc_c2_hh2-rate.fv_c2_hh2)*q; %[1/m]
            rate.g_linewidth_CL2 = rate.C0*rate.oL.Ic2_lh2*0.5*rate.mmE.C2_LH2.*rate.rho_CL2.*(rate.fc_c2_lh2-rate.fv_c2_lh2)*q; %[1/m]
            rate.g_linewidth = rate.g_linewidth_CH1 + rate.g_linewidth_CL1 + rate.g_linewidth_CH2 + rate.g_linewidth_CL2;
        end
        
        function rate = getGain_broad(rate,gamma,ES)
            
            global q
            
            Et_c1_hh1 = rate.hw - ES.Ec1_hh1;
            Et_c1_lh1 = rate.hw - ES.Ec1_lh1;
            Et_c2_hh2 = rate.hw - ES.Ec2_hh2;
            Et_c2_lh2 = rate.hw - ES.Ec2_lh2;
            
            factor = @(Et,Etx) (gamma/pi)./((-Et+Etx).^2+gamma^2);
            Etx = 0 : 0.01 : 5;
            
            f_c1_hh1 = zeros(1,length(Et_c1_hh1));
            f_c1_lh1 = zeros(1,length(Et_c1_lh1));
            f_c2_hh2 = zeros(1,length(Et_c2_hh2));
            f_c2_lh2 = zeros(1,length(Et_c2_lh2));
            
            for i = 1 : length(Et_c1_hh1)
                f_c1_hh1(i) = sum(factor(Et_c1_hh1(i),Etx))/100;
                f_c1_lh1(i) = sum(factor(Et_c1_lh1(i),Etx))/100;
                f_c2_hh2(i) = sum(factor(Et_c2_hh2(i),Etx))/100;
                f_c2_lh2(i) = sum(factor(Et_c2_lh2(i),Etx))/100;
            end
            
            rate.g_linewidth_CH1_broad = rate.C0*rate.oL.Ic1_hh1*1.5*rate.mmE.C1_HH1.*rate.rhoC_ch1_hh1.*(rate.fc_c1_hh1-rate.fv_c1_hh1).*f_c1_hh1*q; %[1/m]
            rate.g_linewidth_CL1_broad = rate.C0*rate.oL.Ic1_lh1*0.5*rate.mmE.C1_LH1.*rate.rhoC_ch1_lh1.*(rate.fc_c1_lh1-rate.fv_c1_lh1).*f_c1_lh1*q; %[1/m]
            rate.g_linewidth_CH2_broad = rate.C0*rate.oL.Ic2_hh2*1.5*rate.mmE.C2_HH2.*rate.rhoC_ch2_hh2.*(rate.fc_c2_hh2-rate.fv_c2_hh2).*f_c2_hh2*q; %[1/m]
            rate.g_linewidth_CL2_broad = rate.C0*rate.oL.Ic2_lh2*0.5*rate.mmE.C2_LH2.*rate.rhoC_ch2_lh2.*(rate.fc_c2_lh2-rate.fv_c2_lh2).*f_c2_lh2*q; %[1/m]
            rate.g_linewidth_broad = rate.g_linewidth_CH1_broad + rate.g_linewidth_CL1_broad + rate.g_linewidth_CH2_broad + rate.g_linewidth_CL2_broad;

        end
        
        
        function rate = getSponRate_broad(rate,gamma,ES)
            global q kT
            
            Et_c1_hh1 = rate.hw - ES.Ec1_hh1;
            Et_c1_lh1 = rate.hw - ES.Ec1_lh1;
            Et_c2_hh2 = rate.hw - ES.Ec2_hh2;
            Et_c2_lh2 = rate.hw - ES.Ec2_lh2;
            
            factor = @(Et,Etx) (gamma/pi)./((-Et+Etx).^2+gamma^2);
            Etx = 0 : 0.01 : 5;
            
            f_c1_hh1 = zeros(1,length(Et_c1_hh1));
            f_c1_lh1 = zeros(1,length(Et_c1_lh1));
            f_c2_hh2 = zeros(1,length(Et_c2_hh2));
            f_c2_lh2 = zeros(1,length(Et_c2_lh2));
            
            for i = 1 : length(Et_c1_hh1)
                f_c1_hh1(i) = sum(factor(Et_c1_hh1(i),Etx))/100;
                f_c1_lh1(i) = sum(factor(Et_c1_lh1(i),Etx))/100;
                f_c2_hh2(i) = sum(factor(Et_c2_hh2(i),Etx))/100;
                f_c2_lh2(i) = sum(factor(Et_c2_lh2(i),Etx))/100;
            end
            
            rate.r_sp_CH1_broad = rate.oL.Ic1_hh1*rate.C0.*rate.C1.*1.5*rate.mmE.C1_HH1.*rate.rhoC_ch1_hh1.*rate.fc_c1_hh1.*(1-rate.fv_c1_hh1).*f_c1_hh1*q*q; % [1/eV-s-m^3]
            rate.r_sp_CL1_broad = rate.oL.Ic1_lh1*rate.C0.*rate.C1.*0.5*rate.mmE.C1_LH1.*rate.rhoC_ch1_lh1.*rate.fc_c1_lh1.*(1-rate.fv_c1_lh1).*f_c1_lh1*q*q;
            rate.r_sp_CH2_broad = rate.oL.Ic2_hh2*rate.C0.*rate.C1.*1.5*rate.mmE.C2_HH2.*rate.rhoC_ch2_hh2.*rate.fc_c2_hh2.*(1-rate.fv_c2_hh2).*f_c2_hh2*q*q;
            rate.r_sp_CL2_broad = rate.oL.Ic2_lh2*rate.C0.*rate.C1.*0.5*rate.mmE.C2_LH2.*rate.rhoC_ch2_lh2.*rate.fc_c2_lh2.*(1-rate.fv_c2_lh2).*f_c2_lh2*q*q;
            
            rate.r_sp_broad = rate.r_sp_CH1_broad + rate.r_sp_CL1_broad + rate.r_sp_CH2_broad + rate.r_sp_CL2_broad;
        end
        
        
    end
    
    
end