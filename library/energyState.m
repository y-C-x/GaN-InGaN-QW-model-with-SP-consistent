%% Class energyState
% this class contains all the calculated information about the energy
% states and wave functions. All energy is in the unit of eV.
classdef energyState
    properties
        Ev0; % normalized
        Ec0; % normalized
        Ev; % denormalized
        Ec; % denormalized
        psi_hh1; % wavefunction for HH1
        psi_lh1; % wavefunction for LH1
        psi_hh2; % wavefunction for HH2
        psi_lh2; % wavefunction for LH2
        psi_e1; % wavefunction for CH1
        psi_e2; % wavefunction for CH2
        CH2;
        CH1;
        HH1;
        LH1;
        HH2;
        LH2;
        Ec1_hh1;
        Ec1_lh1;
        Ec2_hh2;
        Ec2_lh2;
    end
    
    methods
        function info = energyState(m_h,m_c)
            info.Ev0 = eig(m_h.Hm);
            info.Ev0 = info.Ev0(end-3:end);
            info.Ec0 = eig(m_c.Hc);
            info.Ec0 = info.Ec0(1:2);
            
            info.Ev = info.Ev0 * m_h.Enorm;
            info.Ec = info.Ec0 * m_c.Enorm;
            
%             info.Ev = eig(m_h.Hm)*m_h.Enorm;
%             info.Ev = info.Ev(end-3:end);
%             info.Ec = eig(m_c.Hc)*m_c.Enorm;
%             info.Ec = info.Ec(1:2);
            
            info.psi_hh1 = zeros(1,length(m_c.Ct_map));
            info.psi_lh1 = zeros(1,length(m_c.Ct_map));
            info.psi_hh2 = zeros(1,length(m_c.Ct_map));
            info.psi_lh2 = zeros(1,length(m_c.Ct_map));
            info.psi_e1 = zeros(1,length(m_c.Ct_map));
            info.psi_e2 = zeros(1,length(m_c.Ct_map));
            
            info.CH1 = info.Ec(1);
            info.CH2 = info.Ec(2);
            
            info.HH1 = info.Ev(4);
            info.LH1 = info.Ev(3);
            info.HH2 = info.Ev(2);
            info.LH2 = info.Ev(1);
        end
        
        function info = getWave(info,m_h,m_c) 
            CH1_0 = info.CH1/m_c.Enorm;
            CH2_0 = info.CH2/m_c.Enorm;
            HH1_0 = info.HH1/m_h.Enorm;
            LH1_0 = info.LH1/m_h.Enorm;
            HH2_0 = info.HH2/m_h.Enorm;
            LH2_0 = info.LH2/m_h.Enorm;
            info.psi_e1 = findWaveFxE(m_c.Hc,CH1_0)';
            info.psi_e2 = findWaveFxE(m_c.Hc,CH2_0)';
            
            info.psi_hh1 = findWaveFx(m_h.Hm,HH1_0,1)';
            info.psi_lh1 = findWaveFx(m_h.Hm,LH1_0,2)';
            info.psi_hh2 = findWaveFx(m_h.Hm,HH2_0,1)';
            info.psi_lh2 = findWaveFx(m_h.Hm,LH2_0,2)';
        end 
        
        function info = getTransEnergy(info)
           info.Ec1_hh1 = info.CH1 - info.HH1;
           info.Ec1_lh1 = info.CH1 - info.LH1;
           info.Ec2_hh2 = info.CH2 - info.HH2;
           info.Ec2_lh2 = info.CH2 - info.LH2;
        end
            
            
    end
end