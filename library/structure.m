%% Class structure
% this class define all relavent parameters in a QW structure 
classdef structure
    properties
        dz = 1e-10; % the step scale
        delz; % normalized step
        Ltot; % total length of the structure
        Lb; % length of the barrier
        Lw; % length of the well
        Ev0; % valence band energy level (not used in the calculation)
        Evz0; % valence band energy level with Pz
        Evz; % valence band energy level for calculation
        Ec0; % conduction band energy level
        Ecz0; % conduction band energy level with Pz
        Ecz; % conduction band energy level for calculation
        Eoff_v; % valence band energy offset
        Eoff_c; % conduction band energy offset
        
        % this three index can change due to extra QW layers
        % if a more complex layer is added or defined, one should update
        % the parameter according in the following functions as well as 
        % the FDM matrices and all other classes except the material
        % parameters class
        
        index_b1; % index of the first barrier end in the structure
        index_b2; % index of the second barrier end in the structure
        index_w1; % index of the first well end in the structure
    end
    
    methods
        function str = structure(Lb,Lw)
            str.Lb = Lb;
            str.Lw = Lw;
            str.Ltot = 2*Lb + Lw;
            str.index_b1 = str.Lb*1e10;
            str.index_w1 = (str.Lb+str.Lw)*1e10;
            str.index_b2 = str.Ltot*1e10;
            str.delz = str.dz/str.Ltot;
        end
        
        function str = setEoff(str,VBO,GaN,InGaN) % add valence band offset energy
            str.Eoff_v = VBO * abs(GaN.Eg-InGaN.Eg);
            str.Eoff_c = (1-VBO)*abs(GaN.Eg-InGaN.Eg);
        end
        
        function str = initialEcv(str,GaN,InGaN)
            str.Ev0 = zeros(1,str.Ltot/str.dz);
            str.Ec0 = str.Ev0 + GaN.Eg;
            str.Ev0(str.index_b1+1:str.index_w1) = str.Ev0(str.index_b1+1:str.index_w1)+str.Eoff_v;
            str.Ec0(str.index_b1+1:str.index_w1) = str.Ev0(str.index_b1+1:str.index_w1)+InGaN.Eg;
            
            str.Evz = str.Ev0;
            str.Ecz = str.Ec0;
        end
        
        function str = addStrain(str,GaN,InGaN,strain)
            str.Ecz = str.Evz + GaN.Eg + GaN.del1 + GaN.del2 - strain.theta_eps_b - strain.lemda_eps_b;
            Eg_shift = InGaN.Eg + InGaN.del1 + InGaN.del2 - strain.theta_eps - strain.lemda_eps + strain.Pce;
            str.Ecz(str.index_b1+1:str.index_w1) = str.Evz(str.index_b1+1:str.index_w1)+Eg_shift;
        end
        
        function str = addPz(str,Pz)
            for i = 2 : str.index_b1
                str.Evz(i) = str.Evz(i-1) + str.dz*Pz.E_b;
                str.Ecz(i) = str.Ecz(i-1) + str.dz*Pz.E_b;
            end
            
            str.Evz(str.index_b1+1) = str.Evz(str.index_b1)+str.Eoff_v;
            str.Ecz(str.index_b1+1) = str.Ecz(str.index_b1)-str.Eoff_c;
            
            for i = str.index_b1 + 2 : str.index_w1
                str.Evz(i) = str.Evz(i-1) + str.dz*Pz.E_w;
                str.Ecz(i) = str.Ecz(i-1) + str.dz*Pz.E_w;
            end
            
            str.Evz(str.index_w1+1) = str.Evz(str.index_w1)-str.Eoff_v;
            str.Ecz(str.index_w1+1) = str.Ecz(str.index_w1)+str.Eoff_c;
            
            for i = str.index_w1 + 2 : str.index_b2
                str.Evz(i) = str.Evz(i-1) + str.dz*Pz.E_b;
                str.Ecz(i) = str.Ecz(i-1) + str.dz*Pz.E_b;
            end
        end
        
        function str = saveOrig(str)
            str.Ecz0 = str.Ecz;
            str.Evz0 = str.Evz;
        end
        
        function str = updateStr(str,Vsc)
            str.Evz = str.Evz - Vsc;
            str.Ecz = str.Ecz - Vsc;
        end
        
        function plotStr0(str)
            figure
            hold on
            plot(str.Ev0)
            plot(str.Ec0)
        end
        
        function plotStr1(str)
            figure
            hold on
            plot(str.Evz)
            plot(str.Ecz)
        end
        
        function plotStrZ0(str)
            figure
            hold on
            plot(str.Evz0)
            plot(str.Ecz0)
        end
        
    end
end