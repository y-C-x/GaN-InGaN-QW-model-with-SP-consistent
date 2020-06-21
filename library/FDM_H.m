classdef FDM_H
    properties
        A1_map;
        A3_map;
        A13_map;
        A2_map;
        A4_map;
        A24_map;
        A5_map;
        A6_map;
        VHH_map;
        VLH_map;
        VCH_map;
        del_map;
        Enorm;
        Ltot;
        delz;
        Hm;
    end
    
    methods
        function m = FDM_H(str,GaN)
            global h2m0
            global q
            m.A1_map = ones(1,str.index_b2)*GaN.A1;
            m.A3_map = ones(1,str.index_b2)*GaN.A3;
            m.A2_map = ones(1,str.index_b2)*GaN.A2;
            m.A4_map = ones(1,str.index_b2)*GaN.A4;
            m.A5_map = ones(1,str.index_b2)*GaN.A5;
            m.A6_map = ones(1,str.index_b2)*GaN.A6;
            
            m.Enorm = h2m0/str.Ltot^2/q;
            m.Hm = zeros(3*str.index_b2);
            m.Ltot = str.Ltot;
            m.delz = str.delz;
        end
        
        function m = initialFDM(m,str,strain,GaN,InGaN)
            m.A1_map(str.index_b1+1:str.index_w1) = InGaN.A1;
            m.A2_map(str.index_b1+1:str.index_w1) = InGaN.A2;
            m.A3_map(str.index_b1+1:str.index_w1) = InGaN.A3;
            m.A4_map(str.index_b1+1:str.index_w1) = InGaN.A4;
            m.A5_map(str.index_b1+1:str.index_w1) = InGaN.A5;
            m.A6_map(str.index_b1+1:str.index_w1) = InGaN.A6;
            m.A13_map = m.A1_map + m.A3_map;
            m.A24_map = m.A2_map + m.A4_map;
            
            VHH_GaN = (str.Evz + GaN.del1 + GaN.del2 + strain.lemda_eps_b + strain.theta_eps_b)/m.Enorm;
            VLH_GaN = (str.Evz + GaN.del1 - GaN.del2 + strain.lemda_eps_b + strain.theta_eps_b)/m.Enorm;
            VCH_GaN = (str.Evz + strain.lemda_eps_b)/m.Enorm;
            
            VHH_InGaN = (str.Evz + InGaN.del1 + InGaN.del2 + strain.lemda_eps + strain.theta_eps)/m.Enorm;
            VLH_InGaN = (str.Evz + InGaN.del1 - InGaN.del2 + strain.lemda_eps + strain.theta_eps)/m.Enorm;
            VCH_InGaN = (str.Evz + strain.lemda_eps)/m.Enorm;
            
            m.VHH_map = VHH_GaN;
            m.VLH_map = VLH_GaN;
            m.VCH_map = VCH_GaN;
            m.del_map = ones(1,str.index_b2)*GaN.del3*sqrt(2)/m.Enorm;
            
            m.VHH_map(str.index_b1+1:str.index_w1) = VHH_InGaN(str.index_b1+1:str.index_w1);
            m.VLH_map(str.index_b1+1:str.index_w1) = VLH_InGaN(str.index_b1+1:str.index_w1);
            m.VCH_map(str.index_b1+1:str.index_w1) = VCH_InGaN(str.index_b1+1:str.index_w1);
            m.del_map(str.index_b1+1:str.index_w1) = InGaN.del3*sqrt(2)/m.Enorm;
        end
        
        function m = buildMatrix(m,kt)
            dz = m.delz;
            k = kt*m.Ltot*1e10;
            
            Mz = @(F,l) [F 0 0;0 F 0;0 0 l];
            Mt = @(F,K,l) [F K 0; K F 0; 0 0 l];
            Mc = @(hh,lh,ch,del) [hh 0 0;0 lh del;0 del ch]; % uintless
            Mii = @(Mz,Mt,Mc) Mz + Mt + Mc;
            Mij = @(Mz,Mt) Mz + Mt;
            % ddz operator
            Ax_1 = @(A_1,A0) (A_1+A0)/2/dz^2;
            Ax0 = @(A_1,A0,A1) -(A_1+2*A0+A1)/2/dz^2;
            Ax1 = @(A0,A1) (A0+A1)/2/dz^2;
            
            % dz operator
            Bx_1 = @(B_1,B0) -(B0+B_1)/4/dz;
            Bx1 = @(B0,B1) (B0+B1)/4/dz;
            
            i = 1;
            j = 1;
            while i < length(m.Hm)
                while j < length(m.Hm)
                    if i == j
                        index = (i+2)/3;
                        
                        if index == 1 || index == length(m.Hm)/3
                            Fz0 = -Ax0(m.A13_map(index),m.A13_map(index),m.A13_map(index));
                            lz0 = -Ax0(m.A1_map(index),m.A1_map(index),m.A1_map(index));
                        else
                            Fz0 = -Ax0(m.A13_map(index-1),m.A13_map(index),m.A13_map(index+1));
                            lz0 = -Ax0(m.A1_map(index-1),m.A1_map(index),m.A1_map(index+1));
                        end
                        
                        Mz0 = Mz(Fz0,lz0);
                        
                        Ft0  = m.A24_map(index)*k^2;
                        Kt0 = m.A5_map(index)*k^2;
                        lt0 = m.A2_map(index)*k^2;
                        
                        Mt0 = Mt(Ft0,Kt0,lt0);
                        
                        MC_b = Mc(m.VHH_map(index),m.VLH_map(index),m.VCH_map(index),m.del_map(index));
                        
                        m.Hm(i:i+2,j:j+2) = Mii(Mz0,Mt0,MC_b);
                        
                        if index == 1
                            Fz1 = -Ax1(m.A13_map(index),m.A13_map(index+1));
                            lz1 = -Ax1(m.A1_map(index),m.A1_map(index+1));
                            Mz1 = Mz(Fz1,lz1);
                            
                            H1 = k*Bx1(m.A6_map(index),m.A6_map(index+1));
                            Mt1 = [0 0 -H1;0 0 -H1;H1 H1 0];
                            
                            m.Hm(i:i+2,j+3:j+5) = Mij(Mz1,Mt1);
                        elseif index == length(m.Hm)/3
                            Fz_1 = -Ax_1(m.A13_map(index-1),m.A13_map(index));
                            lz_1 = -Ax_1(m.A1_map(index-1),m.A1_map(index));
                            Mz_1 = Mz(Fz_1,lz_1);

                            H_1 = k*Bx_1(m.A6_map(index-1),m.A6_map(index));

                            Mt_1 = [0 0 -H_1;0 0 -H_1;H_1 H_1 0];
 
                            m.Hm(i:i+2,j-3:j-1) = Mij(Mz_1,Mt_1);                           
                        else
                            Fz1 = -Ax1(m.A13_map(index),m.A13_map(index+1));
                            lz1 = -Ax1(m.A1_map(index),m.A1_map(index+1));
                            Mz1 = Mz(Fz1,lz1);
                            
                            Fz_1 = -Ax_1(m.A13_map(index-1),m.A13_map(index));
                            lz_1 = -Ax_1(m.A1_map(index-1),m.A1_map(index));
                            Mz_1 = Mz(Fz_1,lz_1);
                            
                            H1 = k*Bx1(m.A6_map(index),m.A6_map(index+1));
                            H_1 = k*Bx_1(m.A6_map(index-1),m.A6_map(index));
                            
                            Mt1 = [0 0 -H1;0 0 -H1;H1 H1 0];
                            Mt_1 = [0 0 -H_1;0 0 -H_1;H_1 H_1 0];
                            
                            m.Hm(i:i+2,j+3:j+5) = Mij(Mz1,Mt1);
                            m.Hm(i:i+2,j-3:j-1) = Mij(Mz_1,Mt_1);
                        end
                    end
                    
                    j = j+3;
                end
                i = i + 3;
                j = 1;
            end 
        end 
    end
end