classdef FDM_C
    properties
        Cz_map;
        Ct_map;
        Vc_map;
        Enorm;
        Ltot;
        delz;
        Hc;
    end
    
    methods
        function m = FDM_C(str,GaN)
            global h2m0
            global q
            m.Cz_map = ones(1,str.index_b2)*1/GaN.me_z;
            m.Ct_map = ones(1,str.index_b2)*1/GaN.me_t;
            
            m.Enorm = h2m0/str.Ltot^2/q;
            m.Hc = zeros(str.index_b2);
            m.Ltot = str.Ltot;
            m.delz = str.delz;
        end
        
        function m = initialFDM(m,str,strain,GaN,InGaN)
            m.Cz_map(str.index_b1+1:str.index_w1) = 1/InGaN.me_z;
            m.Ct_map(str.index_b1+1:str.index_w1) = 1/InGaN.me_t;
            
            VHH_GaN = (str.Evz + GaN.del1 + GaN.del2 )/m.Enorm;
            VHH_InGaN = (str.Evz + InGaN.del1 + InGaN.del2 + strain.Pce)/m.Enorm;
            
            m.Vc_map = VHH_GaN+GaN.Eg/m.Enorm;
            m.Vc_map(str.index_b1+1:str.index_w1) = VHH_InGaN(str.index_b1+1:str.index_w1)+InGaN.Eg/m.Enorm;
        end
        
        function m = buildMatrix(m,kt)
            dz = m.delz;
            k = kt*m.Ltot*1e10;
            
            Mii = @(Mz,Mt,Mc) Mz + Mt + Mc;
            % ddz operator
            Ax_1 = @(A_1,A0) (A_1+A0)/2/dz^2;
            Ax0 = @(A_1,A0,A1) -(A_1+2*A0+A1)/2/dz^2;
            Ax1 = @(A0,A1) (A0+A1)/2/dz^2;
            
            i = 1;
            j = 1;
            while i < length(m.Hc) + 1
                while j < length(m.Hc) + 1
                    if i == j
                        index = i;
                        
                        if index == 1 || index == length(m.Hc)
                            Mz0 = -Ax0(m.Cz_map(index),m.Cz_map(index),m.Cz_map(index));
                        else
                            Mz0 = -Ax0(m.Cz_map(index-1),m.Cz_map(index),m.Cz_map(index+1));
                        end
                        
                        Mt0 = m.Ct_map(index)*k^2;
                        MC_b = m.Vc_map(index);

                        m.Hc(i,j) = Mii(Mz0,Mt0,MC_b);
                        
                        if index == 1
                            Cz1 = -Ax1(m.Cz_map(index),m.Cz_map(index+1));
                            Mz1 = Cz1;
                            m.Hc(i,j+1) = Mz1;
                        elseif index == length(m.Hc)
                            Cz_1 = -Ax_1(m.Cz_map(index-1),m.Cz_map(index));
                            Mz_1 = Cz_1;
                            m.Hc(i,j-1) = Mz_1;                         
                        else
                            Cz1 = -Ax1(m.Cz_map(index),m.Cz_map(index+1));
                            Mz1 = Cz1;
                            m.Hc(i,j+1) = Mz1;
                            
                            Cz_1 = -Ax_1(m.Cz_map(index-1),m.Cz_map(index));
                            Mz_1 = Cz_1;
                            m.Hc(i,j-1) = Mz_1; 
                        end
                    end
                    
                    j = j + 1;
                end
                i = i + 1;
                j = 1;
            end 
        end 
    end
end