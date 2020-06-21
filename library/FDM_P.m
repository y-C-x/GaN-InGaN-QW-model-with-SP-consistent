classdef FDM_P
    properties
        eps_map;
        Ltot;
        dz;
        Hp;
        Vsc;
    end
    
    methods
        function m = FDM_P(str,GaN)
            m.eps_map = ones(1,str.index_b2)*GaN.eps;
            m.Ltot = str.Ltot;
            m.dz = str.dz;
            m.Hp = zeros(str.index_b2);
        end
        
        function m = initialFDM(m,str,InGaN)
            m.eps_map(str.index_b1+1:str.index_w1) = InGaN.eps;
        end
        
        function m = buildMatrix(m)
            dz = m.dz;

            % ddz operator
            Ax_1 = @(A_1,A0) (A_1+A0)/2/dz^2;
            Ax0 = @(A_1,A0,A1) -(A_1+2*A0+A1)/2/dz^2;
            Ax1 = @(A0,A1) (A0+A1)/2/dz^2;
            
            i = 1;
            j = 1;
            while i < length(m.Hp) + 1
                while j < length(m.Hp) + 1
                    if i == j
                        index = i;
                        
                        if index == 1 || index == length(m.Hp)
                            Mz0 = Ax0(m.eps_map(index),m.eps_map(index),m.eps_map(index));
                        else
                            Mz0 = Ax0(m.eps_map(index-1),m.eps_map(index),m.eps_map(index+1));
                        end
                        
                        m.Hp(i,j) = Mz0;
                        
                        if index == 1
                            Mz1 = Ax1(m.eps_map(index),m.eps_map(index+1));
                            m.Hp(i,j+1) = Mz1;
                        elseif index == length(m.Hp)
                            Mz_1 = Ax_1(m.eps_map(index-1),m.eps_map(index));
                            m.Hp(i,j-1) = Mz_1;
                        else
                            Mz1 = Ax1(m.eps_map(index),m.eps_map(index+1));
                            m.Hp(i,j+1) = Mz1;
                            
                            Mz_1 = Ax_1(m.eps_map(index-1),m.eps_map(index));
                            m.Hp(i,j-1) = Mz_1;
                        end
                    end
                    
                    j = j + 1;
                end
                i = i + 1;
                j = 1;
            end
        end
        
        function m = solveVsc(m,rho)
            m.Vsc  = (m.Hp\-rho')';
        end
    end
    
end