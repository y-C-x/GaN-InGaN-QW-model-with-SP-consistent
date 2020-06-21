classdef fermi_level
   properties
       mhh1;
       mlh1;
       mhh2;
       mlh2;
       
       me1;
       me2;
       Fc;
       Fv;
       
       n;
       p;
       
       rho;
   end
   
   methods
       function Fermi = fermi_level(InGaN,strain,n)
           Fermi.n = n;
           Fermi.p = n;
           
           Fermi.mhh1 = 1/(-(InGaN.A2+InGaN.A4));
           Fermi.mhh2 = Fermi.mhh1;
           
           del_const = (InGaN.del1-InGaN.del2+strain.theta_eps)/2;
           E0_2 = del_const+strain.lemda_eps+sqrt(del_const^2+2*InGaN.del3^2);
           E0_3 = del_const+strain.lemda_eps-sqrt(del_const^2+2*InGaN.del3^2);
           
           Fermi.mlh1 = -1/(InGaN.A2+(E0_2-strain.lemda_eps)/(E0_2-E0_3)*InGaN.A4);
           Fermi.mlh2 = Fermi.mlh1;
           
           Fermi.me1 = InGaN.me_t;
           Fermi.me2 = Fermi.me1;
       end
       
       function Fermi = solveFermi(Fermi,str,ES)
           global m0 kT q h_b
           syms Fcx Fvx;
           Nvc = @(m) m*m0*kT*q/pi/h_b^2/str.Lw; % all SI unit [1/m^3]
           
%            mec = (1/Fermi.me1+1/Fermi.me2)^-1;
%            nc1 = Nvc(mec);
%            nc2 = Nvc(mec);
           
%            mhc = (1/Fermi.mhh1+1/Fermi.mhh2+1/Fermi.mlh1+1/Fermi.mlh2)^-1;
%            nhh1 = Nvc(mhc);
%            nlh1 = Nvc(mhc);
%            nhh2 = Nvc(mhc);
%            nlh2 = Nvc(mhc);           
           
           nc1 = Nvc(Fermi.me1);
           nc2 = Nvc(Fermi.me2);
           nhh1 = Nvc(Fermi.mhh1);
           nlh1 = Nvc(Fermi.mlh1);
           nhh2 = Nvc(Fermi.mhh2);
           nlh2 = Nvc(Fermi.mlh2);
           
           eqnc1 = nc1*log(1+exp((Fcx-ES.CH1)/kT));
           eqnc2 = nc2*log(1+exp((Fcx-ES.CH2)/kT));
           eqnN = eqnc1+eqnc2 == Fermi.n;
           
           Fc = double(solve(eqnN,Fcx));
           
           eqnhh1 = nhh1*log(1+exp((ES.HH1-Fvx)/kT)); % m-3
           eqnlh1 = nlh1*log(1+exp((ES.LH1-Fvx)/kT));
           eqnhh2 = nhh2*log(1+exp((ES.HH2-Fvx)/kT));
           eqnlh2 = nlh2*log(1+exp((ES.LH2-Fvx)/kT));
           eqnP = eqnhh1+eqnhh2+eqnlh1+eqnlh2 == Fermi.p;
           
           Fv = double(solve(eqnP,Fvx));
           
           ns_e1 = nc1*log(1+exp((Fc-ES.CH1)/kT))*ES.psi_e1.^2; % m-3
           ns_e2 = nc2*log(1+exp((Fc-ES.CH2)/kT))*ES.psi_e2.^2;
           
           ps_hh1 = nhh1*log(1+exp((ES.HH1-Fv)/kT))*ES.psi_hh1.^2;
           ps_lh1 = nlh1*log(1+exp((ES.LH1-Fv)/kT))*ES.psi_lh1.^2;
           ps_hh2 = nhh2*log(1+exp((ES.HH2-Fv)/kT))*ES.psi_hh2.^2;
           ps_lh2 = nlh2*log(1+exp((ES.LH2-Fv)/kT))*ES.psi_lh2.^2;
           
           ns_e = ns_e1+ns_e2;
           ps_h = ps_hh1+ps_lh1+ps_hh2+ps_lh2;
           
           % this equation is given in both Hongping's paper and Chuang's book in ch.4 I
           % guess.
           rho = q*(ps_h-ns_e); % C/m^3
           
           Fermi.Fc = Fc;
           Fermi.Fv = Fv;
           Fermi.rho = rho;
       end
       
   end
end