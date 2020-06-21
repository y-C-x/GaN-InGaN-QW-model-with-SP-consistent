function Mb = Mb2(me,Ec,Ev,delta)
    global m0 q
    Eg = Ec-Ev; %eV
    Mb = (1/me-1)*m0*Eg*(Eg+delta)/6/(Eg+2/3*delta);
end