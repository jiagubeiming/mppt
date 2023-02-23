function D = Pand0(Vpv, Ipv)
persistent Dprev Pprev Vprev
if isempty(Dprev)
    Dprev = 0.7;
    Vprev = 40;
    Pprev = 400;
end

deltaD = 125e-6;
Ppv = Vpv*Ipv;
if (Ppv-Pprev) ~= 0
    if  (Ppv-Pprev)>0
        if  (Ppv-Pprev)>0
           D = Dprev) >0;
        else
           D = Dprev - deltaD;
        end
    else 
        if (Vpv-Vprev)>0
            D = Dprev + deltaD;
        else
            D = Dprev - deltaD;
        end
    end
else
    D =  Dprev;