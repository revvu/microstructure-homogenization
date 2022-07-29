Cstar11=h.*l.*(C11.*A11+C12.*A21+C13.*A31);
Cstar12=h.*l.*(C11.*A12+C12.*A22+C13.*A32);
Cstar13=h.*l.*(C11.*A13+C12.*A23+C13.*A33);
Cstar14=h.*l.*(C11.*A14+C12.*A24+C13.*A34);
Cstar15=h.*l.*(C11.*A15+C12.*A25+C13.*A35);
Cstar16=h.*l.*(C11.*A16+C12.*A26+C13.*A36);

Cstar21=h.*l.*(C12.*A11+C22.*A21+C23.*A31);
Cstar22=h.*l.*(C12.*A12+C22.*A22+C23.*A32);
Cstar23=h.*l.*(C12.*A13+C22.*A23+C23.*A33);
Cstar24=h.*l.*(C12.*A14+C22.*A24+C23.*A34);
Cstar25=h.*l.*(C12.*A15+C22.*A25+C23.*A35);
Cstar26=h.*l.*(C12.*A16+C22.*A26+C23.*A36);

Cstar31=h.*l.*(C13.*A11+C23.*A21+C33.*A31);
Cstar32=h.*l.*(C13.*A12+C23.*A22+C33.*A32);
Cstar33=h.*l.*(C13.*A13+C23.*A23+C33.*A33);
Cstar34=h.*l.*(C13.*A14+C23.*A24+C33.*A34);
Cstar35=h.*l.*(C13.*A15+C23.*A25+C33.*A35);
Cstar36=h.*l.*(C13.*A16+C23.*A26+C33.*A36);

Cstar41=h.*l.*C44.*A41;
Cstar42=h.*l.*C44.*A42;
Cstar43=h.*l.*C44.*A43;
Cstar44=h.*l.*C44.*A44;
Cstar45=h.*l.*C44.*A45;
Cstar46=h.*l.*C44.*A46;

Cstar51=h.*l.*C55.*A51;
Cstar52=h.*l.*C55.*A52;
Cstar53=h.*l.*C55.*A53;
Cstar54=h.*l.*C55.*A54;
Cstar55=h.*l.*C55.*A55;
Cstar56=h.*l.*C55.*A56;

Cstar61=h.*l.*C66.*A61;
Cstar62=h.*l.*C66.*A62;
Cstar63=h.*l.*C66.*A63;
Cstar64=h.*l.*C66.*A64;
Cstar65=h.*l.*C66.*A65;
Cstar66=h.*l.*C66.*A66;


Cstar=[sum(Cstar11) sum(Cstar12) sum(Cstar13) sum(Cstar14) sum(Cstar15) sum(Cstar16);....
sum(Cstar21) sum(Cstar22) sum(Cstar23) sum(Cstar24) sum(Cstar25) sum(Cstar26);....    
sum(Cstar31) sum(Cstar32) sum(Cstar33) sum(Cstar34) sum(Cstar35) sum(Cstar36);....    
sum(Cstar41) sum(Cstar42) sum(Cstar43) sum(Cstar44) sum(Cstar45) sum(Cstar46);....    
sum(Cstar51) sum(Cstar52) sum(Cstar53) sum(Cstar54) sum(Cstar55) sum(Cstar56);....    
sum(Cstar61) sum(Cstar62) sum(Cstar63) sum(Cstar64) sum(Cstar65) sum(Cstar66)];    


Gamma11=C11.*alpha11+C12.*alpha22+C13.*alpha33;
Gamma21=C12.*alpha11+C22.*alpha22+C23.*alpha33;
Gamma31=C13.*alpha11+C23.*alpha22+C33.*alpha33;

Gammastar11=h.*l.*(A11.*Gamma11+A21.*Gamma21+A31.*Gamma31);
Gammastar21=h.*l.*(A12.*Gamma11+A22.*Gamma21+A32.*Gamma31);
Gammastar31=h.*l.*(A13.*Gamma11+A23.*Gamma21+A33.*Gamma31);
Gammastar41=h.*l.*(A14.*Gamma11+A24.*Gamma21+A34.*Gamma31);
Gammastar51=h.*l.*(A15.*Gamma11+A25.*Gamma21+A35.*Gamma31);
Gammastar61=h.*l.*(A16.*Gamma11+A26.*Gamma21+A36.*Gamma31);

Gammastar=[sum(Gammastar11);sum(Gammastar21);sum(Gammastar31);...
    sum(Gammastar41);sum(Gammastar51);sum(Gammastar61)];
