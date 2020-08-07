function Sigma=f_getSigmaT(Rs,lon1,lat1)
global lonGrid;
global latGrid;

Re=6.371e6;			%Earth radius

cosbeta=cos(lat1)*cos(lon1)*cos(lonGrid).*cos(latGrid)+cos(lat1)*sin(lon1)*sin(lonGrid).*cos(latGrid)+sin(lat1)*sin(latGrid);
				%beta refers to the angel between the spatial diretion S and U

A_sigma_RU=Rs*cosbeta-Re;		%Abnormalized cos angle between Earth Unit and Radius
A_sigma_RU(find(A_sigma_RU<0))=0;	%Check if blocked by the ground
					
A_sigma_RS=Rs-Re*cosbeta;		%Abnormalized cos angle between Spatial Direction and Radius

d2=Rs*Rs+Re*Re-2*Rs*Re*cosbeta;

Sigma=A_sigma_RU.*A_sigma_RS./d2./d2/pi;

end
