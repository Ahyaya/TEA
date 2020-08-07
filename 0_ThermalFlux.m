global lonGrid;
global latGrid;
Re=6.371e6;
ThermalFlux=zeros(100,1);			%Record data in size of 100

lonGrid=importdata('lonGrid');
latGrid=importdata('latGrid');
lonGrid=lonGrid*pi/180;latGrid=latGrid*pi/180;	%Convert degree into rad
[lonGrid,latGrid]=meshgrid(lonGrid,latGrid);	%Convert axis vector into matrix

SOL_0=importdata('SOL_2015May10');
t_0=1084147264.184;				%getTT(2034 May 10 @ 12:0:0.000)
SOL_1=importdata('SOL_2015Jun10');
t_1=1086825664.184;				%getTT(2034 Jun 10 @ 12:0:0.000)

Flight=importdata('GCS_SC1');
						%Flight(:,1)=time(s)
						%Flight(:,2)=R(km)
						%Flight(:,3)=longitude(rad,E)
						%Flight(:,4)=latitude(rad,N)
Flight(:,2)=1000*Flight(:,2);			%Convert into Flight(:,2)=R(m)

S_unit=pi*pi*Re*Re/518400*cos(latGrid);		%Get area of Sphere Unit (1440*720)

for pf=1:10					%Loop size depends on flight dataset
t=Flight(pf,1);Rs=Flight(pf,2);lons=Flight(pf,3);lats=Flight(pf,4);
Sigma=f_getSigmaT(Rs,lons,lats);		%Get geometric coefficient
SOL=(t_1-t)/(t_1-t_0)*SOL_0+(t-t_0)/(t_1-t_0)*SOL_1;	%Interpolate SOL
ThermalFlux(pf)=sum(sum(SOL'.*S_unit.*Sigma));
end

save ATF_result ThermalFlux
