function prob2_2_gs5660()
lat=[.00001; 10; 45; 70; 89.99999]*pi/180;
lon=-83*pi/180;
h=[1;5;10;20;50;100;300]*1000;
[a,b,e2]=refell('grs80');
stock=zeros(size(lat,1),size(h,1)*3);
stockbor=zeros(size(lat,1),size(h,1)*2);
for j=1:size(h,1)
    for i=1:size(lat,1)    
        xyz= blh2xyz([lat(i), lon, h(j)], a, e2);
        [latb,lonb,hb]=xyz2ell_borkowski(xyz(1,1),xyz(1,2),xyz(1,3),a,b,e2);
        [blhh,iter] = xyz2blh_Hirvonen(xyz, a, e2, 1e-4/3600);
        stock(i,j*3-2)=iter;
        stock(i,j*3-1)=blhh(1,3)-h(j,1);
        stock(i,j*3)=blhh(1,1)-lat(i,1);
        stockbor(i,j*2-1)=latb-lat(i,1);
        stockbor(i,j*2)=hb-h(j,1);
    end
end
format longg
disp(stock)
disp(stockbor)