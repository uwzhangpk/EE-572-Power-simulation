clc,clear,close all
% prompt = 'Enter the filename located in the same file location of this driver:';
% str = input(prompt,'s');
[V, Vangle] = termV1("ieee57cdf.txt");
sizeofV = size(V,1);
sizeofAngle = size(Vangle,1);
for i = 1:sizeofV
    fprintf("Final voltage (%d) = %f with angle = %f\n",i,V(i),Vangle(i));
end
s = "ieee57cdf.txt";
fid = fopen(s);
size = textscan(fid,"%s%s%s%f%s%s%s%d",1);
busSize = size{1,8};
base = size{1,4};
fileType = size{1,1};
textBusData = textscan(fid,'%d%11c%s%d%d%d%f%f%f%f%f%f%f%f%f%f%f%f%d',busSize,'headerlines',2);

busType = textBusData{1,6};
finalV = textBusData{1,7};
finalVangle = textBusData{1,8};


errorbusV = zeros(busSize,1);
errorbusVangle = zeros(busSize,1);

for i =1:busSize
    errorbusV(i,1) = abs((V(i,1)-finalV(i,1))/finalV(i,1));
    errorbusVangle(i,1) = abs((Vangle(i,1)-finalVangle(i,1))/finalVangle(i,1));
end

fprintf("max error V = %f\n", max(errorbusV));
fprintf("max error Vangle = %f\n", max(errorbusVangle));
