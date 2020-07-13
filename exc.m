clear all; close all; clc;
s="ieee14.txt";
fid = fopen(s);
size = textscan(fid,"%s%s%s%f%s%s%s%d",1);
busSize = size{1,8};
base = size{1,4};
fileType = size{1,1};



% 6   Type (I) *
%                  0 - Unregulated (load, PQ)
%                  1 - Hold MVAR generation within voltage limits, (PQ)
%                  2 - Hold voltage within VAR limits (gen, PV)
%                  3 - Hold voltage and angle (swing, V-Theta) (must always
%                       have one)
% 7   Final voltage, p.u. (F) *
% 8   Final angle, degrees (F) *
% 9   Load MW (F) *
% 10   Load MVAR (F) *
% 11   Generation MW (F) *
% 12   Generation MVAR (F) *
% 13   Base KV (F)
% 14   Desired volts (pu) (F) (This is desired remote voltage if
%                 this bus is controlling another bus.
% 15   Maximum MVAR or voltage limit (F)
% 16   Minimum MVAR or voltage limit (F)
% 17 Shunt conductance G (per unit) (F) *
% 18 Shunt susceptance B (per unit) (F) *
% 19 bus number
textBusData = textscan(fid,'%d%11c%s%d%d%d%f%f%f%f%f%f%f%f%f%f%f%f%d',busSize,'headerlines',2);

% busType = textscan(fid,'%*18c%d',busSize,'headerlines',2);
busType = textBusData{1,6};
finalV = textBusData{1,7};
finalVangle = textBusData{1,8};
loadP = textBusData{1,9};
loadP=loadP/base;
loadQ = textBusData{1,10};
loadQ=loadQ/base;
GenP = textBusData{1,11};
GenP =GenP/base;
GenQ = textBusData{1,12};
GenQ =GenQ/base;

busB = textBusData{1,18};

Qmax = textBusData{1,15};
Qmax = Qmax/base;
Qmin = textBusData{1,16};
Qmin = Qmin/base;
fixedVbus = textBusData{1,14};

lineSkipToBranch = 4+busSize;
lineSkipToBranchSize = lineSkipToBranch-1;

fid = fopen(s);
brachDataLine = textscan(fid,"%s%s%s%d",1,'headerlines',lineSkipToBranchSize);
branchSize = brachDataLine{1,4};

fid = fopen(s);
textBranchData = textscan(fid,'%d%d%d%d%d%d%f%f%f%d%d%d%d%d%f%f%f%f%f%f%f',branchSize,'headerlines',lineSkipToBranch);
busFrom = textBranchData{1,1};
busTo = textBranchData{1,2};
branchR =textBranchData{1,7};
branchX =textBranchData{1,8};
branchB =textBranchData{1,9};

Ybus = zeros(busSize,busSize);

fclose(fid);


%find Y(i,j)
for i =1:branchSize
    yi = busFrom(i,1);
    yj = busTo(i,1);
    Ybus(yi,yj) = -1/(branchR(i,1)+(branchX(i,1)*1j));
    Ybus(yj,yi) = Ybus(yi,yj);
end
% count = 0;
%
% for i = 1:14
%     for j=1:14
%
%         if(Ybus(i,j)~=0)
%             count=count+1;
%         end
%     end
% end


% find Y ii
for i = 1:busSize
    branchShuntB = 0;
    for k = 1:branchSize
        if(busFrom(k,1)==i||busTo(k,1)==i)
            branchShuntB = (1/2)*((branchB(k,1)*1j)) + branchShuntB;
        end
        
        
    end
    %Thus we get branchShuntB
    
    
    busShuntB = 1j*busB(i,1);
    
    for k = 1:busSize
        if(k~=i)
            Ybus(i,i) = -Ybus(i,k) + Ybus(i,i);
        end
        
    end
    
    Ybus(i,i) = busShuntB + branchShuntB + Ybus(i,i);
end

%Ybus completed

%Construct YbusAngle
YbusAngle = zeros(busSize,busSize);

for i = 1:busSize
    for j = 1:busSize
        YbusAngle(i,j) = angle(Ybus(i,j));
%         
%         YbusAngle(i,j) = atand(imag(Ybus(i,j)/real(Ybus(i,j))));
%        YbusAngleCompare(i,j) = atand(imag(Ybus(i,j)/real(Ybus(i,j))));
        
        
    end
end
 Ybus = abs(Ybus);
%  for i = 1:busSize
%     for j = 1:busSize
%         Ybus(i,j) = Ybus(i,j)*cos(YbusAngle(i,j))+1j*Ybus(i,j)*sin(YbusAngle(i,j));
% %         
% %         YbusAngle(i,j) = atand(imag(Ybus(i,j)/real(Ybus(i,j))));
% %        YbusAngleCompare(i,j) = atand(imag(Ybus(i,j)/real(Ybus(i,j))));
%         
%         
%     end
% end
% YbusAngle = rad2deg(YbusAngle);
%Find the swing bus
PVType=zeros(busSize,1);
PQType=zeros(busSize,1);
validPV=1;
validPQ=1;
for i = 1:busSize
    if(busType(i,1)==3)
        swingBus = i;
    end
    
    
    
end
% finalV(swingBus,1) = 1;
% finalVangle(swingBus,1) = 0;
for i = 1:busSize
    if(fixedVbus(i,1) ~=0&&i~=swingBus)
        PVType(validPV,1) = i;        
        validPV=validPV+1;
        
        finalV(i,1) = fixedVbus(i,1);
        finalVangle(i,1) = 0;
    elseif(i~=swingBus)
        PQType(validPQ,1) = i;
        validPQ=validPQ+1;
        
        finalV(i,1) = 1;
        finalVangle(i,1) = 0;
    else
        continue
    end
end
staticPVType = PVType;
staticValidPV = validPV;

%adjust the swing bus voltage and magnitude to 1 and 0
% finalV(swingBus,1) = 1;
% finalVangle(swingBus,1) = 0;
%All set up for calculation now

e = 0.001;
mismatch = e+1;
PBus= zeros(busSize,1);
QBus = zeros(busSize,1);
%PBus = P generator - P load
for i = 1:busSize
    PBus(i,1) = GenP(i,1) - loadP(i,1);
end
%QBus = Q generator - Q load
for i = 1:busSize
    QBus(i,1) = GenQ(i,1) - loadQ(i,1);
end
%Starts calculation

P = zeros(busSize,1);
Q = zeros(busSize,1);

%construct x
firstHalf = zeros(busSize-1,1);
secondHalf = ones(busSize-validPV,1);
x = [firstHalf;secondHalf];


preQ = Q;

while mismatch > e
    
    for i = 1:busSize
        sum = 0;
        for j = 1:busSize
            sum = sum +finalV(j,1)*Ybus(i,j)*cos((finalVangle(i,1)-...
                finalVangle(j,1)-YbusAngle(i,j)));
        end
        P(i,1) = finalV(i,1)*sum-PBus(i,1);
        
        
        
    end
    
    %second Q from 1 to busSize
    for i = 1:busSize
        sum = 0;
        for j = 1:busSize
            sum = sum +finalV(j,1)*Ybus(i,j)*sin((finalVangle(i,1)-...
                finalVangle(j,1)-YbusAngle(i,j)));
        end
        
        
%         preQ(i,1) = finalV(i,1)*sum-QBus(i,1);
        preQ(i,1) = finalV(i,1)*sum;
        if(Qmax(i,1)~=0&&Qmin(i,1)~=0)
            check = checkExceed(preQ(i,1),Qmax(i,1),Qmin(i,1));
            if(check == 0)
                Q(i,1) = preQ(i,1)-QBus(i,1);
                
            elseif(check == 1)
                PVType(staticValidPV,1) = i;
                staticValidPV=staticValidPV+1;
                Q(i,1) = Qmax(i,1)-QBus(i,1);
                validPV=validPV+1;
                fixedVbus(i,1) = finalV(i,1);
                staticPVType = PVType;
%                 Q(i,1)=[];
%                 x(size(firstHalf))
            elseif (check == 2)
                PVType(staticValidPV,1) = i;
                
                Q(i,1) = Qmin(i,1)-QBus(i,1);
                validPV=validPV+1;
                fixedVbus(i,1) = finalV(i,1);
                staticPVType = PVType;
                staticValidPV = validPV;
%                 Q(i,1)=[];
            end
        else
            Q(i,1) = preQ(i,1)-QBus(i,1);
            
        end
        
    end
    
    
    %Compute Jacobian Matrix
    
    J1 = zeros(busSize,busSize);
    J2 = zeros(busSize,busSize);
    J3 = zeros(busSize,busSize);
    J4 = zeros(busSize,busSize);
    %Looking for J1
    for k = 1:busSize
        for n = 1:busSize
            if(k==swingBus || n == swingBus)
                continue
            elseif(k~= n)
                J1(k,n) = finalV(k,1)* Ybus(k,n)*finalV(n,1)...
                    *sin((finalVangle(k,1)-finalVangle(n,1)- YbusAngle(k,n)));
            elseif(k==n)
                sum = 0;
                for i = 1:busSize
                    if(i~=k)
                        sum = sum +Ybus(k,i)*finalV(i,1)...
                            *sin((finalVangle(k,1)-finalVangle(i,1)-YbusAngle(k,i)));
                    end
                end
                J1(k,k) = -finalV(k,1)*sum;
            end
            
        end
    end
    
    %looking for J2
    for k = 1:busSize
        for n = 1:busSize
            if(k==swingBus || n == swingBus)
                continue
            elseif(k~= n)
                J2(k,n) = finalV(k,1)* Ybus(k,n)...
                    *cos((finalVangle(k,1)-finalVangle(n,1)- YbusAngle(k,n)));
            elseif(k==n)
                sum = 0;
                for i = 1:busSize
                    sum = sum +Ybus(k,i)*finalV(i,1)...
                        *cos((finalVangle(k,1)-finalVangle(i,1)-YbusAngle(k,i)));
                end
                J2(k,k) = finalV(k,1)*Ybus(k,k)*cos((YbusAngle(k,k)))+sum;
            end
            
        end
    end
    
    %looking for J3
    for k = 1:busSize
        for n = 1:busSize
            if(k==swingBus || n == swingBus)
                continue
            elseif(k~= n)
                J3(k,n) = -finalV(k,1)* Ybus(k,n)*finalV(n,1)...
                    *cos((finalVangle(k,1)-finalVangle(n,1)- YbusAngle(k,n)));
            elseif(k==n)
                sum = 0;
                for i = 1:busSize
                    if(i~=k)
                       sum = sum +Ybus(k,i)*finalV(i,1)...
                         *cos((finalVangle(k,1)-finalVangle(i,1)-YbusAngle(k,i)));
                    end
                end
                J3(k,k) = finalV(k,1)*sum;
            end
            
        end
    end
    
    %looking for J4
    for k = 1:busSize
        for n = 1:busSize
            if(k==swingBus || n == swingBus)
                continue
            elseif(k~= n)
                J4(k,n) = finalV(k,1)* Ybus(k,n)...
                    *sin((finalVangle(k,1)-finalVangle(n,1)- YbusAngle(k,n)));
            elseif(k==n)
                sum = 0;
                for i = 1:busSize
                    sum = sum +Ybus(k,i)*finalV(i,1)...
                        *sin((finalVangle(k,1)-finalVangle(i,1)-YbusAngle(k,i)));
                end
                J4(k,k) = -finalV(k,1)*Ybus(k,k)*sin((YbusAngle(k,k)))+sum;
            end
            
        end
    end
%     
    if(swingBus<PVType(1,1))
        J1(:,swingBus) = [];
        J1(swingBus,:) = [];
        
        J2(:,swingBus) = [];
        J2(swingBus,:) = [];
        J3(:,swingBus) = [];
        J3(swingBus,:) = [];
        J4(:,swingBus) = [];
        J4(swingBus,:) = [];
        
        P(swingBus,:) = [];
        Q(swingBus,:) = [];
        PVType = PVType -1;
    else
        J1(:,swingBus) = [];
        J1(swingBus,:) = [];  
        P(swingBus,:) = [];
        for i = 1:validPV
            if((swingBus>PVType(i,1))&&(swingBus<PVType(i+1,1)))
                PVType = [PVType(1:i,1);swingBus;PVType(i+1:end,1)];
                validPV = validPV+1;
                break
            end
        end
        if(PVType(i,1)==PVType(validPV,1))
            PVType(validPV+1,1) = swingBus;
            validPV = validPV+1;
        end
    end
    
 
    process = 1;
    count = 0;
    while(PVType(process,1)>0)

       
        J2(:,PVType(process,1))=[];
        J3(PVType(process,1),:)=[];
        J4(:,PVType(process,1))=[];
        J4(PVType(process,1),:)=[];
        
        
        Q(PVType(process,1),:)=[];
        
        PVType = PVType -1;
       count = count+1;
        
        process=process+1;
    end
    
    
    J = [J1 J2;
         J3 J4];
    %Jacobian Matrix completed
     
    f=[P;Q];
    
    x = x - inv(J)*f;
%     
%     firstHalf = x(1:busSize-1,1);       
%     secondHalf = x(busSize:end,1);
%     firstHalf = [firstHalf(1:swingBus-1,1);0;firstHalf(swingBus:end,1)];
%     secondHalf = [secondHalf(1:swingBus-1,1);1;secondHalf(swingBus:end,1)];
%     for i = 1:busSize
%         finalVangle(i,1) = firstHalf(i,1);
%     end
%     for i = 1:busSize
%         if (i==2)
%             continue
%         end
%         finalV(i,1) = (secondHalf(i,1));
%     end
 

    firstHalf = x(1:busSize-1);
    firstHalf = [firstHalf(1:swingBus-1,1);0;firstHalf(swingBus:end,1)];
    for i = 1:busSize
        finalVangle(i,1) = firstHalf(i,1);
    end
    
    secHalf = x(busSize:end,1);
    secSize = length(secHalf);
    
    for i = 1:secSize
        finalV(PQType(i,1),1) = secHalf(i,1);
    end
    PVType = staticPVType;
    validPV = staticValidPV;
    mismatch = max(abs((f)));
end

finalVangle = rad2deg(finalVangle);

for i = 1:busSize
    fprintf("Final voltage (%d) = %f with angle = %f\n",i,finalV(i),finalVangle(i));
end

fprintf("bus8Q %fMVAR\n",preQ(8,1)*100);


function check = checkExceed(P,Pmax,Pmin)
check = 0;
if P>=Pmax 
    check = 1;
elseif P <= Pmin
    check =2;
end




end
