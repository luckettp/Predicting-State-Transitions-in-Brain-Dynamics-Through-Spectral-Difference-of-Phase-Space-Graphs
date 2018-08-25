
function a = Copy_of_FWEIGtest(RAW)
tic
DIM = 2;
SYM = 3;
LAG = 50;
MU = 46;
HW = 28;
CS = 20676;
BC = 9;
socFlag = 5;
[rawSize,~] = size(RAW);
filtered = zeros((CS-2*HW),(round(rawSize/CS))-1);
BCDM = zeros(nchoosek(BC,2),1);
NPS = (CS-2*HW)-((DIM-1)*LAG);
NPC = NPS-MU;
nodes = zeros(NPS,DIM);
nodesID = zeros(NPS,1);
for i = 1:(rawSize(1,1)/CS)
    filtered(:,i) = ArtFil(RAW(((i-1)*CS+1):((i)*CS),1),HW);
end
mx = max(filtered(:,1));
mn = min(filtered(:,1));
scale = (SYM/((mx-mn)));
filtered = filtered - mn;
filtered = max(0,min(SYM-1,int8(floor(scale*filtered))));
[~,A] = size(filtered);
eigVal = zeros(A,SYM^DIM);
TCDM = zeros(A,1);
adjSize = SYM^DIM;
for j = 1:A
    for i = 1:NPS
        nodes(i,DIM:-1:1) = filtered(i:LAG:i+(LAG*DIM-1),j);
    end
    tmp = unique(nodes,'rows');
    [s1,~] = size(tmp);
    for i = 1:NPS
        for k = 1:s1
            if(nodes(i,:) == tmp(k,:))
                nodesID(i,1)=k;
                break;
            end
        end
    end
    links = [nodesID(1:NPC,1),nodesID(MU+1:NPS,1)];
    links = unique(links,'rows');
    ADJM = zeros(adjSize,adjSize);
    for i = 1:size(links)
        if(links(i,1)~=links(i,2))
            ADJM(links(i,1),links(i,2)) = 1;
            ADJM(links(i,2),links(i,1)) = 1;
        end
    end
    tmpEIG = eig(ADJM)';
    eigVal(j,1:length(tmpEIG)) = tmpEIG;
end
count = 1;
for i = 1:BC-1
    for j = (i+1):BC
        BCDM(count,1) = i;
        BCDM(count,2) = j;
        BCDM(count,3) = pdist([eigVal(i,:);eigVal(j,:)]);
        count = count + 1;
    end
end
EIGMeanD = mean(BCDM(:,3));
EIGSTDD =std(BCDM(:,3));
for i = 1:BC
    u = 0;
    j = 1;
    for z = 1:size(BCDM)
        if ((BCDM(j,1) == i)||(BCDM(j,2) == i))
            u = u + BCDM(j,3);
        end
        j = j+1;
    end
    TCDM(i,1) = u;
end
for j = BC+1:A
    for k = 1:BC
        TCDM(j,1) = TCDM(j,1) + pdist([eigVal(j,:);eigVal(k,:)]);
    end
end
% Normalize data
TCDM(1:BC,1) = TCDM(1:BC,1)/(BC-1);
TCDM((BC+1):A,1) = TCDM((BC+1):A,1)/BC;
TCDM(:,1) = (TCDM(:,1)-EIGMeanD)/EIGSTDD;

mi = min(TCDM);
mx = max(TCDM);
TCDM = TCDM-mi;
TCDM = TCDM./(mx-mi);

SOAT = zeros(A,1);
threshold = 0.366914;
status = 0;
FWCS = 0;
for i = 1:A
    if(TCDM(i,1)>threshold)
        SOAT(i,1) = 1;
    end
end
count = 0;
for i = 1:size(TCDM)
    if(SOAT(i,1)>=1)
        count = count + 1;
    else
        count = 0;
    end
    if(count >= socFlag)
        status = 1;
        FWCS = i;
        break;
    end
end
a = [status,FWCS];
end

function VF = ArtFil(ICUT,W)
NP = size(ICUT);
NN= 3*(3*W*W+3*W-1);
RD = (2*W+1)*(4*W*W+4*W-3);
N1 = W+1;
CF1 = NN/RD;
CF2 = 15.0D0/RD;
F0  = 0.0;
F1  = 0.0;
F2  = 0.0;
VU = ICUT;
for I = N1:(NP-W)
    if(I==N1)
        for J =-W:W
            F0 =F0+VU(I+J);
            F1 =F1+VU(I+J)*J;
            F2 =F2+VU(I+J)*J*J;
        end
    else
        G0  =F0-VU(I-N1)+VU(I+W);
        G1  =F1+VU(I-N1)*N1+VU(I+W)*W-F0;
        G2  =F2-VU(I-N1)*N1*N1+VU(I+W)*W*W-2*F1+F0;
        F0  =G0;
        F1  =G1;
        F2  =G2;
    end
    VA=CF1*F0-CF2*F2;
    VF(I-W)=VU(I)-VA;
end
end

    
