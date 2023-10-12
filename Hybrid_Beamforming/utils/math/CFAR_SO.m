% input
%   inData: data sieze [1, N]
%
%output
%   N_obj: number of objects detected
%   Ind_obj: 2D bin index of the detected object
%   noise_obj: noise value at the detected point after integration
% 
function [N_obj, Ind_obj, noise_obj, CFAR_SNR] = CFAR_SO(inData, cellNum, gapNum, K0, maxEnable)

gaptot = gapNum + cellNum;
discardCellLeft = 10;
discardCellRight = 10;
N_obj=0;
Ind_obj=[];
noise_obj = [];
CFAR_SNR = [];

M_samp=length(inData);

vec = inData(discardCellLeft+1:M_samp-discardCellRight);
vecLeft = vec(1:(gaptot));
vecRight = vec(end-(gaptot)+1:end);
vec = [vecLeft vec vecRight];  % boundary filling

for j=1:(M_samp-discardCellLeft-discardCellRight)
    cellInd=[j-gaptot: j-gapNum-1 j+gapNum+1:j+gaptot];
    cellInd=cellInd + gaptot;
    cellInda=[j-gaptot: j-gapNum-1];
    cellInda=cellInda + gaptot;
    cellIndb=[j+gapNum+1 : j+gaptot];
    cellIndb=cellIndb + gaptot;
    
    cellave1a =sum(vec(cellInda))/(cellNum);
    cellave1b =sum(vec(cellIndb))/(cellNum);
    cellave1 = min(cellave1a,cellave1b);
    if maxEnable==1 %check if it is local maximum peak
        maxInCell = max(vec([cellInd(1):cellInd(end)]));
        if (vec(j+gaptot)>K0*cellave1 && ( vec(j+gaptot)>=maxInCell))
            N_obj=N_obj+1;
            Ind_obj(N_obj)=j+discardCellLeft;
            noise_obj(N_obj) = cellave1; %save the noise level
            CFAR_SNR(N_obj) = vec(j+gaptot)/cellave1;
        end
    else
        if vec(j+gaptot)>K0*cellave1
            N_obj=N_obj+1;
            Ind_obj(N_obj)=j+discardCellLeft;
            noise_obj(N_obj) = cellave1; %save the noise level
            CFAR_SNR(N_obj) = vec(j+gaptot)/cellave1;
        end
    end        
end

%get the noise variance for each antenna
for i_obj = 1:N_obj
    ind_range = Ind_obj;
    if ind_range<= gaptot
        %on the left boundary, use the right side samples twice
        cellInd=[ind_range+gapNum+1:ind_range+gaptot ind_range+gapNum+1:ind_range+gaptot];
    elseif ind_range>=M_samp-gaptot+1
        %on the right boundary, use the left side samples twice
        cellInd=[ind_range-gaptot: ind_range-gapNum-1 ind_range-gaptot: ind_range-gapNum-1];
    else
        cellInd=[ind_range-gaptot: ind_range-gapNum-1 ind_range+gapNum+1:ind_range+gaptot];
        
    end
    
end



end