clear all;
%%%%%%%%%%%%%%%%%%%%%
% MAIN PARAMETERS
%%%%%%%%%%%%%%%%%%%%%

% INPUT FOLDER should contain .txt files exported from ImageJ
% with 2 columns: 
% 1) distance along the axon in pixels
% 2) intensity values

%number of averaging points for smoothing 
%bigger values correspond to more smoothing 
nAverPoints = 50; 
%pixel resolution. number of pixels in one micron.
nPixMkm = 15.8;
%threshold values of normalized intensity to determine 
%Start and End values (in range from 0 to 1)
dThreshold = 0.33;
%including tails or not in case of center moving average
bTails = true;

% END OF PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%


% open text file with intensity profile
[filename,filepath] = uigetfile({'*.txt';'*.*'},'Select txt file with distance and intensity profile...');
%if cancel button is not pressed
if ~isequal(filename, 0)
        disp('Script started.');
        disp('Please wait......');
        filenamein = sprintf('%s%s',filepath,filename);    
        %read intensity trace
        inttrace = dlmread(filenamein,'\t');
        sz = size(inttrace);
        %total number of points
        nTotPoints = sz(1);
        
        %array containing results
        norm_int = zeros(nTotPoints,2,'double');
        nonnorm_int = zeros(nTotPoints,2,'double');
        norm_int(:,1)=inttrace(:,1)*(1/nPixMkm);
        nonnorm_int(:,1)=norm_int(:,1);
        nonnorm_int(:,2)=inttrace(:,2);
        
        if (rem(nAverPoints,2)==0)
            nAverPoints=nAverPoints+1;
        end
        nHalf = (0.5*(nAverPoints-1));
        
        %SMOOTHING
        %including 'tails' in the beginning and end of the track
        %depending on the flag
        if(bTails)
            %beginning
            for i=1:nHalf
                norm_int(i,2) = mean(inttrace(1:(i+nHalf),2));
            end
            %end
            for i=(nTotPoints-nHalf+1):nTotPoints
                norm_int(i,2) = mean(inttrace((i-nHalf):nTotPoints,2));
            end
        end
        %smoothing middle part
        for i=(nHalf+1):(nTotPoints-nHalf)
            norm_int(i,2) = mean(inttrace((i-nHalf):(i+nHalf),2));
        end
        
        %NORMALIZATION
        %min intensity
        minInt = min(norm_int(:,2));
        %max-min intensity, normalization factor
        dSpan = max(norm_int(:,2)) - minInt;
        %normalize intensity
        norm_int(:,2) = (norm_int(:,2) - minInt) / dSpan;
        
                    %SEARCH FOR REFERENCE POINTS
                    %Max
                    [~, MaxInd] = max(norm_int(:,2));
                    MaxVal = norm_int(MaxInd,1);
                    
                    %Start 1 (from left most point)
                    i=0;
                    dSpan = 0;
                    while dSpan < dThreshold
                        i = i + 1;
                        dSpan = norm_int(i, 2);
                    end
                    if i==1
                        i = 2;
                    end
                    StartVal1 = norm_int(i-1, 1);
                    StartInd1 = i-1;

                    %Start 2 (from maximum)
                    i=MaxInd;
                    dSpan = MaxVal;
                    while (dSpan > dThreshold) && (i>1)
                        i = i - 1;
                        dSpan = norm_int(i, 2);
                    end
                    StartVal2 = norm_int(i, 1);
                    StartInd2 = i;
                    
                    %End 1 (from right most point)
                    i=nTotPoints+1;
                    dSpan = 0;
                    while dSpan < dThreshold
                        i = i - 1;
                        dSpan = norm_int(i, 2);
                    end
                    EndVal1 = norm_int(i, 1);
                    End1Ind = i;
                    
                    %End 2 (from maximum)
                    i=MaxInd;
                    dSpan = 2;
                    while (dSpan > dThreshold) && (i < nTotPoints)
                        i = i + 1;
                        dSpan = norm_int(i, 2);
                    end
                    EndVal2 = norm_int(i-1, 1);
                    End2Ind = i-1;
                    
                    

                    %Calculate AIS lengths
                    LengthVal1 = EndVal1 - StartVal1;
                    LengthVal2 = EndVal2 - StartVal1;                                            
                    LengthVal3 = EndVal1 - StartVal2;
                    LengthVal4 = EndVal2 - StartVal2;
             
                    %INTENSITY QUANTIFICATION
                    bLeftTail1 = sum(nonnorm_int(1:StartInd1,2))/nonnorm_int(StartInd1,1);
                    bLeftTail2 = sum(nonnorm_int(1:StartInd2,2))/nonnorm_int(StartInd2,1);
                    bMeanInt11 = sum(nonnorm_int((StartInd1+1):End1Ind,2))/(nonnorm_int(End1Ind,1)-nonnorm_int(StartInd1,1));
                    bMeanInt21 = sum(nonnorm_int((StartInd2+1):End1Ind,2))/(nonnorm_int(End1Ind,1)-nonnorm_int(StartInd2,1));
                    bMeanInt12 = sum(nonnorm_int((StartInd1+1):End2Ind,2))/(nonnorm_int(End2Ind,1)-nonnorm_int(StartInd1,1));
                    bMeanInt22 = sum(nonnorm_int((StartInd2+1):End2Ind,2))/(nonnorm_int(End2Ind,1)-nonnorm_int(StartInd2,1));
                    bRightTail1 = sum(nonnorm_int((End1Ind+1):nTotPoints,2))/(nonnorm_int(nTotPoints,1)-nonnorm_int(End1Ind,1));
                    bRightTail2 = sum(nonnorm_int((End2Ind+1):nTotPoints,2))/(nonnorm_int(nTotPoints,1)-nonnorm_int(End2Ind,1));
                   
                    
        
        %Make plot
        plot((norm_int(:,1)),(norm_int(:,2)),'Color','b','LineWidth',2);
        line([StartVal1,StartVal1],[0,1],'LineStyle',':', 'Color','b');
        line([EndVal1,EndVal1],[0,1],'Linestyle',':','Color','b');
        line([MaxVal,MaxVal],[0,1],'LineStyle','--', 'Color','r');
        line([StartVal2,StartVal2],[0,1],'LineStyle','--', 'Color','g');
        line([EndVal2,EndVal2],[0,1],'Linestyle','--','Color','g');
        hleg1=legend('normalized intensity','start from left', 'end from right','max','start from max', 'end from max');


        %saving results data
        filenamein =  strcat(filenamein(1:length(filenamein)-4),'_results.xls');
        
        headers = {'length,mkm', 'normalized intensity', 'Start from left (1), mkm', 'Start from max value(2), mkm', 'Max, mkm', 'End from right end (1), mkm','End from max value (2), mkm', 'AISlength (End1-Start1), mkm', 'AISlength (End2-Start1), mkm', 'AISlength (End1-Start2), mkm', 'AISlength (End2-Start2), mkm', '# of smooth points', 'intensity threshold', 'resolution (pixels in mkm)', 'Tails', 'Int_Left_tail1, a.u./mkm', 'Int_Left_tail2, a.u./mkm','Int_Middle_(Start1_End1)','Int_Middle_(Start2_End1)', 'Int_Middle_(Start1_End2)', 'Int_Middle_(Start2_End2)', 'Int_Right_tail1, a.u./mkm', 'Int_Right_tail2, a.u./mkm'};
        xlswrite(filenamein,headers,'Sheet1','A1');
        params = [StartVal1, StartVal2, MaxVal, EndVal1, EndVal2, LengthVal1, LengthVal2, LengthVal3, LengthVal4, nAverPoints, dThreshold, nPixMkm, bTails, bLeftTail1, bLeftTail2, bMeanInt11, bMeanInt21, bMeanInt12, bMeanInt22, bRightTail1, bRightTail2];
        xlswrite(filenamein,params,'Sheet1','C2');
        xlswrite(filenamein,norm_int,'Sheet1','A2');
        disp('Done.');
end