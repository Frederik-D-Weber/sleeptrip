function [REMInBurst funName remEpochs]= allNightREMsPrcInMidBurstsFromExactVct(hip, epochLength, exactREMsVct, sampling)


funName = 'allNightAR%IMidB';
%miniEpochSize = 1; % how long our miniepochs should be?
%REMsInMiniEpochs = computeMiniEpochsWithREMs(exactREMsVct, sampling, miniEpochSize);

distanceWithinBurst = 2; % how much time (in sec) can be btween REMs in order to classify them as one burst?
distanceWithinBurst = distanceWithinBurst*sampling;

minimumREMsInBurst = 3; %how many REMs must be in burst, so we say it is a burst


hipLen = length(hip);


if hipLen > length(exactREMsVct)/sampling/epochLength
    
   fprintf('Vector is shorter than hipnogram!!!\n');     
   hipLen = floor(length(exactREMsVct)/sampling/epochLength);
   
end
%hipLen


remEpochs = 0;
allREMsNb = 0;

allREMsInBurst = 0;

previousREMPlace = -1000;
burstIsOver = 0;

%duringREM = 0;
%FollowMinis = 1000;


%minimumREMsInBurst
nbOfREMsInCurrentBurst = 0;


for i = 1:1:hipLen

    if hip(i) == 5


        placeInExtVctBeg = (i-1)*sampling*epochLength + 1;
        placeInExtVctEnd = i*sampling*epochLength;


        currentFragment = exactREMsVct(placeInExtVctBeg:placeInExtVctEnd);
        %REMsInCurrFrag = 0;
        duringREM = 0;



        for jj = 1:1:length(currentFragment)


            if currentFragment(jj) > 0

                if ~duringREM 
                    
                    % its a new REM

                    duringREM = 1;
                    allREMsNb = allREMsNb + 1;


                    if jj + (i-1)*sampling*epochLength - previousREMPlace <= distanceWithinBurst 

                        % we are in fact in burst
                        nbOfREMsInCurrentBurst = nbOfREMsInCurrentBurst + 1;
                        if nbOfREMsInCurrentBurst == 1  % if result is one, it means that REM before should also be included
                            nbOfREMsInCurrentBurst = 2;
                        end
                        
                        burstIsOver = 0;
                        
                        % is burst long enough?
                        
                        
%                         allREMsInBurst = allREMsInBurst + 1;
% 
%                         if lastInBurst == 0
% 
%                             allREMsInBurst = allREMsInBurst + 1;
% 
%                         end
% 
%                         lastInBurst = 1;

                    else % this REM is not close enough to previous one

%                        burstIsOver = 1;
                        nbOfREMsInCurrentBurst = 0;

                    end

                end % if ~duringREM 

            elseif duringREM == 1

                duringREM = 0;                
                previousREMPlace = jj + (i-1)*sampling*epochLength;

%                 if burstIsOver == 1
% 
%                     nbOfREMsInCurrentBurst = 0;
%                     
%                 end

                
                if nbOfREMsInCurrentBurst == minimumREMsInBurst % include all burst

                    allREMsInBurst = allREMsInBurst + nbOfREMsInCurrentBurst;
                    
                elseif nbOfREMsInCurrentBurst > minimumREMsInBurst

                    allREMsInBurst = allREMsInBurst + 1;

                end



            end % if currentFragment(jj) > 0

            
        end % for jj = 1:1:length(currentFragment)

      
    end % if hip(i) == 5
    
end % for i = 1:1:hipLen


if allREMsNb
   
    REMInBurst = 100*(allREMsInBurst / allREMsNb);
    
else
    
    REMInBurst = nan;
    
end


