function [TrueTroughVal, TrueTroughLoc, PeakVal, PeakLoc] = PeaksAndTroughs_TH(x, threshPk, threshDist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Marsa Taheri and Gregory Handy, 2016
% This code was used to simulate the mathematical model of Astrocyte 
% IP3-dependent Ca responses in 2 papers submitted in Nov 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Uses the 'threshPk' value to determine all the Peaks seen in trace 'x'.
%Also finds the last trough before each peak and uses the 'threshDist' to
%determine whether it is a True Trough or not.
%The output is all the True Trough values and locations as well as all the
%Peak values and locations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PeakVal, PeakLoc]=findpeaks(x, 'MINPEAKHEIGHT', threshPk);
[~, locn]=findpeaks(-x);

count=1; L=length(PeakLoc);
troughLoc=zeros(L); troughVal=zeros(L);
if L==0
    TrueTroughVal=[]; TrueTroughLoc=[];
    PeakVal=[]; PeakLoc=[];
else
    for i=1:L
        if isempty(locn) || locn(1)>PeakLoc(1)%If there's No trough before 1st peak
            startResp = find(x>1.05*x(1), 1, 'first');
            if size(locn,2)==1
                locn=[startResp; locn];
            else
                locn=[startResp, locn];
            end
        end
        troughLoc(i)=locn(find(locn<PeakLoc(i), 1, 'last'));
        troughVal(i)=x(troughLoc(i));
        if PeakVal(i) > (troughVal(i) + threshDist)
            TrueTroughLoc(count)=troughLoc(i);
            TrueTroughVal(count)=troughVal(i);
            count=count+1;
        end
    end
end
