function [RespType, CaDur, CaAmount, CaLatency, StartOfResp, EndOfResp]= FourCaResponseTypes_TH(TrueTroughVal,...
    TrueTroughLoc, PeakVal, PeakLoc, x, t, t_inp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Marsa Taheri and Gregory Handy, 2016
% This code was used to simulate the mathematical model of Astrocyte 
% IP3-dependent Ca responses in 2 papers submitted in Nov 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code finds the model's Ca response type, duration, total amount (i.e. 
%area under the curve), latency (i.e. time from the start of IP3 input to 
%the start of the Ca response), and start and end time points of the Ca
%response. 
%
%The possible resulting response types are SP (Single-Peak), MP (Multi-Peak), 
%Plat (Plateau), and LL (Long-Lasting). However, the response might also be 
%identified as too small, too long, or too large to check whether or not 
%it's a reasonable response.
%
%The user needs to input the following information: 
%All the response troughs and peaks (as defined in the plotCa_TH.m
%function) including their values and locations, the actual Ca trace, the 
%corresponding time vector, and the IP3/stimulus input time.
%
%Notes: 1) This code is only valid if the 1st peak is the highest peak;
%otherwise, the code will need to be modified.
%2) The code uses another .m file in some places which is titled 
%'RunParallelFunc_RespType.m'
%
%(Last code update: Feb 19, 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


LL_Dur = 70; %The duration (sec) of Ca elevation that determines whether the 
%signal is LL or not

P=length(PeakVal); T=length(TrueTroughLoc);
if P==0
    display('No Response')
    RespType='NR';
    CaDur=0; CaAmount=0; CaLatency=NaN;
    StartOfResp=0; EndOfResp=0;
else

stepsize=diff(t); stepsize=stepsize(1);
stepsizeScale=stepsize/0.01; %adjusts for different dt value used to generate Ca trace

indx2 = find(t>=t_inp,1);
CaLatIndx = find(x(indx2:end) > x(indx2-1), 1,'first');
b4Resp = x(indx2 + CaLatIndx - 3); %before Ca starts to increase
h = max(x) - b4Resp;%max height


%%%%%%Ca Response Start and End Points:
AboveBaseline=x(2)*1.4; %This defines whether the response started or not

StartIndex = find(x>AboveBaseline,1,'first'); 
time_of_max = PeakLoc(PeakVal==max(PeakVal));
EndIndex = find(x(time_of_max:end)>AboveBaseline,1,'last')+time_of_max-1; 

CaDur=t(EndIndex)-t(StartIndex);
StartOfResp = t(StartIndex);
EndOfResp = t(EndIndex);


%%%%%%Ca Amount:
CaResp = x(StartIndex:EndIndex);
CaAmountTemp = cumsum(CaResp);
CaAmount = CaAmountTemp(end)/length(x)*t(end); %Needs to be length of all x 

%%%%%%Ca Latency:
CaLatency = t(StartIndex) - t_inp;


%%%%%%Ca Response Type Algorithm:
if CaDur>200 %Lasts > 200 sec
    display('Other (Too long)')
    RespType='O (Too long)';
elseif max(PeakVal)<0.4 %if the LARGEST peak is too small (<0.4 microM)
    display('Other (too small)')
    RespType='O (SP/NR)'; %between Single-Peak and No Response
elseif sum(PeakVal>3.5)>0 %if ANY peak is too large (>3.0 or >3.5 microM)
    display('Other (too large)')
    RespType='O (too large)';
    
else    
    LastTrough = EndIndex;
    TrueTroughLoc(1) = StartIndex;
   
    [~, DerivMinLoc]=findpeaks(-diff(x(PeakLoc(1):end)), 'MINPEAKHEIGHT', stepsizeScale*0.0001);
    Decreasing=PeakLoc(1) + DerivMinLoc(1); %[First, Not Max] negative slope of Ca2+ after its peak
    %Find where the deriv. is first close to zero after the Ca2+ decreases
    %the [first, Not most] after the Peak:
    ZeroDerLoc=Decreasing + find((diff(x(Decreasing:end))).^2<stepsizeScale.^2*1e-7, 1,'first');
    if P==1; [bumpVal, ~]=findpeaks(x(ZeroDerLoc:end), 'MINPEAKHEIGHT', 0.1*(PeakVal(1)-b4Resp));
    else bumpVal=[]; end
    %NOTE: the bumpy response is valid & makes sense only if there was one
    %detectable peak (P=1)
    
    countTroughMP = 1;
    MPTrough=[]; PeakVal2Include=[];
    if P>1 && T>1 
        if P>T
            for i=1:T
                if i==T; EndTrough=LastTrough; else EndTrough=TrueTroughLoc(i+1);end
                Cond1 = PeakLoc>TrueTroughLoc(i); Cond2 = PeakLoc<EndTrough;
                PeakVal2Include(i) = max(PeakVal(Cond1 == Cond2)); %find where both conditions hold
            end
        else %P==T, but still they're 2 or more
            PeakVal2Include=PeakVal;
        end
        
        for k=2:T
            HalfHeight = 0.5*(max(PeakVal2Include(k-1), PeakVal2Include(k)) - b4Resp);
            if (TrueTroughVal(k) - b4Resp)<HalfHeight
                MPTrough(countTroughMP) = TrueTroughLoc(k);
                countTroughMP = countTroughMP + 1;
            end
        end
    end
    Segments = [TrueTroughLoc(1), MPTrough, LastTrough]; %Whether "MPTrough" is empty or not
    
    if max(diff(t(Segments))) > LL_Dur
        display('Long-Lasting')
        RespType='LL';
        
    elseif isempty(MPTrough)==0 %There ARE troughs that are low enough to make response MP
        NotMPUnlessBump=1;
        for k=2:T %Checks if h2 (or h1) is > 5% of h1 (or h2) for each pair of heights
            height1 = PeakVal2Include(k-1) - TrueTroughVal(k-1);
            height2 = PeakVal2Include(k) - TrueTroughVal(k);
            AllHeights = [height1, height2]; 
            if max(AllHeights)*0.05<=min(AllHeights)
                display('Multi-Peak')
                RespType='MP';
                return
            end
        end
        if P<3 && NotMPUnlessBump==1 && length(MPTrough)==T-1
            %If there are no more than 2 peaks, and all troughs go below 
            %50% of h, but none had a height of >5% adjacent peak height,
            %then check for two more things- Trough Loc and Distances:
            if (x(MPTrough)-b4Resp)<0.2*HalfHeight
                if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                    display('Long-Lasting')
                    RespType='LL';
                    return
                else
                    display('Single-Peak')
                    RespType='SP';
                    return
                end
            else
                PlatPoint = min(MPTrough, ZeroDerLoc);
                PlatToEnd = LastTrough - PlatPoint;
                Trough1ToPlat = PlatPoint - TrueTroughLoc(1);
                if PlatToEnd > 0.5*Trough1ToPlat %Back to original condition
                    if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                        display('Long-Lasting')
                        RespType='LL';
                        return
                    else
                        display('Plateau')
                        RespType='Plat';
                        return
                    end
                elseif length(DerivMinLoc)>3 %4 or more times
                    if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                        display('Long-Lasting')
                        RespType='LL';
                        return
                    else
                        display('Plateau')
                        RespType='Plat';
                        return
                    end
                else
                    if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                        display('Long-Lasting')
                        RespType='LL';
                        return
                    else
                        display('Single-Peak')
                        RespType='SP';
                        return
                    end
                end
            end
        else 
            RunParallelFunc_TH %uses another .m file
        end
    
    elseif isempty(bumpVal)==0 %There is a bump after the (one) peak [see NOTE above]
        height1 = PeakVal(1) - b4Resp;
        height2 = bumpVal - x(ZeroDerLoc); %Should not use "b4Resp"
        AllHeights = [height1, height2];
        if max(AllHeights)*0.05<=min(AllHeights) && (x(ZeroDerLoc) - b4Resp) < 0.5*h
            display('Multi-Peak')
            RespType='MP';
            return
        else %Now, only the d1 and d2 matter to determine if it's a SP or Plat
            %There's 1 peak and a bump that is too small to be MP, or the "trough" is at >50% of h
            EndOfResponse = LastTrough;
            PlatToEnd = EndOfResponse - ZeroDerLoc;
            Trough1ToPlat = ZeroDerLoc - TrueTroughLoc(1);
            if (x(ZeroDerLoc)-b4Resp)>=0.1*h && PlatToEnd>0.5*Trough1ToPlat
                if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                    display('Long-Lasting')
                    RespType='LL';
                    return
                else
                    display('Plateau')
                    RespType='Plat';
                    return
                end
            elseif length(DerivMinLoc)>3 %4 or more times
                if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                    display('Long-Lasting')
                    RespType='LL';
                    return
                else
                    display('Plateau')
                    RespType='Plat';
                    return
                end
            else
                if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                    display('Long-Lasting')
                    RespType='LL';
                    return
                else
                    display('Single-Peak')
                    RespType='SP';
                    return
                end
            end
            
        end
    
    %From here on, it must be a SP, LL, or Plateau response (though had
    %some SP & Plat cases above too), though not MP:
    elseif isempty(MPTrough) && T>1 %The trough was >50% of adjacent heights
        if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
            display('Long-Lasting')
            RespType='LL';
            return
        else
            display('Plateau')
            RespType='Plat';
            return
        end
        
    elseif P>2
        %Still need to check two things- Distance and TroughLoc:
        EndOfResponse = PeakLoc(2) + find(x(PeakLoc(2):end) > 1.15*x(1), 1, 'last');
        PlatToEnd = EndOfResponse - ZeroDerLoc;
        Trough1ToPlat = ZeroDerLoc - TrueTroughLoc(1);
        %Note: This part is only valid if the 1st peak is the highest peak;
        %otherwise, the code will need to be modified.
        if (x(ZeroDerLoc)-b4Resp)>=0.1*(PeakVal(1)-b4Resp) && PlatToEnd>0.5*Trough1ToPlat
            if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                display('Long-Lasting')
                RespType='LL';
                return
            else
                display('Plateau')
                RespType='Plat';
                return
            end
        else
            if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                display('Long-Lasting')
                RespType='LL';
                return
            else
                display('Single-Peak')
                RespType='SP';
            end
        end
        
        %NOTE: The part below assumes the 1st peak of 2 peaks is larger; If
        %this is not the case, then the code needs to be adjusted to consider 
        %the opposite case too.
    elseif P==2 %&& T==1, at this point this must be true too
        EndOfResponse = PeakLoc(2) + find(x(PeakLoc(2):end) > 1.15*x(1), 1, 'last');
        PlatToEnd = EndOfResponse - ZeroDerLoc;
        Trough1ToPlat = ZeroDerLoc - TrueTroughLoc(1);
        if (x(ZeroDerLoc)-b4Resp)>=0.1*(PeakVal(1)-b4Resp) && PlatToEnd>0.5*Trough1ToPlat
            if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                display('Long-Lasting')
                RespType='LL';
                return
            else
                display('Plateau')
                RespType='Plat';
                return
            end
        else
            if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                display('Long-Lasting')
                RespType='LL';
                return
            else
                display('Single-Peak')
                RespType='SP';
                return
            end
        end
        
    elseif P==1 %The only case remaining
        if length(DerivMinLoc)==1
            if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                display('Long-Lasting')
                RespType='LL';
                return
            else
                display('Single-Peak')
                RespType='SP';
                return
            end
        elseif length(DerivMinLoc)>3 %4 or more times
            if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                display('Long-Lasting')
                RespType='LL';
                return
            else
                display('Plateau')
                RespType='Plat';
                return
            end
        else %When the deriv decreases 2 or 3 times
            diffTrace=diff(x);
            [~,l]= findpeaks((diff(diff(x(Decreasing:LastTrough))).^2),'SORTSTR', 'descend');
            if isempty(l)
                PlatPoint=ZeroDerLoc;
            else
                PlatPointTemp = Decreasing + l(1); %l(1) is max peak, since findpeaks is in descending order
                PlatPoint = min(PlatPointTemp, ZeroDerLoc);
                if abs(diffTrace(PlatPoint))>stepsizeScale*5e-4 %Checks if the deriv at the chosen PlatPoint 
                    %is not close enough to zero
                    PlatPoint = max(PlatPointTemp, ZeroDerLoc);
                end
            end
            EndOfResponse = LastTrough;
            PlatToEnd = EndOfResponse - PlatPoint; 
            Trough1ToPlat = PlatPoint - TrueTroughLoc(1);
            if (x(PlatPoint)-b4Resp)>=0.1*h && PlatToEnd > 0.5*Trough1ToPlat
                if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                    display('Long-Lasting')
                    RespType='LL';
                    return
                else
                    display('Plateau')
                    RespType='Plat';
                    return
                end
            else
                if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now
                    display('Long-Lasting')
                    RespType='LL';
                    return
                else
                    display('Single-Peak')
                    RespType='SP';
                    return
                end
            end
        end
    end
end
end
