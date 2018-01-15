%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Marsa Taheri and Gregory Handy, 2016
% This code was used to simulate the mathematical model of Astrocyte 
% IP3-dependent Ca responses in 2 papers submitted in Nov 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The script here is used around line 188 of "FourCaResponseTypes.m" 
%function and is essentially the same as the end of the code in
%"FourCaResponseTypes.m"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(bumpVal)==0 %There is a bump after the (one) peak
    height1 = PeakVal(1) - b4Resp;
    height2 = bumpVal - x(ZeroDerLoc); %Should not use "b4Resp"
    AllHeights = [height1, height2];
    if max(AllHeights)*0.05<=min(AllHeights) && (x(ZeroDerLoc) - b4Resp) < 0.5*h
        display('Multi-Peak')
        RespType='MP';
        return
    else %Now, only the d1 and d2 matter to determine if it's a SP or Plat
        %There's 1 peak and a bump that is too small to be MP OR the "trough" is at >50% of h
        EndOfResponse = LastTrough1; 
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
            if t(LastTrough)-t(TrueTroughLoc(1)) > LL_Dur %since the whole response is one segment now!
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
        [~,l]= findpeaks((diff(diff(x(Decreasing:LastTrough1))).^2),'SORTSTR', 'descend');
        if isempty(l)
            PlatPoint=ZeroDerLoc;
        else
            PlatPointTemp = Decreasing + l(1); %l(1) is max peak since findpeaks in descending order
            PlatPoint = min(PlatPointTemp, ZeroDerLoc);
            if abs(diffTrace(PlatPoint))>5e-4 %Checks if the deriv at the chosen PlatPoint
                %is not close enough to zero
                PlatPoint = max(PlatPointTemp, ZeroDerLoc);
            end
        end
        EndOfResponse = LastTrough1;
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