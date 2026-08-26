function ret = phase_jump(Vout, Vin,time,Theta)

Vin_n = Vin - Vout*Theta;
Vout_n = Vout + Vin*Theta;

index = 1;
while (time(index) < 0)
    index = index + 1;
end

FirstHalf = Vout_n(1:index);
SecondHalf = Vout_n(index:index+index);

% FirstHalf = Vout_n(1:index-3);
% SecondHalf = Vout_n(index+3:index+index);

ret = abs(mean(FirstHalf)-mean(SecondHalf));



