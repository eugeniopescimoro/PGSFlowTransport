% Read inputs
fileID=fopen("/data/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp_4/TS1/LOGs/logTime",'r');
time = fscanf(fileID,'%f');
fclose(fileID);
fileID=fopen("/data/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp_4/TS1/LOGs/logFlux",'r');
conc = fscanf(fileID,'%f');
fclose(fileID);

time = time(1:100:end);
conc = conc(1:100:end);

u0 = 7.48813e-06;
D0 = 3.50486648940331e-05;
C = @(T) 0.5*erfc((2-u0*T)./(2*sqrt(D0*T)))+0.5*exp(u0*2/D0)*erfc((2+u0*T)./(2*sqrt(D0*T)))+0.5*(2+u0*2/D0+u0^2*T/D0).*exp(u0*2/D0).*erfc((2+u0*T)./(2*sqrt(D0*T)))-(u0^2*T/(pi*D0)).^0.5.*exp(u0*2/D0-(2+u0*T).^2./(4*D0*T));
f = chebfun(C, [min(time) max(time)]);
plot(f)
hold on
plot (time, conc, 'ro')

c = @(param, time) 0.5*erfc((2-param(1)*time)./(2*sqrt(param(2)*time)))+0.5*exp(param(1)*2/param(2))*erfc((2+param(1)*time)./(2*sqrt(param(2)*time)))+0.5*(2+param(1)*2/param(2)+param(1)^2*time/param(2)).*exp(param(1)*2/param(2)).*erfc((2+param(1)*time)./(2*sqrt(param(2)*time)))-(param(1)^2*time/(pi*param(2))).^0.5.*exp(param(1)*2/param(2)-(2+param(1)*time).^2./(4*param(2)*time));
param0 = [u0, D0];
[param,resnorm,~,exitflag,output] = lsqcurvefit(c,param0,time,conc)

hold on 
plot(time, c(param, time))
hold off

% t = chebfun(@(t) t, [0 max(conc)]);
% c = @(param, T) 0.5*erfc((2-param(1)*T)./(2*sqrt(param(2)*T)))+0.5*exp(param(1)*2/param(2))*erfc((2+param(1)*T)./(2*sqrt(param(2)*T)))+0.5*(2+param(1)*2/param(2)+param(1)^2*T/param(2)).*exp(param(1)*2/param(2)).*erfc((2+param(1)*T)./(2*sqrt(param(2)*T)))-(param(1)^2*T/(pi*param(2))).^0.5.*exp(param(1)*2/param(2)-(2+param(1)*T).^2./(4*param(2)*T));
% opt_par = fminsearch(@(param) immse(conc, c(param, time)), [1e-3,1e-3]);
% print(opt_par)