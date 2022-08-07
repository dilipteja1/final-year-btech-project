function file_name=optical2()
addpath('D:\m\codes\');
warning('on','MATLAB:RandStream:ActivatingLegacyGenerators') 
warning('on','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState')
clear all
NgType=1;
if NgType==1, nt='CP'; elseif NgType==2, nt='ZP'; end
Ch=0;
if Ch==0, chType='AWGN'; Target_neb=100; else chType='CH'; Target_neb=500; end
figure(Ch+1), clf
Nbps=2; M=2^Nbps; 
Nfft=128; 

%Ng=Nfft/4;
Ng=0;
Nsym=Nfft+Ng; % Symbol duration
EbN0=0:1:30; % EbN0
N_iter=1000; % Number of iterations for each EbN0
Nframe=1; % Number of symbols per frame
sigPow=0; 

file_name=['OFDM_BER_' chType '_' nt '_' 'GL' num2str(Ng) '.dat'];
fid=fopen(file_name, 'w+');

const=qammod([0:M-1],M);
scale=modnorm(const,'avpow',1);

for i=1:length(EbN0)
randn('state',0); rand('state',0); 
Neb=0; Ntb=0; 
for m=1:N_iter
% Tx
X= randi(M-1,1,((Nfft/2)-1)*Nframe);
Xmod= qammod(X,M,0)*scale;

if NgType~=2, x_GI=zeros(1,Nframe*Nsym);
elseif NgType==2, x_GI= zeros(1,Nframe*Nsym+Ng);
% Extend an OFDM symbol by Ng zeros
end
%kk1=1:Nfft/2; kk2=Nfft/2+1:Nfft; kk3=1:Nfft; kk4=1:Nsym;
kk1=1:(Nfft/2)-1;kk4=1:Nsym;
for k=1:Nframe
%X_shift= [Xmod(kk2) Xmod(kk1)];
X_shift=Xmod(kk1);
X_shift2=[0,X_shift,0,fliplr(conj(X_shift))];


x= ifft(X_shift2);
x_GI(kk4)= guard_interval(Ng,Nfft,NgType,x);
%kk1=kk1+Nfft; kk2= kk2+Nfft; kk3=kk3+Nfft; kk4=kk4+Nsym;
kk1=kk1+((Nfft/2)-1);kk4=kk4+Nsym;
end

bdc=12;
clip=sqrt((10.^(bdc/10))-1);
x1_GI=sum(x_GI.^2);
x11_GI=x1_GI/length(x_GI);
b=clip*sqrt(x11_GI);

x2_GI=x_GI+b;
for l=1:length(x2_GI)
if x2_GI(l)<0,
    x2_GI(l)=0;
end
end


y= x2_GI; % No channel
%if i==0 % Only to measure the signal power for adding AWGN noise
%y1=y(1:Nframe*Nsym); sigPow = sigPow + y1*y1'; continue;
%end
sigPow=sum(y.^2)/length(y);
% Add AWGN noise
snr = EbN0(i)+10*log10(Nbps); % SNR vs. Eb/N0 
noise_mag = sqrt((10.^(-snr/10))*sigPow);
y_GI = y + noise_mag*(randn(size(y)));
% Rx
kk1=(NgType==2)*Ng+[1:Nsym];
kk2=1:Nfft;
kk3=1:(Nfft/2)-1;
%kk4=Nfft/2+1:Nfft; kk5=[1:Nfft/2];
kk4=2:Nfft/2;
Y=zeros;
Xmod_r=zeros;
for k=1:Nframe
Y(kk2)= fft(remove_GI(Ng,Nsym,NgType,y_GI(kk1)));
Y_shift=Y(kk4);
Xmod_r(kk3) = Y_shift;
kk1=kk1+Nsym; kk2=kk2+Nfft; kk3=kk3+(Nfft/2)-1; kk4=kk4+Nfft;
%kk5=kk5+Nfft;
end
X_r=qamdemod(Xmod_r/scale,M,0);

Neb=Neb+sum(sum(de2bi(X_r,Nbps)~=de2bi(X,Nbps)));
%Ntb=Ntb+Nfft*Nframe*Nbps; %[Ber,Neb,Ntb]=ber(bit_Rx,bit,Nbps);
Ntb=Ntb+((Nfft/2)-1)*Nframe*Nbps;
if Neb>Target_neb, break; end
end
Ber = Neb/Ntb;
fprintf('EbN0=%3d[dB], BER=%4d/%8d =%11.3e\n', EbN0(i), Neb,Ntb,Ber)
fprintf(fid, '%d\t%11.3e\n', EbN0(i), Ber);
if Ber<1e-8, break; end
end


%if (fid~=0), fclose(fid); end
%plot_ber_op(file_name);


meanSquareValue = x*x'/length(x);
peakValue = max(x.*conj(x));
paprSymbol1 = peakValue/meanSquareValue; 


paprSymboldB1= 10*log10(paprSymbol1);
