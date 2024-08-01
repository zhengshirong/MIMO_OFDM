clear 

Ngtype = 1;%前缀
if Ngtype==1
    nt = 'CP';
elseif Ngtype==2
    nt = 'ZP';
end

Ch = 1;%信道
if Ch==0
    chtype  = 'AWGN';
%     targetneb = 100;
else
    chtype = 'CH';
%     targetneb = 500;
end
% figure(Ch+1),clf

powerdB = [0 -8 -17 -21 -25];
power = 10.^(powerdB/10);
Ntap=length(power);
norms=[1 sqrt(2) 0 sqrt(10) 0 sqrt(42)];     % BPSK 4-QAM 16-QAM

Delay = [0 3 5 6 8];
Lch=Delay(end)+1;
Nbps=4;
M=2^Nbps;%调制
%%%%%%%
N_iter = 1e1;%迭代
Nframe = 3 %符号数
Nfft =64 ;%=载波
Ng = Nfft/4;
Nsym = Nfft+Ng;%符号
Nvc= 0;%Nfft/4;%虚拟载波
Nuse = Nfft-Nvc;
sigp=0;         % Signal power initialization
EbN0 = [0:5:20];
for i = EbN0
%     ber()
Neb=0;%error
Ntb=0;%total
for m = 1:N_iter
%TX
x = randi([0 M-1],1,Nuse*Nframe); %3*64 = 192
xmod = qammod(x,M,'gray');
if Ngtype==1                  %CP
    xGI = zeros(1,Nframe*Nsym);%3*80=240
elseif Ngtype ==2             %ZP
    xGI = zeros(1,Nframe*Nsym+Ng);%3*80+16=256
end

%Nuse Nfft(Ng) Nsym 
kk1 = [1:Nuse/2];% 1 32
kk2 = [Nuse/2+1 :Nuse]; %33 64 
kk3 = [1:Nfft];    %1 64 
kk4 = [1:Nsym];    % 1 80

for k =1:Nframe%符号
    if Nvc ~=0
        x_shift = [0 xmod(kk2) zeros(1,Nvc-1) xmod(kk1)];
    else
        x_shift = [xmod(kk2) xmod(kk1)];
    end
    Xft = ifft(x_shift);
    xGI(kk4) = guardinterval(Ng,Nfft,Ngtype,Xft);%加保护间隔 %240/256
    kk1 = kk1+Nuse;
    kk2 = kk2+Nuse;
    kk3 = kk3+Nfft;
    kk4 = kk4+Nsym;
end

if Ch ==0 
    y = xGI;
else
    %多径信道
    channel = (randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(power/2);
    h = zeros(1,Lch);
    h(Delay+1)=channel;
    y = conv(xGI,h);
end
if i ==0
    y1 = y(1:Nframe*Nsym);
    sigp = y1*y1' ;
    continue;%下一次循环
end
%AWGN
snr = i +10*log10(Nbps*(Nuse/Nfft));
noisemag = sqrt( (10.^(-snr/10)) * sigp/2 );
yGI = y + noisemag*( randn(size(y))+j*randn(size(y))) ;

%RX
kk1 = (Ngtype==2)*Ng+[1:Nsym];
kk2 = [1:Nfft]; 
kk3 = [1:Nuse];   
kk4 = [Nuse/2+Nvc+1 : Nfft] ;   
kk5 = (Nvc~=0)+[1:Nuse/2];

if Ch ==1
    H = fft( [h zeros(1,Nfft-Lch) ] );%信道
    H_shift (kk3) = [H(kk4) H(kk5)] ;
end

for k  = 1:Nframe
    Y(kk2) = fft( removeGI(Ng,Nsym,Ngtype,yGI(kk1)) );
    Yshift = [Y(kk4) Y(kk5)];
    if Ch ==0
        xmodr(kk3) = Yshift;
    else
        xmodr(kk3) = Yshift./H_shift ;
    end
    kk1 = kk1+Nsym;
    kk2 = kk2+Nfft;
    kk3 = kk3+Nuse;
    kk4 = kk4+Nfft;
    kk5 = kk5+Nfft;
end
X_r = qamdemod(xmodr*norms(Nbps),M,'gray');
Neb = Neb+ sum( sum( de2bi( X_r,Nbps) ~= de2bi( x,Nbps)));
Ntb = Ntb +Nuse*Nframe*Nbps ;
if Neb >Targetneb,break ;end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%
function y =guardinterval(Ng,Nfft,Ngtype,ofdmsym)
if Ngtype ==1
    y = [ofdmsym(Nfft-Ng+1:Nfft) ofdmsym(1:Nfft)];
elseif Ngtype ==2
    y = [zeros(1,Ng) ofdmsym(1:Nfft)];
end
end

function y = removeGI(Ng,Nsym,Ngtype,ofdmsym)
    if Ngtype == 1          %CP
        y = ofdmsym(Ng+1:Nsym);
    elseif Ngtype ==2       %ZP
        y=ofdmsym(1:Nsym-Ng)+[ofdmsym(Nsym-Ng+1:Nsym) zeros(1,Nsym-2*Ng)];
    end
end



