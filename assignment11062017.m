% to simulate the symbol error 
% performance and the bit error rate performance of the 8PSK 
close all;
clear all;
clc;

num = 1e4; % the number of symbols
M = 8; % for MPSK
N_bpsym = log2(M); % for MPSK, kbits/symbol
flag_pic = 00; % set 0 2 4 6 8 10 to show the figure when EbN0 = 0 2 4 6 8 10 
s_w = randi([0 M-1],1,num); % generate the omega of the signal symbols, 0~M-1
s_s = exp(1j*s_w*2*pi/M); % generate the signal of symbols

t_s = linspace(1,num,num); % time of symbols
t_b = 1:1/N_bpsym:num+1;
t_b = t_b(1:num*N_bpsym); % time of bits
t_p(1,:) = ones(1,num)*10; % the pack line
t_p(2,:) = ones(1,num)*-10; % the pack line
y_lable_coe = rats(0:2/M:2*(M-1)/M); % the coefficient of y lable
y_lable_coe = reshape(y_lable_coe,14,M).';
y_lable = cell(1,M); % % the text of y lable
for i = 1:M
    y_lable(i) = {[strip(y_lable_coe(i,:)),'*pi']}; % strip: eliminate extra space
end

ref = 0:M-1;
map = bitxor(ref,floor(ref/2)); % map for Gray code

s_wg = map(s_w+1); % after gray code
s_bg = dec2bin(s_wg)-'0';
s_b =  reshape(s_bg.',1,N_bpsym*num);

EbN0_db = 0:2:10;
BER = zeros(1,length(EbN0_db));
SER = zeros(1,length(EbN0_db));
Eb = 1/N_bpsym; % Eb is the average energy per bit, Es=abs(s_s)=1=log2(M)*Eb,for PSK,Es=1
N0_sd = Eb./10.^(EbN0_db/10); % SNR is Eb/N0, N0_sd is noise power spectral density
for m = 1:1:length(EbN0_db)
    N0 = sqrt(N0_sd(m)/2)*(randn(1,length(s_s))+1j*randn(1,length(s_s))); % N0_sd/2 is two-sided power spectral density of the noise
    s_sn = s_s+N0; % symbol signal with noise
    
%     s_wn = round((angle(s_sn))/(2*pi/M)); % the angle of symbol signal with noise,-pi~pi/(2*pi/M))
%     s_w1 = mod(s_wn,M);
    s_wn = mod(angle(s_sn)/(2*pi/M),M); % the angle of symbol signal with noise,-pi~pi/(2*pi/M))
    s_w1 = mod(round(s_wn),M);
    SER(m) = length(find(s_w1-s_w))./length(s_w1); % SER of symbol Signal with AWGN
    
    s_w1g = map(s_w1+1); % after gray code
    s_b1g = dec2bin(s_w1g)-'0'; % decimal number to binary number
    s_b1 =  reshape(s_b1g.',1,N_bpsym*num);
    BER(m) = length(find(s_b1-s_b))./length(s_b1); % BER of bit Signal with AWGN
    
    if EbN0_db(m) == flag_pic
        figure;
        stem(t_s,s_w*2*pi/M,'*'); % Original Symbol-Signal
        hold on;
        stairs(t_b,s_b,'b'); % Bit Signal
        hold on;
        stem(t_s,t_p(1,:),'--y'); % pack line
        hold on;
        stem(t_s,t_p(2,:),'--y'); % pack line
        legend('Original Phase of Symbol Signal',...
            'Bit Signal');
        hold off;
        set(gca,'ytick', 0:2*pi/M:2*pi*(M-1)/M);
        set(gca,'yticklabel', y_lable);
        xlim([1 20]);
        ylim([-0.25 2*pi+0.25]);
        figure;
        stem(t_s,s_w*2*pi/M,'*'); % Original Symbol-Signal
        hold on;
        stairs(t_b,s_b,'b'); % Bit Signal
        hold on;
        stem(t_s,s_wn*2*pi/M,'p'); % Symbol Signal+noise
        hold on;
        stem(t_s,t_p(1,:),'--y'); % pack line
        hold on;
        stem(t_s,t_p(2,:),'--y'); % pack line
        legend('Original Phase of Symbol Signal',...
            'Bit Signal',...
            ['Phase of Symbol Signal with ',num2str(EbN0_db(m)),'db AWGN']);
        hold off;
        set(gca,'ytick', 0:2*pi/M:2*pi*(M-1)/M);
        set(gca,'yticklabel', y_lable);
        xlim([1 20]);
        ylim([-0.25 2*pi+0.25]);
        figure;
        stem(t_s,s_w1*2*pi/M,'ro'); % Symbol Signal with noise after quantization
        hold on;
        stairs(t_b,s_b1,'r'); % Bit Signal with noise after quantization
        hold on;
        stem(t_s,t_p(1,:),'--y'); % pack line
        hold on;
        stem(t_s,t_p(2,:),'--y'); % pack line
        legend('Phase of Symbol Signal with noise after quantization',...
            'Bit Signal with noise after quantization');
        hold off;
        set(gca,'ytick', 0:2*pi/M:2*pi*(M-1)/M);
        set(gca,'yticklabel', y_lable);
        xlim([1 20]);
        ylim([-0.25 2*pi+0.25]);
%         figure;
%         stem(t_s,s_w*2*pi/M,'*'); % Original Symbol Signal
%         hold on;
%         stem(t_s,s_w1*2*pi/M,'o'); % Symbol Signal with noise after quantization
%         hold on;
%         stem(t_s,t_p(1,:),'--y'); % pack line
%         hold on;
%         stem(t_s,t_p(2,:),'--y'); % pack line
%         legend('Original Phase of Symbol Signal',...
%             'Phase of Symbol Signal with noise after quantization');
%         hold off;
%         set(gca,'ytick', 0:2*pi/M:2*pi*(M-1)/M);
%         set(gca,'yticklabel', y_lable);
%         xlim([1 20]);
%         ylim([-0.25 2*pi+0.25]);
%         figure;
%         stairs(t_b,s_b,'b'); % Original Bit Signal
%         hold on;
%         stairs(t_b,s_b1,'r--'); % Bit Signal with noise after quantization
%         hold on;
%         stem(t_s,t_p(1,:),'--y'); % pack line
%         hold on;
%         stem(t_s,t_p(2,:),'--y'); % pack line
%         legend('Original Bit Signal',...
%             'Bit Signal with noise after quantization');
%         hold off;
%         xlim([1 20]);
%         ylim([-0.25 1.25]);

    end
end

EbN0_db = linspace(0,10,num);
EbN0 = 10.^(EbN0_db/10);
EsN0 = N_bpsym*EbN0;
BER_t_bpsk = qfunc(sqrt(2*EbN0)); % rho = -1,BPSK
SER_t = 2*qfunc(sqrt(2*EsN0)*sin(pi/M));
BER_t = SER_t/N_bpsym; % coherence detected MPSK
% [BER_t,SER_t] = berawgn(EbN0_db,'psk',M,'nondiff'); % No results for coherent detection of differentially encoded PSK with M > 4.
figure;
semilogy(EbN0_db, SER_t,'b');
hold on;
semilogy(EbN0_db, BER_t,'r--');
hold on;
semilogy(EbN0_db, BER_t_bpsk,'r-.');
for m = 1:1:m
    hold on;
    plot(2*(m-1),SER(m),'b*');
    hold on;
    plot(2*(m-1),BER(m),'ro');
end
grid on;
ylabel('P');
xlabel('E_b/N_0 (dB)');
legend(['theoretical SER in ',num2str(M),'PSK'],...
    ['theoretical BER in ',num2str(M),'PSK'],...
    'theoretical BER in BPSK',...
    ['simulated SER in ',num2str(M),'PSK'],...
    ['simulated BER in ',num2str(M),'PSK']);
hold off;
return






