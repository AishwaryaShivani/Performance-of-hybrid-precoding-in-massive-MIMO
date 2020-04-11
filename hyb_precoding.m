clc
clear all
close all
warning('off')

% System Parameters 
SNR_dB_range=-10:5:60;  % SNR in dB
Num_of_users=4; % Number of users
Trans_ant=64; %Number of UPA TX antennas
Trans_ant_w=sqrt(Trans_ant); % width
Trans_ant_h=sqrt(Trans_ant); % hieght 
ind_TX_w=reshape(repmat([0:1:Trans_ant_w-1],Trans_ant_h,1),1,Trans_ant_w*Trans_ant_h);
ind_TX_h=repmat([0:1:Trans_ant_h-1],1,Trans_ant_w);

rev_ant=4; %Number of UPA RX antennas
rev_ant_w=sqrt(rev_ant); % width 
rev_ant_h=sqrt(rev_ant); % hieght
ind_RX_w=reshape(repmat([0:1:rev_ant_w-1],rev_ant_h,1),1,rev_ant_w*rev_ant_h);
ind_RX_h=repmat([0:1:rev_ant_h-1],1,rev_ant_w);

% Channel Paramters
Num_paths=50; % Number of channel paths

Analog_BS=zeros(1,length(SNR_dB_range));% Carrying the rate of analog-only beamsteering

Hyb_ZF=zeros(1,length(SNR_dB_range)); % Carrying the rate of hybrid ZF precoder

Hyb_MMSE=zeros(1,length(SNR_dB_range)); % Carrying the rate of the hybrid MMSE precoder

Full_ZF=zeros(1,length(SNR_dB_range)); % Carrying the rate of the fully digital hybrid ZF precoder

Full_MMSE=zeros(1,length(SNR_dB_range)); % Carrying the rate of the fully digital hybrid MMSE precoder

Num_of_iter=500; % Number of iterations
    

for i=1:1:Num_of_iter
    
    % Generating user channels 
    [H,a_TX,a_RX]=generate_channels(Num_of_users,Trans_ant_w,Trans_ant_h,rev_ant_w,rev_ant_h,Num_paths); 
   
    % H is a 3-dimensional matrix having with Num_users,RX_ant,TX_ant dimensions

        
    % Analog beamforming precoding
    BRF=zeros(Trans_ant,Num_of_users); % Base station RF precoders 
    MRF=zeros(rev_ant,Num_of_users); % Mobile Station RF precoders 
    
    for j=1:1:Num_of_users
        BRF(:,j)=a_TX(:,j);
        MRF(:,j)=a_RX(:,j);
    end      
    
    % Construction of effective channels
    for j=1:1:Num_of_users
        Channel=zeros(rev_ant,Trans_ant);
        Channel(:,:)= H(j,:,:);
        He(j,:)=MRF(:,j)'*Channel*BRF ;    
    end
    
    % Effective channel for fully digital precoding
    for j=1:1:Num_of_users
        Channel=zeros(rev_ant,Trans_ant);
        Channel(:,:)= H(j,:,:);
        He_full(j,:)=MRF(:,j)'*Channel;    
    end
    
    % Baseband zero-forcing precoding
    ZF_hyb=He'*(He*He')^(-1); 
    for j=1:1:Num_of_users % Normalization of the hybrid precoders
        ZF_hyb(:,j)=ZF_hyb(:,j)/sqrt((BRF*ZF_hyb(:,j))'*(BRF*ZF_hyb(:,j)));
    end
    
    % Fully-digital zero-forcing precoding
    ZF_full=He_full'*pinv(He_full*He_full');
    for j=1:1:Num_of_users % Normalization of the hybrid precoders
        ZF_full(:,j)=ZF_full(:,j)/sqrt((ZF_full(:,j))'*(ZF_full(:,j)));
    end 
              
    % Spectral efficiency calculations
    count=0;
    for SNR_dB=SNR_dB_range
        count=count+1;
        SNR=10^(.1*SNR_dB)/Num_of_users; % SNR value
        sigma=1/SNR;
        
        % MMSE baseband precoder
        
        MMSE_Hyb=inv(He'*He+Num_of_users*sigma*BRF'*BRF)*He';
        
        for j=1:1:Num_of_users % Normalization of the hybrid precoders
            MMSE_Hyb(:,j)=MMSE_Hyb(:,j)/sqrt((BRF*MMSE_Hyb(:,j))'*(BRF*MMSE_Hyb(:,j)));
        end
        
        MMSE_Full=inv(He_full'*He_full+Num_of_users*sigma*eye(Trans_ant))*He_full';
        for j=1:1:Num_of_users % Normalization of the hybrid precoders
            MMSE_Full(:,j)=MMSE_Full(:,j)/sqrt((MMSE_Full(:,j))'*(MMSE_Full(:,j)));
        end
                
        for j=1:1:Num_of_users
            Int_set=[]; % interference index
            for i=1:1:Num_of_users
                if(i~=j)
                    Int_set=[Int_set i]; 
                end
            end
            Channel=zeros(rev_ant,Trans_ant);
            Channel(:,:)= H(j,:,:);
            [U_channel S_channel V_channel]=svd(Channel);
            
           
            
%           Analog-only beamforming
            SINR_BS=(SNR*(abs(MRF(:,j)'*Channel*BRF(:,j)).^2))/(SNR*sum((abs(MRF(:,j)'*Channel*BRF(:,Int_set)).^2))+1);
            Analog_BS(count)=Analog_BS(count)+log2(1+SINR_BS)/(Num_of_users*Num_of_iter);            
        end
    
        % ZF Hybrid Precoding
        Hyb_ZF(count)=Hyb_ZF(count)+log2(det(eye(Num_of_users)+SNR*(He*(ZF_hyb*ZF_hyb')*He')))/(Num_of_users*Num_of_iter);
            
        % Hybrid Precoding MMSE
        Hyb_MMSE(count)=Hyb_MMSE(count)+log2(det(eye(Num_of_users)+SNR*(He*(MMSE_Hyb*MMSE_Hyb')*He')))/(Num_of_users*Num_of_iter);
        
        % ZF fully digital precoding
         Full_ZF(count)=Full_ZF(count)+log2(det(eye(Num_of_users)+SNR*(He_full*(ZF_full*ZF_full')*He_full')))/(Num_of_users*Num_of_iter);
        
        % MMSE fully digital precoding
           Full_MMSE(count)=Full_MMSE(count)+log2(det(eye(Num_of_users)+SNR*(He_full*(MMSE_Full*MMSE_Full')*He_full')))/(Num_of_users*Num_of_iter);
    end % End of SNR loop
end % End of ITER loop

% Plotting the spectral efficiency
    figure(1)
    plot(SNR_dB_range,Hyb_ZF,'-mo','LineWidth',1.5);
    hold on; 
    plot(SNR_dB_range,Full_ZF,'-rs','LineWidth',1.5);
    hold on; 
    plot(SNR_dB_range,Hyb_MMSE,'-k*','LineWidth',1.5);
    hold on;    
   plot(SNR_dB_range,Full_MMSE,'--c','LineWidth',2);
if Num_paths==1
  hold on; 
     plot(SNR_dB_range,Analog_BS,'-r^','LineWidth',1.5);
    legend('ZF Hybrid Precoding','ZF Fully-Digital Precoding','MMSE Hybrid Precoding','MMSE Fully-Digital Precoding','Analog-only Beamsteering');

else
    hold on;
    plot(SNR_dB_range,Analog_BS,'-r^','LineWidth',1.5);
    legend('ZF Hybrid Precoding','ZF Fully-Digital Precoding','MMSE Hybrid Precoding','MMSE Fully-Digital Precoding','Analog-only Beamsteering');
end
xlabel('SNR (dB)','FontSize',12);
ylabel('Spectral Efficiency (bps/ Hz)','FontSize',12);
grid;
