%Giant Fiber membrane potential linear integration model
%written for ModelDb by Catherine von Reyn
%Card Lab, Janelia Research Campus, Ashbun VA
%06/2017
%Model and fits are from von Reyn et al., 2017.
%Parameters were found as described in von Reyn et al., 2017.
%Full citation:
% Feature Integration Drives Probabilistic Behavior in the Drosophila Escape Response.
% von Reyn CR, Nern A, Williamson WR, Breads P, Wu M, Namiki S, Card GM.
% Neuron. 2017 Jun 21;94(6):1190-1204.e6. doi: 10.1016/j.neuron.2017.05.036.
% PMID:  28641115 

function []=GF_loom_model()

clear all
close all

%neural delays estimted from data in paper
e_delay_non_LC4=11.9;  %in ms
i_delay_non_LC4=37.5; %in ms

e_delay_LC4=10.8; %in ms
i_delay_LC4=24.8; %in ms

sampling_frequency=20000; %in Hz
step_size=1/sampling_frequency; %in seconds
ms_step=step_size*1000; %in ms

% looming (r/v) stimuli
rv=[10 20 40 80]; %in ms
deg2rad=pi()/180;
start_size_rv=10 * deg2rad; %in rad to calculate looms
stop_size_rv=90 * deg2rad; %in rad to calculate looms

%imput stimulus parameters for
%constant angular velocity stimuli
ang_vel=[100 500 1000 2000]; %in deg/sec
start_size_vel=0; %in deg
stop_size_vel=90; %in deg 

%Non-LC4 components
%excatatory non_LC4 component parameters and equation (angular size)
Ce_nLC4=[1.7 42 0.52];
fe_nLC4=@(Ce_nLC4,x) Ce_nLC4(1)*exp(-(log(x)-log(Ce_nLC4(2))).^2/(2*Ce_nLC4(3)^2)); %log transformed gaussian
%inhibitory non_LC4 component parameters and equation (angular size)
Ci_nLC4=[-0.53 0.59 66 -11];
fi_nLC4 = @(Ci_nLC4,x) Ci_nLC4(1) + Ci_nLC4(2) ./ (1 + exp(-(x-Ci_nLC4(3))/Ci_nLC4(4)));
%for combined non_LC4 components
sf = 2.27; %inhibitory scale factor (level of inhibiton needed to overcome excitation at 90 deg)

%LC4_components
%excitatory LC4 component parameters and equation (angular velocity)
Ce_LC4=[-0.68 2.6 2.8 0.48];
fe_LC4=@(Ce_LC4,x) (Ce_LC4(1)+(Ce_LC4(2)./(1 + exp(-Ce_LC4(3)*(x - Ce_LC4(4))))));
%inhibitory LC4 component parameters and equation (angular size)
Ci_LC4=[-0.52 26 7.8];
fi_LC4=@(Ci_LC4,x) Ci_LC4(1)*exp(-(x-Ci_LC4(2)).^2./(2*Ci_LC4(3)^2));

for i=1:length(rv)
    figure('Name',sprintf('Looming stimulus rv = %d ms',rv(i)));
    %calculate theta and theta dot for this rv
    [time, theta, theta_dot] = calculate_loom(rv(i),start_size_rv,stop_size_rv,ms_step);
    
    %calculate Vm from non_LC4 component
    Vm_non_LC4=fe_nLC4(Ce_nLC4,neural_delay(theta,e_delay_non_LC4,ms_step))+fi_nLC4(Ci_nLC4,neural_delay(theta,i_delay_non_LC4,ms_step))*sf;
    %calculate Vm from LC4 component
    Vm_LC4=fe_LC4(Ce_LC4,neural_delay(theta_dot,e_delay_LC4,ms_step))+fi_LC4(Ci_LC4,neural_delay(theta,i_delay_LC4,ms_step));
    %calculate linear combination
    Vm=Vm_non_LC4+Vm_LC4;
   
    %plot non_LC4, LC4, and linearly combined components
    a1=subplot(4,1,1);
    plot(time,Vm_non_LC4,'r')
    title('Non LC4 component')
    ylabel('Vm (mV)')
    a2=subplot(4,1,2);
    plot(time,Vm_LC4,'b')
    title('LC4 component')
    ylabel('Vm (mV)')
    a3=subplot(4,1,3);
    plot(time,Vm,'k')
    title('Linear integration of both components')
    ylabel('Vm (mV)')
    linkaxes([a1,a2,a3],'y')
    ylim([-1 4])
    a4=subplot(4,1,4);
    plot(time,theta,'k');
    ylabel('Theta (degrees)')
    xlabel('Time (ms)')
    linkaxes([a1,a2,a3,a4],'x')
    xlim([-1000 500])
end


for i=1:length(ang_vel)
    figure('Name',sprintf('Constant angular velocity stimulus %d deg/sec',ang_vel(i)));
    %calculate theta and theta dot for this ang_vel
    [time, theta, theta_dot] = calculate_ang_vel(ang_vel(i),start_size_vel,stop_size_vel,ms_step);
    
    %calculate Vm from non_LC4 component
    Vm_non_LC4=fe_nLC4(Ce_nLC4,neural_delay(theta,e_delay_non_LC4,ms_step))+fi_nLC4(Ci_nLC4,neural_delay(theta,i_delay_non_LC4,ms_step))*sf;
    %calculate Vm from LC4 component
    Vm_LC4=fe_LC4(Ce_LC4,neural_delay(theta_dot,e_delay_LC4,ms_step))+fi_LC4(Ci_LC4,neural_delay(theta,i_delay_LC4,ms_step));
    %calculate linear combination
    Vm=Vm_non_LC4+Vm_LC4;
   
    %plot non_LC4, LC4, and linearly combined components
    a1=subplot(4,1,1);
    plot(time,Vm_non_LC4,'r')
    title('Non LC4 component')
    ylabel('Vm (mV)')
    a2=subplot(4,1,2);
    plot(time,Vm_LC4,'b')
    title('LC4 component')
    ylabel('Vm (mV)')
    a3=subplot(4,1,3);
    plot(time,Vm,'k')
    title('Linear integration of both components')
    ylabel('Vm (mV)')
    linkaxes([a1,a2,a3],'y')
    ylim([-1 4])
    a4=subplot(4,1,4);
    plot(time,theta,'k');
    ylabel('Theta (degrees)')
    xlabel('Time (ms)')
    linkaxes([a1,a2,a3,a4],'x')
    xlim([-500 1000])
    
end
end


function [time, theta, theta_dot] = calculate_loom(rv,start_size_rv,stop_size_rv,time_step)
deg2rad=pi()/180;
start_time=-rv/(tan(start_size_rv/2));
stop_time=-rv/(tan(stop_size_rv/2));
expansion_time=fliplr(stop_time:-time_step:start_time);
expansion_size=2*atan(rv./-expansion_time);
% add one second start to the begining and one second hold to the end
pre_expansion_time=(expansion_time(1)-time_step-1000:time_step:expansion_time(1)-time_step);
pre_expansion_size=zeros(size(pre_expansion_time));
post_expansion_time=(expansion_time(end)+time_step:time_step:expansion_time(end)+1000+time_step);
post_expansion_size=stop_size_rv.*ones(size(post_expansion_time));
time=[pre_expansion_time expansion_time post_expansion_time];
theta_rad = [pre_expansion_size expansion_size post_expansion_size];
theta=theta_rad./deg2rad;
theta_dot=[0 diff(theta)./(time_step)];
end

function [time, theta, theta_dot] = calculate_ang_vel(ang_vel,start_size_vel,stop_size_vel,time_step)
ang_vel=ang_vel/1000; %deg/ms
start_time=start_size_vel/ang_vel;
stop_time=stop_size_vel/ang_vel;
expansion_time=fliplr(stop_time:-time_step:start_time);
expansion_size=ang_vel*expansion_time;
% add one second start to the begining and one second hold to the end
pre_expansion_time=(expansion_time(1)-time_step-1000:time_step:expansion_time(1)-time_step);
pre_expansion_size=zeros(size(pre_expansion_time));
post_expansion_time=(expansion_time(end)+time_step:time_step:expansion_time(end)+1000+time_step);
post_expansion_size=stop_size_vel.*ones(size(post_expansion_time));
time=[pre_expansion_time expansion_time post_expansion_time];
theta = [pre_expansion_size expansion_size post_expansion_size];
theta_dot=[0 diff(theta)./(time_step)];
end

function [delayed_x]=neural_delay(x,delay,ms_step)
%delay and step needs to be in ms
step_shift=round(delay/ms_step);
delayed_x=[zeros(1,step_shift) x(1:end-step_shift)];
end