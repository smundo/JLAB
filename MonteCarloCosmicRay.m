% clear
t=[];
tphys=[]; % tphys=di/vi
tcon=[];  % tcon=dti
D=[];     % D=di
c=physconst('LightSpeed')/100000000; 
h=103;    % Seperation between detectors.
N=10000;  % Number of trials.
theta=2*pi*rand(N,1);
tnormpd=truncate(makedist('norm',0,1.3/sqrt(2)),-pi/2,pi/2);
phi=random(tnormpd,[N,1]);
% phi=zeros(N,1);
x=60*rand(N,1);
y=40*rand(N,1);
dx=diag(h*diag(abs(tan(phi)))*diag(cos(theta)));
dy=diag(h*diag(abs(tan(phi)))*diag(sin(theta)));
v=(2.9980-0.2998)*rand(N,1)+0.29980*ones(N,1); % Distribution of muons momentums.
% v=2.998*ones(N,1);
dt=-dx*1.5/c; % Signals travel time inside detectors.
xf=x+dx;
yf=y+dy;
j=0;
for i=1:N
    if 0<=xf(i) & xf(i)<=60 & 0<=yf(i) & yf(i)<=40
%         m{3*i}=[x(i) xf(i)];
%         m{3*i+1}=[y(i) yf(i)];
%         m{3*i+2}=[h 0];
        di=sqrt(dx(i)^2+dy(i)^2+h^2);
        t=[t di/(10*v(i)) + dt(i)/10];
        tphys=[tphys di/(10*v(i))];
        tcon=[tcon dt(i)/10];
        D=[D di];
        j=j+1;
    end
end
% for i=1:N
%         m{3*i}=[x(i) xf(i)];
%         m{3*i+1}=[y(i) yf(i)];
%         m{3*i+2}=[h 0];
%         j=j+1;
%     
% end

n=m(find(~cellfun(@isempty,m)));
% r=r(find(~cellfun(@isempty,r)));
% t=t(find(~cellfun(@isempty,t)));
% D(1:length(r))=r{1:length(r)};
% T(1:length(t))=t{1:length(t)};
% figure(1)
% plot3(n{1:length(n)})
figure(1)
hist(t,length(t))
% plot3(s{1:24},'*')
j/N
