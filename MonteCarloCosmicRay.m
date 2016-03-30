clear
t=[];
c=physconst('LightSpeed')/100000000;
h=100;
N=10000;
theta=2*pi*rand(N,1);
phi=random('norm',0.0000000,1.3/sqrt(2),[N,1]);
x=60*rand(N,1);
y=40*rand(N,1);
dx=diag(h*diag(abs(tan(phi)))*diag(cos(theta)));
dy=diag(h*diag(abs(tan(phi)))*diag(sin(theta)));
v=(2.998-0.2998)*rand(N,1)+0.2998*ones(N,1);
dt=-dx/c;
xf=x+dx;
yf=y+dy;
j=0;
for i=1:N
    if 0<=xf(i) & xf(i)<=60 & 0<=yf(i) & yf(i)<=40
        m{3*i}=[x(i) xf(i)];
        m{3*i+1}=[y(i) yf(i)];
        m{3*i+2}=[h 0];
        di=sqrt(dx(i)^2+dy(i)^2+h^2);
        r{i}=di;
        t=[t di/(10*v(i)) + dt(i)/10 + 24];
        
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
r=r(find(~cellfun(@isempty,r)));
% t=t(find(~cellfun(@isempty,t)));
% D(1:length(r))=r{1:length(r)};
% T(1:length(t))=t{1:length(t)};
figure(1)
plot3(n{1:length(n)})
figure(2)
hist(t,length(t))
% plot3(s{1:24},'*')
j/N
