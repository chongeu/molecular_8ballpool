function CreateTable

% CREATED BY LEE CHONG EU (C) 2012
%
% This script create a pool table configuration to be play on
% EightBallPool.m
%

tableName = input('Key in table name (e.g. TriangleTable): ','s');

figure(1)
clf
set(gca,'xlim',[0 1],'ylim',[0 1])
aa = axis;
grid on
daspect([1 1 1])

disp('Draw the table (clockwise direction). Click outside axis when done.')
B = [];
while true
    figure(1)
    bufB = ginput(1);
    if bufB(1)<aa(1) || bufB(1)>aa(2) || bufB(2)<aa(3) || bufB(2)>aa(4)
        break
    else
        B = [B;bufB];
        figure(1); hold on
        plot(B(:,1),B(:,2),'-k');
        daspect([1 1 1])

    end
    
end
text(B(:,1),B(:,2),num2str((1:size(B,1))'));

bx=B(:,1)';
by=B(:,2)';

% clf
% figure(1); hold on
% plot(bx,by,'-k');
% daspect([1 1 1])
% text(bx,by,num2str((1:length(bx))'));

Npocket = input('Key in number of pocket(s): ');
pocketLine = zeros(Npocket,2);
pocketX = zeros(1,Npocket);
pocketY = zeros(1,Npocket);
pocketRadius = zeros(1,Npocket);
disp('Each hole is construct from any two points of table.');
disp('Input two edge number for each hole.');
for i=1:Npocket
    pocketLine(i,1) = input(['Pocket ' num2str(i) ' point 1: ']);
    pocketLine(i,2) = input(['Pocket ' num2str(i) ' point 2: ']);
    
    pocketX(i) = mean(bx(pocketLine(i,:)));
    pocketY(i) = mean(by(pocketLine(i,:)));

    pocketRadius(i) = sqrt(diff(bx(pocketLine(i,:)))^2+diff(by(pocketLine(i,:)))^2)/2;
end

% pocketLine
bx = fliplr(bx);
by = fliplr(by);

% CLOSE UP END TO END LINE
bx = [bx bx(2:3)];
by = [by by(2:3)];

radius = 0.0123;

innerB = [];

for i=1:numel(bx)-2
    r1 = [bx(i) by(i)];
    r2 = [bx(i+1) by(i+1)];
    r3 = [bx(i+2) by(i+2)];
    
    a = r2-r1;
    b = r3-r2;
    
    b = b/sqrt(b*b');
    a = a/sqrt(a*a');
    
    n = (b-a);
    n = n/sqrt(n*n');
    
    nsign = sign(a(1)*b(2)-a(2)*b(1));
    if nsign<0
        at = [a(2) -a(1)];
        bt = [b(2) -b(1)];
        
        % MAKE SURE THEY ARE POINTINT INWARD
        if at*n'>0
            at = -at;
        end
        
        if bt*n'>0
            bt = -bt;
        end
        
%         ra = r2 + radius*at;
%         rb = r2 + radius*bt;
        
        alpha = atan2(at(2),at(1));
        beta = atan2(bt(2),bt(1));
        rc = MakeCirclePoint(alpha,beta,r2(1),r2(2),radius);
        
        innerB = [innerB;rc];
    else
        
        r = r2 + n*nsign*radius*sqrt(2/(1+a*b'));
        innerB = [innerB; r];
    end
    
end
hold on
plot(innerB(:,1),innerB(:,2),'-g');
plot(bx,by,'-b');
for i=1:Npocket
    DrawCircle(pocketX(i),pocketY(i),pocketRadius(i),'-k');
end
daspect([1 1 1]);

disp('Click initial ball position')
bufIni=ginput(1);
ballInitialX = bufIni(1);
ballInitialY = bufIni(2);

save([tableName '.mat'],'innerB','bx','by','radius','pocketX','pocketY','pocketRadius','ballInitialX','ballInitialY')
disp([tableName '.mat saved!'])


function rc = MakeCirclePoint(a,b,x,y,r)
if a<0
    a=2*pi+a;
end
if b<0
    b=2*pi+b;
end

if abs(b-a)>pi
    b = b - 2*pi;
end

NOP=round(32*abs(a-b)/(2*pi));
p=linspace(a,b,NOP);
q=ones(1,NOP)*r;
[A,B] = pol2cart(p,q);
A=A+x;
B=B+y;
rc = [A' B'];

function rc = MakeCirclePoint2(a,b,x,y,r)
if a<0
    a=2*pi+a;
end
if b<0
    b=2*pi+b;
end

NOP=round(32*abs(a-b)/(2*pi));
p=linspace(a,b,NOP);
q=ones(1,NOP)*r;
[A,B] = pol2cart(p,q);
A=A+x;
B=B+y;
rc = [A' B'];

function DrawCircle(x,y,r,str)
NOP=100;
p=linspace(0,2*pi,NOP);
q=ones(1,NOP)*r;
[A,B] = pol2cart(p,q);
A=A+x;
B=B+y;
plot(A,B,str);