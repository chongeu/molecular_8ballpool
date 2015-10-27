function EightBallPool

% Created by LEE CHONG EU (c) 2012
%
% Credit to Prof. CHEONG SIEW ANN and LEAW JIA NING
%
% Event Driven Molecular Dynamic Eight Ball Pool
% Run this filename in MATLAB to play the game
%
%
% ball diameter: 2.25 in
% 
% ball mass: 6 oz
% 
% ball mass moment of inertia: 2/5 mR2
% 
% ball-ball coefficient of friction (?): 0.03-0.08
% 
% ball-ball coefficient of restitution (e): 0.92-0.98
% 
% ball-cloth coefficient of rolling resistance (?): 0.005 - 0.015
% 
% ball-cloth coefficient of sliding friction (?): 0.15-0.4 (typical value: 0.2)
% 
% ball-cloth spin deceleration rate: 5-15 rad/sec2
% 
% ball-rail coefficient of restitution (e): 0.6-0.9
% 
% ball-table coefficient of restitution (e): 0.5
% 
% cue-tip-ball coefficient of friction (?): 0.6
% 
% cue-tip-ball coefficient of restitution (e): 0.71-0.75 (leather tip), 0.81-0.87 (phenolic tip)

clear all
clc

% load OfficialTable.mat
load haha.mat
% load DadaTable.mat
% load Hahah.mat
% load LShapeTable.mat
% % load PentagonTable.mat
% % load TriangleTable.mat
% % load SquareTable.mat

table.boundaryX = bx;
table.boundaryY = by;
table.innerBoundaryX = innerB(:,1);
table.innerBoundaryY = innerB(:,2);

table.right = max(table.boundaryX);
table.left = min(table.boundaryX);
table.upper = max(table.boundaryY);
table.lower = min(table.boundaryY);

table.pocketX = pocketX;
table.pocketY = pocketY;
table.pocketRadius = pocketRadius;

table.ballInitialX = ballInitialX;
table.ballInitialY = ballInitialY;

ball = initializeBall(table);

ball.isFirstTimePotted = 1;
ball.player = randperm(2);  % CREATE PLAYER INSTANCE
ball.currentPlayer = ball.player(1);    % PICK FIRST PLAYER
ball.nextPlayer = ball.player(2);       % SECOND PLAYER
ball.playerBall = [nan nan -ones(1,7) ones(1,7)];
ball.pottedBall = [];

disp(['Player ' num2str(ball.currentPlayer) ' is breaking. Good luck!'])


ball.plotrate = 0;    %s   % FRAME RATE
ball.cueBallv0 = 0.5;     %ms-1     % INITIAL CUE BALL VELOCITY
ball.frictionCoefficient = -0.2;  %ms-2     % BALL-CLOTH SLIDING FRICTION
ball.constant.epsilon = 1e-5;

ball.servoMode = 1;


drawTable(table,ball)
drawBalls(ball)

% SET CUE BALL POSITION
ball = setCueBallPosition(ball);


if ball.servoMode==1
    fhandle.ball = ball;
    fhandle.table = table;
    set(gcf,'UserData',fhandle);
    set(gcf,'WindowButtonMotionFcn',@getFutureCueBall);
else
    ball = getFutureCueBall(ball,table);
    updateSystem(ball,table);
end


function ball = setPlayerBall(ball,pottedBall)

solidIndex = find(ball.playerBall<0);
stripeIndex = find(ball.playerBall>0);

firstPottedBall = pottedBall(1);
if firstPottedBall>=3 && firstPottedBall<=9
    ball.currentPlayerColour = 'Solid/Blue';
    ball.nextPlayerColour = 'Stripe/Red';
    ball.playerBall(solidIndex) = ball.currentPlayer*ones(1,length(solidIndex));
    ball.playerBall(stripeIndex) = ball.nextPlayer*ones(1,length(stripeIndex));
    
elseif firstPottedBall>=10 && firstPottedBall<=16
    ball.currentPlayerColour = 'Stripe/Red';
    ball.nextPlayerColour = 'Solid/Blue';
    ball.playerBall(stripeIndex) = ball.currentPlayer*ones(1,length(stripeIndex));
    ball.playerBall(solidIndex) = ball.nextPlayer*ones(1,length(solidIndex));
    
end

disp(['Player ' num2str(ball.currentPlayer) ' is ' ball.currentPlayerColour])

ball.isFirstTimePotted = 0;


function ball = switchPlayer(ball)
buf = ball.currentPlayer;
ball.currentPlayer = ball.nextPlayer;
ball.nextPlayer = buf;

if ball.isFirstTimePotted==0
    buf1 = ball.currentPlayerColour;
    ball.currentPlayerColour = ball.nextPlayerColour;
    ball.nextPlayerColour = buf1;
end

disp(['Player ' num2str(ball.currentPlayer) ' now has the ball in hand.'])


function ball = setCueBallPosition(ball)
%% PLACE CUE BALL
ball.gameFlag = 0;
disp('Click cue ball position.')
ball = setPosition(ball,1,ginput(1));
drawBall(ball,1);


function getFutureCueBall(varargin)
fhandle = get(gcf,'UserData');
ball = fhandle.ball;
table = fhandle.table;

% function ball = getFutureCueBall(ball,table)
% GET VIRTUAL BALL POSITION


% SET CUE BALL VELOCITY
if ball.servoMode==1    
    MousePos = get(gca,'CurrentPoint');
else
    MousePos = ginput(1);
end

einput = MousePos(1,1:2);

ri0 = getPosition(ball,1);

e = einput - ri0;
e = e/sqrt(e*e'); 

vi0 = ball.cueBallv0*e;

ball = setVelocity(ball,1,vi0);

[b,delta_v,colT] = findFirstCollision(table,ball);

if isinf(colT)
    ball = updatePosition(ball);
    ball = updateVelocity(ball);

else
    ball = updatePosition(ball,colT);
    ball = updateVelocity(ball,colT);

    for i=1:length(b)
        ball.onTablev(b(i),:) = getVelocity(ball,b(i)) + delta_v(i,:);
    end
    
end

if ball.servoMode~=1
    ball = updatePottedBall(ball,table);
end

clf
drawTable(table,ball)
drawBalls(ball)

if ball.gameFlag==1
    disp('Opponent can has cue ball in hand.')
    ball = setCueBallPosition(ball);
    ball = getFutureCueBall(ball,table);
    
elseif ball.gameFlag==0
    drawCue(ri0,vi0,ball)
    DrawCircle(ri0(1),ri0(2),ball.radius,'-m')
    if numel(b)>1
        ri = getPosition(ball,1);
        vi = getVelocity(ball,1);
        drawVirtualBall(ball,ri0,ri,vi)
        rj = getPosition(ball,b(b>1));
        vj = getVelocity(ball,b(b>1));
        drawVirtualBall(ball,rj,rj,vj)
    else
        ri = getPosition(ball,1);
        vi = getVelocity(ball,1);
        drawVirtualBall(ball,ri0,ri,vi)
    end
end

if ball.servoMode==1
    set(gcf,'WindowButtonDownFcn',@getFutureCueBallServo);
end


function getFutureCueBallServo(varargin)
set(gcf,'WindowButtonMotionFcn','');
fhandle = get(gcf,'UserData');
ball = fhandle.ball;
table = fhandle.table;

MousePos = get(gca,'CurrentPoint');

einput = MousePos(1,1:2);

ri0 = getPosition(ball,1);

e = einput - ri0;
e = e/sqrt(e*e'); 

vi0 = ball.cueBallv0*e;

ball = setVelocity(ball,1,vi0);


[b,delta_v,colT] = findFirstCollision(table,ball);

if isinf(colT)
    ball = updatePosition(ball);
    ball = updateVelocity(ball);

else
    ball = updatePosition(ball,colT);
    ball = updateVelocity(ball,colT);

    for i=1:length(b)
        ball.onTablev(b(i),:) = getVelocity(ball,b(i)) + delta_v(i,:);
    end
    
end

ball = updatePottedBall(ball,table);

clf
drawTable(table,ball)
drawBalls(ball)

if ball.gameFlag==1
    disp('Opponent can has cue ball in hand.')
    ball = setCueBallPosition(ball);
    ball = getFutureCueBall(ball,table);
    
elseif ball.gameFlag==0
%     hold on
    drawCue(ri0,vi0,ball)
    DrawCircle(ri0(1),ri0(2),ball.radius,'-m')
    if numel(b)>1
        ri = getPosition(ball,1);
        vi = getVelocity(ball,1);
        drawVirtualBall(ball,ri0,ri,vi)
        rj = getPosition(ball,b(b>1));
        vj = getVelocity(ball,b(b>1));
        drawVirtualBall(ball,rj,rj,vj)
    else
        ri = getPosition(ball,1);
        vi = getVelocity(ball,1);
        drawVirtualBall(ball,ri0,ri,vi)
    end
end

updateSystem(ball,table);

    
function [b,delta_v,colT] = findFirstCollision(table,ball)
%% GET COLLISION TIME MATRIX
Nballs = size(ball.onTabler,1);

colTmat = nan(Nballs);
colTmatB = zeros(Nballs,1);
stopTmat = inf(Nballs,1);

for i=1:Nballs
    for j=i+1:Nballs
        colT = getCollisionTime(ball,i,j);

        colTmat(i,j) = colT;
        
    end
    [colTmatB(i),ball] = getCollisionTimeBoundary(table,ball,i);
    stopTmat(i) = getStopTime(ball,i);
end

colTmat = [colTmat colTmatB stopTmat];

[colT,bj] = min(min(colTmat));
[~,arraybi] = min(colTmat);
bi = arraybi(bj);


if bi<=Nballs && bj>=Nballs+2
    % BI IS STOPPING
%     disp([num2str(ball.surfaceNumber(bi)) '-stop'])
    delta_v = [0 0];
    b = bi;

elseif bi<=Nballs && bj==Nballs+1
    % BI IS HITTING A BOUNDARY
%     disp([num2str(ball.surfaceNumber(bi)) '-boundary'])
    delta_v = ball.delta_v{bi};
    b = bi;
    
else
    % BI AND BJ IS COLLIDING
%     disp([num2str(ball.surfaceNumber(bi)) '-' num2str(ball.surfaceNumber(bj))])
    [delta_vi,delta_vj] = getVelocityChange(ball,bi,bj,colT);
    b =[bi bj];
    delta_v = [delta_vi;delta_vj];
    
end


function [delta_vi,delta_vj] = getVelocityChange(ball,bi,bj,colT)
%% BALL-BALL VELOCITY CHANGE
rij = getRij(ball,bi,bj,colT);

ui = getVelocity(ball,bi);
uj = getVelocity(ball,bj);
ai = getAcceleration(ball,bi);
aj = getAcceleration(ball,bj);

ui = ui + ai*colT;
uj = uj + aj*colT;

uij = ui - uj;
bij = rij*uij';
delta_vi = -rij*bij/(4*ball.radius^2);
delta_vj = -delta_vi;


function rc = MakeCirclePoint2(a,b,x,y,r)
%% MAKE CIRCLE INTERSECTION
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


function colT = getCollisionTime(ball,i,j)
%% BALL-BALL COLLISION TIME

colT = inf;
% BALL-BALL COLLISION
rij = getPosition(ball,i) - getPosition(ball,j);    %getRij(ball,i,j,0);
uij = getVelocity(ball,i) - getVelocity(ball,j);    %getUij(ball,i,j,0);
aij = getAcceleration(ball,i)-getAcceleration(ball,j);

if sqrt(uij*uij')>ball.constant.epsilon && uij*rij'<0
    k = roots([.25*(aij*aij') uij*aij' rij*aij'+uij*uij' 2*rij*uij' rij*rij'-4*ball.radius*ball.radius]);
else
    k = [];
end

mink = min(k(k>ball.constant.epsilon));

if ~isempty(mink) && isreal(mink) && mink>ball.constant.epsilon
    colT = mink;
end


function [colT,ball] = getCollisionTimeBoundary(table,ball,i)
%% BALL-BOUNDARY COLLISION TIME

colT = inf;
ball.delta_v{i}=[0 0];

ri = getPosition(ball,i);
vi = getVelocity(ball,i);
ai = getAcceleration(ball,i);

stopT = getStopTime(ball,i);

ex = [ri(1) ri(1)+vi(1)*stopT+.5*ai(1)*stopT*stopT]';
ey = [ri(2) ri(2)+vi(2)*stopT+.5*ai(2)*stopT*stopT]';

[rjx, rjy] = polyxpoly(table.innerBoundaryX',table.innerBoundaryY',ex,ey, 'unique');

if ~isempty(rjx)
    
    % DISTANCE TO BOUNDARY
    s = sqrt((ri(1)-rjx).^2+(ri(2)-rjy).^2);

    % ELIMINATE REACHED NODE
    rjx(s<=ball.constant.epsilon)=[];
    rjy(s<=ball.constant.epsilon)=[];
    s(s<=ball.constant.epsilon)=[];
    
    % ELIMINATE THOSE THAT HAVE COLLIDED
    [s,indMin] = min(s);
    rjx = rjx(indMin);
    rjy = rjy(indMin);
    
    if numel(s)>0

        k = roots([.5*(-sqrt(ai*ai')) sqrt(vi*vi') -s]);

        mink = min(k(k>ball.constant.epsilon));
        
        if isreal(mink) %&& mink>ball.constant.epsilon

            colT = mink;
            
            % CALCULATE DELTA V_I
            rr = MakeCirclePoint2(0,2*pi,rjx,rjy,ball.radius+1e-3);
            [xj, yj] = polyxpoly(rr(:,1)', rr(:,2)',table.boundaryX,table.boundaryY, 'unique');
            
            xj = mean(xj);
            yj = mean(yj);
            
            rij = [rjx-xj rjy-yj];  % VECTOR NORMAL TO BOUNDARY
            vi = vi + ai*colT;
            rij = rij/sqrt(rij*rij');

            ball.delta_v{i} = -2*(vi*rij')*rij;

        end
    end
end


function t = getStopTime(ball,i)
%% GET STOP TIME
t = inf;
ui = getVelocity(ball,i);
ai = getAcceleration(ball,i);

if sqrt(ui*ui')>ball.constant.epsilon
    t = sqrt(ui*ui')/sqrt(ai*ai');
end


function ball = updateSystem(ball,table)
%% UPDATE SYSTEM
[b,delta_v,colT] = findFirstCollision(table,ball);

if isinf(colT)
    ball = updatePosition(ball);
    ball = updateVelocity(ball);

else
    ball = updatePosition(ball,colT);
    ball = updateVelocity(ball,colT);

    for i=1:length(b)
        ball.onTablev(b(i),:) = getVelocity(ball,b(i)) + delta_v(i,:);
    end
    
end

ball = updatePottedBall(ball,table);

pause(ball.plotrate)
drawTable(table,ball)
drawBalls(ball)

ball = gameEndEventHandler(ball,table,colT);


function ball = gameEndEventHandler(ball,table,colT)
%% GAME END EVENT HANDLER
if isinf(colT)
    if isempty(ball.pottedBall)
        disp('No ball potted.')
        ball = switchPlayer(ball);
        fhandle.ball = ball;
        fhandle.table = table;
        set(gcf,'UserData',fhandle);
        set(gcf,'WindowButtonMotionFcn',@getFutureCueBall);

    else

        cueBallPotted = ismember(1,ball.pottedBall);
        blackBallPotted = ismember(2,ball.pottedBall);
        firstPottedBall = ball.pottedBall(1);
        legalBallPotted = ball.playerBall(firstPottedBall)==ball.currentPlayer;

        % FIRST BALL OF PLAYER
        if ball.isFirstTimePotted==1 && ~cueBallPotted && ~blackBallPotted
            ball = setPlayerBall(ball,ball.pottedBall);
        end 
        
        ball.pottedBall(ball.pottedBall==1)=[];
        ball.pottedBall(ball.pottedBall==2)=[];
        
        % CLEARED THE POTTED BALL
        ball.playerBall(ball.pottedBall) = []; 
        ball.onTabler(ball.pottedBall,:) = [];
        ball.onTablev(ball.pottedBall,:) = [];
        ball.surfaceNumber(ball.pottedBall) = [];
        ball.surfaceColor(ball.pottedBall) = [];
        
        % RESET POTTED BASKET
        ball.pottedBall = [];
        set(gcf,'WindowButtonMotionFcn','');
        set(gcf,'WindowButtonDownFcn','');
        
        if blackBallPotted
            disp('Black potted.')
            if isempty(find(ball.playerBall==ball.currentPlayer))
                disp(['Player ' num2str(ball.currentPlayer) ' wins.'])
            else
                disp(['Player ' num2str(ball.nextPlayer) ' wins.'])
            end
            disp('Game over.')

        elseif cueBallPotted
            disp('Cue ball potted.')
            ball = switchPlayer(ball);
            ball = setCueBallPosition(ball);
            fhandle.ball = ball;
            fhandle.table = table;
            set(gcf,'UserData',fhandle);
            set(gcf,'WindowButtonMotionFcn',@getFutureCueBall);
            
        elseif ~legalBallPotted
            disp('Illegal ball potted.')
            ball = switchPlayer(ball);
            ball = setCueBallPosition(ball);
            fhandle.ball = ball;
            fhandle.table = table;
            set(gcf,'UserData',fhandle);
            set(gcf,'WindowButtonMotionFcn',@getFutureCueBall);
            
        else
            fhandle.ball = ball;
            fhandle.table = table;
            set(gcf,'UserData',fhandle);
            set(gcf,'WindowButtonMotionFcn',@getFutureCueBall);

        end
        
    end
else
    ball = updateSystem(ball,table);
end


function ball = updatePottedBall(ball,table)
%% UPDATE POTTED BALL

for i=1:size(ball.onTabler,1)
    ri = getPosition(ball,i);
    rpx = table.pocketX;
    rpy = table.pocketY;

    rip = sqrt((ri(1)-rpx).^2+(ri(2)-rpy).^2);

    if sum(rip<table.pocketRadius)>0
        if i>2
            disp([num2str(ball.surfaceNumber(i)) ' potted'])    
        end
        ball.pottedBall = [ball.pottedBall i];
        
        ball.onTabler(i,:) = [nan nan];
        ball.onTablev(i,:) = [0 0];
    end
end


function ball = updatePosition(ball,t)
%% UPDATE POSITION
Nballs = size(ball.onTabler,1);

for i=1:Nballs
    ri = getPosition(ball,i);
    ui = getVelocity(ball,i);
    ai = getAcceleration(ball,i);

    if nargin==1
        t = getStopTime(ball,i);
        if isinf(t)
            t = 0;
        end
    end

    ball.onTabler(i,:) = ri + ui*t + .5*ai*t*t;
end


function ball = updateVelocity(ball,t)
%% UPDATE VELOCITY

Nballs = size(ball.onTabler,1);
for i=1:Nballs
    ui = getVelocity(ball,i);
    ai = getAcceleration(ball,i);
    
    if nargin==1
        t = getStopTime(ball,i);
        if isinf(t)
            t = 0;
        end
    end
    vi = ui + ai*t;
    ball.onTablev(i,:) = vi;
% 	ball.onTablev(i,:) = round(vi*10000)/10000;
end


function rij = getRij(ball,i,j,t)
%% GET R_ij vector

ri = getPosition(ball,i);
ui = getVelocity(ball,i);
ai = getAcceleration(ball,i);

ri = ri + ui*t + .5*ai*t*t;

rj = getPosition(ball,j);
uj = getVelocity(ball,j);
aj = getAcceleration(ball,j);

rj = rj + uj*t + .5*aj*t*t;

rij = ri - rj;


function ball = setPosition(ball,i,ri)
%% SET POSITION
ball.onTabler(i,:) = ri;


function ri = getPosition(ball,i)
%% GET POSITION
ri = ball.onTabler(i,:);


function ball = setVelocity(ball,i,vi)
%% SET POSITION
ball.onTablev(i,:) = vi;


function vi = getVelocity(ball,i)
%% GET VELOCITY
vi = ball.onTablev(i,:);


function ai = getAcceleration(ball,i)
ui = getVelocity(ball,i);

if sqrt(ui*ui')>ball.constant.epsilon
    ai = ball.frictionCoefficient*ui;
else
    ai = [0 0];
end


function ball = initializeBall(table)
%% INITIALIZE BALL
setPosition = [0,0];

Nball = 15;

i = 0;
% CREATE THE TRIANGLE MESH
while Nball>0
    i = i+1;
    Nball = Nball-i;
    k = (1:i)';
    y = k-mean(k);
    x = (i-1)*sqrt(1-0.5^2)*ones(size(y));
    setPosition = [setPosition;[x y]];  % POSITION OF BALLS
end
   
ball.radius = 0.0123; %0.028575; % M
ball.mass = 0.163;  % KG

randomNoiseTol = 0.005;

% RESCALED INITIAL POSITION OF BALLS
% setPosition = (2*ball.radius+0.005*rand)*setPosition(2:16,:)+[ones(15,1)*(table.right-20*ball.radius) ones(15,1)*table.upper/2];
setPosition = (2*ball.radius + randomNoiseTol*rand)*setPosition(2:16,:)+[ones(15,1)*table.ballInitialX ones(15,1)*table.ballInitialY];

% BALL POSITION INDEX
stripePositionIndex = [1 3 4 8 10 11 13];
solidPositionIndex = [2 6 7 9 12 14 15];
blackPositionIndex = 5;

ball.surfaceNumber = [0 8 1 2 3 4 5 6 7 9 10 11 12 13 14 15];
ball.surfaceColor = {'-m';'-k';'-b';'-b';'-b';'-b';'-b';'-b';'-b';'-r';'-r';'-r';'-r';'-r';'-r';'-r'};

% % SHUFFLED BALL POSITION
ball.whiter = initilizeCueBall([0.5 1 0.5 1]);
ball.blackr = setPosition(blackPositionIndex,:);
ball.solidr = setPosition(solidPositionIndex(randperm(7)),:);
ball.striper = setPosition(stripePositionIndex(randperm(7)),:);

% % BALL INITIAL VELOCITY
ball.whitev = [0 0];
ball.blackv = [0 0];
ball.solidv = zeros(7,2);
ball.stripev = zeros(7,2);

% PUT BALL ON TABLE
ball.onTabler = [[nan nan]; ball.blackr; ball.solidr; ball.striper];
ball.onTablev = [ball.whitev; ball.blackv; ball.solidv; ball.stripev];


function drawBalls(ball)
%% DRAW 
hold on
for i=1:size(ball.onTabler,1)
    text(ball.onTabler(i,1),ball.onTabler(i,2),num2str(ball.surfaceNumber(i)),'HorizontalAlignment','center');
    DrawCircle(ball.onTabler(i,1),ball.onTabler(i,2),ball.radius,ball.surfaceColor{i});
end


function drawTable(table,ball)
%% DRAW TABLE
figure(1)
clf
hold on

if ball.isFirstTimePotted==0
    title(['Player ' num2str(ball.currentPlayer) '''s turn. (' ball.currentPlayerColour ')'],'fontsize',20)
end

% DRAW TABLE BED BOUNDARY
plot(table.boundaryX, table.boundaryY,'-b');

% DRAW INNER BOUNDARY LINE
% plot(table.innerBoundaryX,table.innerBoundaryY,'-g');

% DRAW TABLE BOUNDARY
% plot([table.left table.right table.right table.left table.left],[table.lower table.lower table.upper table.upper table.lower],'-k')

% DRAW POKET LINE
for i=1:length(table.pocketX)
    DrawCircle(table.pocketX(i),table.pocketY(i),table.pocketRadius(i),'-k');
end
set(gca,'yticklabel','','xticklabel','','ytick',[],'xtick',[]);
axis([table.left table.right table.lower table.upper]);
daspect([1 1 1]);


function drawBall(ball,i)
%% DRAW BALL
text(ball.onTabler(i,1),ball.onTabler(i,2),num2str(ball.surfaceNumber(i)));
DrawCircle(ball.onTabler(i,1),ball.onTabler(i,2),ball.radius,ball.surfaceColor{i});


function drawTrace(r1,r2,str)
%% DRAW TRACE
plot([r1(1) r2(1)],[r1(2) r2(2)],str);


function drawCue(r,v,ball)
normv = sqrt(v*v');
L = 0.4;
plot([r(1)+ball.radius*(-v(1)/normv) r(1)+L*(-v(1)/normv)],[r(2)+ball.radius*(-v(2)/normv) r(2)+L*(-v(2)/normv)],'-b','linewidth',2);


function drawVector(r,v,str)
lengthFactor = 3;   % REDUCE THIS FACTOR TO GET LONGER LENGTH
normv = sqrt(v*v')/lengthFactor;
%% DRAW VECTOR LINE
plot([r(1) r(1)+v(1)*normv],[r(2) r(2)+v(2)*normv],str);

% TO STATE THE VELOCITY
% text(r(1)+v(1)/(factor*normv),r(2)+v(2)/(factor*normv),num2str(sqrt(v*v')))


function drawVirtualBall(ball,ri0,ri,vi)
%% DRAW VIRTUAL BALL
hold on
DrawCircle(ri(1),ri(2),ball.radius,'-m')
drawTrace(ri0,ri,':m');
drawVector(ri,vi,'-k');


function r = initilizeCueBall(rr)
%% INITIALIZED CUE BALL
rx = rr(1) + (rr(2)-rr(1))*rand;
ry = rr(3) + (rr(4)-rr(3))*rand;
r = [rx ry];


function DrawCircle(x,y,r,str)
NOP=100;
p=linspace(0,2*pi,NOP);
q=ones(1,NOP)*r;
[A,B] = pol2cart(p,q);
A=A+x;
B=B+y;
plot(A,B,str);


