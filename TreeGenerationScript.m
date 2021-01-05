clear
tic

% Set Starting Conditions
MaxPoints = 100;                    % Sets the number of points generated
MaxSegment = (MaxPoints*2) + 1;     % Total Number of segments
MaxPoints = MaxPoints+1;            % For Caluclating r
Aperf = pi()*0.05^2;                % Area To be filled
load TestDataPatient                % Loads Patient Data
Seg = TestDataPatient;              % Assigns data to Seg structure
Points = 200;                       % Random Points created for each iteration
BiPoints = 200;                     % Number of bifurcation Points Tested

u = 3.6*10^-3;                      % Blood Viscosity
p = ((1.33*10^4)-(7.98*10^3));      % Blood Pressure
Qtot = 8.33*10^-6;                  % Total Flow Rate to LAD
Q = Qtot/ (MaxPoints+1);            % Set the flow for each terminal end
K = 21;                             % Start after the initial 21 Segments

while K <= MaxSegment-2             % -2 as it is K + 2
    
    P = (K-1)/2 + 2;                   % P is number of points currently
    r = CircBound(Aperf,P,MaxPoints);  % Calculate r of current Circle
    Seg(:,1:4) = Seg(:,1:4)*r;         % Scale Segments
    
    dthresh = sqrt(pi()*(r^2)/P);      % Set Distance threshold
    dcrit = 0;                         % Reset Critical Distance
    N = 1;                             % Start Counter for Points
    RandPoint = r*2*rand(Points,2);    % Set random points
    
    % Continue to loop until a satisfactory point is found
    while dcrit < dthresh
        
        C = sqrt((RandPoint(N,1)-r)^2 + (RandPoint(N,2)-r)^2);% Check it lies within circle
        
        if C < r
            [Seg,dcrit] = pointcheck(Seg,RandPoint,K,N,dcrit);
        end
        N = N+1;
        
        if N >= Points              % If no random point passes
            N = 1;                  % Reset Count
            dthresh = dthresh*0.9;  % Lower dthresh
        end
    end
    
    % The new Point has been found and is now going to be connected
    
    % Loop through every possible segment connection
    % M is the connection currently being tested
    
    % Set a high Vmin to begin with
    % if no point is found Vmin will still be 1000 at the end
    Vmin = 1000;
    
    for M = 1:K
        
        Int = 0;     % Set it to not initially intersect
        
        % If it is connecting to the original Segments it may only
        % bifurcate within a restricted range
        
        if M <= 21
            minx = min([Seg(M,1);Seg(M,3)]);
            maxx = max([Seg(M,1);Seg(M,3)]);
            
            miny = min([Seg(M,2);Seg(M,4)]);
            maxy = max([Seg(M,2);Seg(M,4)]);
            
            % Murrays law is given a large leeway for these due to the
            % restriction in bifurcation points
            r1murrayscale = 1.5;
            r2murrayscale = 0.5;
        else
            
            % Regular Murrays Law leniency
             r1murrayscale = 1.2;
             r2murrayscale = 0.8;
            % Finds the area for the new bifurcation points
            minx = min([Seg(M,1);Seg(M,3);Seg(K+2,3)]);
            maxx = max([Seg(M,1);Seg(M,3);Seg(K+2,3)]);
            
            miny = min([Seg(M,2);Seg(M,4);Seg(K+2,4)]);
            maxy = max([Seg(M,2);Seg(M,4);Seg(K+2,4)]);
            
        end
        
        % Generate Random Bifurcation Points within allocated area
        RandBiPoint = rand(BiPoints,2);
        RandBiPoint(:,1) = minx + (maxx - minx)*RandBiPoint(:,1); 
        RandBiPoint(:,2) = miny + (maxy - miny)*RandBiPoint(:,2);
        
        % Loop through bifurcation points finding the optimal within
        % constraints
        
        for  F = 1:BiPoints
            
            % Check it lies within circle
            C = sqrt((RandBiPoint(F,1)-r)^2 + (RandBiPoint(F,2)-r)^2);
            
            if C < r
                
                % Calculate l and r to find total volume of connection
                l1 = sqrt((RandBiPoint(F,1)-Seg(M,1))^2 + (RandBiPoint(F,2)-Seg(M,2))^2);
                l2 = sqrt((RandBiPoint(F,1)-Seg(M,3))^2 + (RandBiPoint(F,2)-Seg(M,4))^2);
                l3 = sqrt((RandBiPoint(F,1)-Seg(K+2,3))^2 + (RandBiPoint(F,2)-Seg(K+2,4))^2);
                
                % Calculates using Poiseuille
                r1 = ((Seg(M,8)+1)*Q*l1*8*u/(p*pi()))^(1/4);
                r2 = (Seg(M,8)*Q*l2*8*u/(p*pi()))^(1/4);
                r3 = (Q*l3*8*u/(p*pi()))^(1/4);
                
                Vnew = pi()*((l1*r1^2)+(l2*r2^2)+(l3*r3^2));
                
                % If this is a smaller volume than previously
                % found to be optimal
                
                if Vnew <= Vmin
                    
                    % Find the angle of the bifurcation
                    a = sqrt((Seg(M,3)-Seg(K+2,3))^2 + (Seg(M,4)-Seg(K+2,4))^2);
                    b = sqrt((RandBiPoint(F,1)-Seg(M,3))^2 + (RandBiPoint(F,2)-Seg(M,4))^2);
                    c = sqrt((RandBiPoint(F,1)-Seg(K+2,3))^2 + (RandBiPoint(F,2)-Seg(K+2,4))^2);
                    Angle = acosd((b^2 + c^2 - a^2)/(2*b*c));
                    
                    if Angle <= 80
                        
                        % Find the angle between the bifuraction and parent
                        a = sqrt((Seg(M,1)-Seg(K+2,3))^2 + (Seg(M,2)-Seg(K+2,4))^2);
                        b = sqrt((RandBiPoint(F,1)-Seg(M,1))^2 + (RandBiPoint(F,2)-Seg(M,2))^2);
                        c = sqrt((RandBiPoint(F,1)-Seg(K+2,3))^2 + (RandBiPoint(F,2)-Seg(K+2,4))^2);
                        Angle = acosd((b^2 + c^2 - a^2)/(2*b*c));
                        
                        if Angle >= 130
                            
                            % Find the angle between the bifuraction and parent
                            a = sqrt((Seg(M,1)-Seg(M,3))^2 + (Seg(M,2)-Seg(M,4))^2);
                            b = sqrt((RandBiPoint(F,1)-Seg(M,1))^2 + (RandBiPoint(F,2)-Seg(M,2))^2);
                            c = sqrt((RandBiPoint(F,1)-Seg(M,3))^2 + (RandBiPoint(F,2)-Seg(M,4))^2);
                            Angle = acosd((b^2 + c^2 - a^2)/(2*b*c));
                            
                            if Angle >= 130
                                
                                % Check the new bifurcation falls under Murray's law
                                % and do not grow in radius
                                r1murray = (r2^3 + r3^3)^(1/3)*r1murrayscale;
                                r2murray = (r2^3 + r3^3)^(1/3)*r2murrayscale;
                                rdiff = abs(r1 - r1murray);
                                
                                if r1 <= r1murray && r1>= r2murray &&  r1 > r2 && r1 > r3
                                    
                                    % Checks the parent Segment and its bifurcation obey Murray's law
                                    ParSeg = Seg(M,5);
                                    
                                    % As long as it is not segment 1
                                    if ParSeg ~= 0
                                        
                                        % Finds Segment Data
                                        DauSeg1 = Seg(ParSeg,6);
                                        DauSeg2 = Seg(ParSeg,7);
                                        QDauSeg1 = Seg(DauSeg1,8);
                                        QDauSeg2 = Seg(DauSeg2,8);
                                        l1 = sqrt((Seg(ParSeg,1)-Seg(ParSeg,3))^2 + (Seg(ParSeg,2)-Seg(ParSeg,4))^2);
                                        l2 = sqrt((Seg(DauSeg1,1)-Seg(DauSeg1,3))^2 + (Seg(DauSeg1,2)-Seg(DauSeg1,4))^2);
                                        l3 = sqrt((Seg(DauSeg2,1)-Seg(DauSeg2,3))^2 + (Seg(DauSeg2,2)-Seg(DauSeg2,4))^2);
                                        
                                        % If it is the new segment being tested the flow increases
                                        
                                        if DauSeg1 == M
                                            QDauSeg1 = QDauSeg1+1;
                                            l2 = sqrt((Seg(DauSeg1,1)-RandBiPoint(F,1))^2 + (Seg(DauSeg1,2)-RandBiPoint(F,2))^2);
                                        end
                                        if DauSeg2 == M
                                            QDauSeg2 = QDauSeg2+1;
                                            l3 = sqrt((Seg(DauSeg2,1)-RandBiPoint(F,1))^2 + (Seg(DauSeg2,2)-RandBiPoint(F,2))^2);
                                        end
                                        
                                        % Calculate r and Murrays law
                                        r1 = ((Seg(ParSeg,8)+1)*Q*l1*8*u/(p*pi()))^(1/4);
                                        r2 = (QDauSeg1*Q*l2*8*u/(p*pi()))^(1/4);
                                        r3 = (QDauSeg2*Q*l3*8*u/(p*pi()))^(1/4);
                                        
                                        r1murray = (r2^3 + r3^3)^(1/3)*r1murrayscale;
                                        r2murray = (r2^3 + r3^3)^(1/3)*r2murrayscale;
                                        rdiff = abs(r1 - r1murray);
                                    end
                                    
                                    if  r1 <= r1murray && r1>= r2murray
                                        if r1 > r2 && r1> r3
                                            
                                            % Check the other Segment and its daughters for Murray's
                                            DauSeg1 = Seg(M,6);
                                            DauSeg2 = Seg(M,7);
                                            
                                            % Provded it is not terminal
                                            if DauSeg1 ~=0
                                                l2 = sqrt((Seg(DauSeg1,1)-Seg(DauSeg1,3))^2 + (Seg(DauSeg1,2)-Seg(DauSeg1,4))^2);
                                                l3 = sqrt((Seg(DauSeg2,1)-Seg(DauSeg2,3))^2 + (Seg(DauSeg2,2)-Seg(DauSeg2,4))^2);
                                                l1 = sqrt((RandBiPoint(F,1)-Seg(M,3))^2 + (RandBiPoint(F,2)-Seg(M,4))^2);
                                                
                                                r1 = ((Seg(M,8)+1)*Q*l1*8*u/(p*pi()))^(1/4);
                                                r2 = (Seg(DauSeg1,8)*Q*l2*8*u/(p*pi()))^(1/4);
                                                r3 = (Seg(DauSeg2,8)*Q*l3*8*u/(p*pi()))^(1/4);
                                                
                                                r1murray = (r2^3 + r3^3)^(1/3)*r1murrayscale;
                                                r2murray = (r2^3 + r3^3)^(1/3)*r2murrayscale;
                                                rdiff = abs(r1 - r1murray);
                                            end
                                            
                                            if  r1 <= r1murray && r1>= r2murray
                                                if r1 >r2 
                                                    
                                                    % Intersection Check
                                                    for Z = 1:K
                                                        
                                                        % Skip the intersection check for the point it connects to
                                                        % Checks if it intersects with any of the previous segments
                                                        
                                                        if M ~= Z && Int ==0
                                                            % This checks the new daughter segment
                                                            line1 = [Seg(Z,1), Seg(Z,2); Seg(Z,3), Seg(Z,4)];
                                                            line2 = [RandBiPoint(F,1), RandBiPoint(F,2); Seg(K+2,3), Seg(K+2,4)];
                                                            Int = LineIntersectDau(Int,line1,line2);
                                                            
                                                            % This checks the old daughter segments
                                                            line2(2,1)= Seg(M,3);
                                                            line2(2,2) = Seg(M,4);
                                                            Int = LineIntersectDau(Int,line1,line2);
                                                            
                                                            % This tests the parent segment
                                                            if Seg(Z,1) ~= Seg(M,1)
                                                                if Seg(Z,3) ~= Seg(M,1)
                                                                    Int = LineIntersectPar(Seg(Z,:),Int,RandBiPoint(F,:),Seg(M,1),Seg(M,2));
                                                                end
                                                            end
                                                        end
                                                    end
                                                    
                                                    % If it has passed the
                                                    % checks save it
                                                    if Int == 0
                                                        Vmin = Vnew;
                                                        NewBiffx = RandBiPoint(F,1);
                                                        NewBiffy = RandBiPoint(F,2);
                                                        L = M;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    Int =0;
     
    % This saves all the new data necessary
    if Vmin ~= 1000 && Int ==0
        % Sets the points for the new segments
        Seg(K+1,1) = NewBiffx;
        Seg(K+1,2) = NewBiffy;
        Seg(K+1,3) = Seg(L,3);
        Seg(K+1,4) = Seg(L,4);
        
        % Changes the parent segments end point
        Seg(L,3) = NewBiffx;
        Seg(L,4) = NewBiffy;
        
        % Changes the starting point of the segment 3,4 already saved
        Seg(K+2,1) = NewBiffx;
        Seg(K+2,2) = NewBiffy;
        
        % Assigns the parent to it
        Seg(K+1,5) = L;
        Seg(K+2,5) = L;
        
        % Assigns parents to old daughters
        Daughter1 = Seg(L,6);
        Daughter2 = Seg(L,7);
        
        % K + 1 always connecting segment so becomes new parent
        if Daughter1 >= 1
            Seg(Daughter1,5) = K+1;
            Seg(Daughter2,5) = K+1;
        end
        
        % Assign the daughter to the new connecting segment
        Seg(K+1,6) = Seg(L,6);
        Seg(K+1,7) = Seg(L,7);
        
        % Assign the new daughters to the parent
        Seg(L,6) = K+1;
        Seg(L,7) = K+2;
        
        % Assigns flow to the new point and the connection point
        Seg(K+2,8) = 1;
        Seg(K+1,8) = Seg(L,8);
        
        V = 1;
        
        % Recursively Increases flow up the tree to the root
         while Seg(L,5) >=1
             Seg(L,8) = Seg(L,8) + 1;
             L = Seg(L,5);
         end
                
        % Root's flow always increased by 1 and is not affected by previous
        Seg(1,8) = Seg(1,8) + 1;
        
        % Calculates the radius of the new segments
        l = sqrt((Seg(K,1)-Seg(K,3))^2 + (Seg(K,2)-Seg(K,4))^2);
        dp = p/25;
        Seg(K,9)= ((Seg(K,8)*Q*l*8*u/(p*pi()))^(1/4));
        
        l = sqrt((Seg(K+1,1)-Seg(K+1,3))^2 + (Seg(K+1,2)-Seg(K+1,4))^2);
        Seg(K+1,9)= ((Seg(K+1,8)*Q*l*8*u/(p*pi()))^(1/4));
        K = K + 2; 
    end
    
    Seg(:,1:4) = Seg(:,1:4)/r;       % Rescaled back to neutral between (0-2)
    NewBiffx= 0;
    NewBiffy = 0;
end

% For some reason it is one extra at the end
Seg(1,8) = Seg(1,8) - 1;

%Scale Tree to Final size
P = (K-1)/2 + 2;                   % P is number of points created
r = CircBound(Aperf,P,MaxPoints);  % Calculate r of current Circle
Seg(:,1:4) = Seg(:,1:4)*r;
rsize = r;

% Goes through each segment and determines its bifurcation level

for I = 1:K
    V = 0;
    T = I;
    while Seg(T,5) >= 1
        T = Seg(T,5);
        V = V+1;
    end
    Seg(I,11) = V;
end


MaxBif = max(Seg(:,11));

% Calculates r for every segment based on the length and the flow rate
% through it
VTot = 0;
for E = 1:MaxSegment
    l = sqrt((Seg(E,1)-Seg(E,3))^2 + (Seg(E,2)-Seg(E,4))^2);
    Seg(E,9) = (Seg(E,8)*Q*l*8*u/((p/MaxBif)*pi()))^(1/4);
    VTot = VTot+ (pi()*(l*Seg(E,9)^2));
end

% Scales the radius for the plot so that its line width is not x10^-5
Seg(:,10) = Seg(:,9)*10000;

% Plots visual representation of vessels 
figure(1)
title('Generated Coronary Vessel Tree')
for E = 1:MaxSegment
    x1 = [Seg(E,1);Seg(E,3)];
    x2 = [Seg(E,2);Seg(E,4)];
    r = Seg(E,10);
    plot(x1,x2,'black','Linewidth',r)
    hold on
end
% Scales the plot so that it is centred (can be removed)
xlim([0 0.1044])
ylim([0 0.1044])
h = circle(rsize,rsize,rsize);             % Plot circle boundary


% Plots the average radius against the number of bifuractions down the tree

figure(2)
title('Average Segment Diameter Along Vessel Tree')
xlabel('Bifurcation level')
ylabel('Average Segment diameter (mm)')

% Calculates the average radius for each level
for N = 0:MaxBif
    k = Seg(:,11);
    Location = find(k==N);
    Radius = 0;
    for K = 1:size(Location)
        Radius = Radius + Seg(Location(K),9);
    end
    Radius = Radius/K;
    RadiusPlot(N+1,1)= N;
    RadiusPlot(N+1,2) = Radius;
end

% Plots the graph and converts to diameter and milimeters

RadiusPlot(1,1) = 0;
RadiusPlot(1,2) = Seg(1,9);
RadiusPlot(:,2) = RadiusPlot(:,2)*2000;
plot(RadiusPlot(:,1),RadiusPlot(:,2));

VTot

toc;
Compuationtime = toc
