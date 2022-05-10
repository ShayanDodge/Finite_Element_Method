%Finite Element Method
NI=50;%Iteration
%Nodes
ND=21;
%Global Coordinates
X=[0 0.2 0.4 0.6 0.8 1.0 0 0.2 0.4 0.6 0.8 0 0.2 0.4 0.6 0 0.2 0.4 0 0.2 0];
Y=[0 0 0 0 0 0 0.2 0.2 0.2 0.2 0.2 0.4 0.4 0.4 0.4 0.6 0.6 0.6 0.8 0.8 1.0];
%Fixed Nodes
NP=15;
NDP=[1 2 3 4 5 6 7 11  12 15 16 18 19 20 21];
VAL=[0 0 0 0 0 50 0 100 0 100 0 100 0 100 50];
%Element
NE=25;
NL=[1 2 7;2 8 7;2 3 8;3 9 8;3 4 9;4 10 9;4 5 10;5 11 10;5 6 11;7 8 12;8 13 12
    ;8 9 13;9 14 13;9 10 14;10 15 14;10 11 15;12 13 16;13 17 16;13 14 17;
    14 18 17;14 15 18;16 17 19;17 20 19;17 18 20;19 20 21];
E0=1.0E-9/(36*pi);
ER=ones(1,NE);
%Evaluate coefficient matrix for each element and assemble globally
C=zeros(ND,ND);
for I=1:NE
%Find local coordinates XL,YL for element I
XL=X(NL(I,:));
YL=Y(NL(I,:));
P=[YL(2)-YL(3),YL(3)-YL(1),YL(1)-YL(2)];
Q=[XL(3)-XL(2),XL(1)-XL(3),XL(2)-XL(1)];
Area=0.5*abs(P(2)*Q(3)-Q(2)*P(3));
CE=ER(I)*(P'*P+Q'*Q)/(4*Area);
J=1:3;
L=1:3;
IR=NL(I,J);
IC=NL(I,L);
C(IR,IC)=C(IR,IC)+CE(J,L);
end
LF=setdiff(1:ND,NDP);
V=zeros(1,ND);
V(NDP)=VAL;
NF=length(LF);
for N=1:NI
    for ii=1:NF
        rowC=C(LF(ii),:);
        rowC(LF(ii))=0;
        V(LF(ii))=-1/(C(LF(ii),LF(ii)))*rowC*V';
    end
    figure(1),stem(V),drawnow
end
disp([{'Node','X','Y','Potential'};num2cell([(1:ND)',X',Y',V'])]);
