%Este programa � um software livre; voc� pode redistribu�-lo e/ou 
%
%   modific�-lo dentro dos termos da Licen�a P�blica Geral GNU como 
%
%    publicada pela Funda��o do Software Livre (FSF); na vers�o 2 da 
%
%    Licen�a, ou (na sua opini�o) qualquer vers�o.
%
%
%
%   Este programa � distribu�do na esperan�a de que possa ser �til, 
%
%    mas SEM NENHUMA GARANTIA; sem uma garantia impl�cita de ADEQUA��O a qualquer
%
%    MERCADO ou APLICA��O EM PARTICULAR. Veja a
%
%    Licen�a P�blica Geral GNU para maiores detalhes.
%
%
%
%    Voc� deve ter recebido uma c�pia da Licen�a P�blica Geral GNU
%
%    junto com este programa, se n�o, escreva para a Funda��o do Software
%
%    Livre(FSF) Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


%ALGOR�TIMO PARA AN�LISE DE P�RTICO PLANO - PROCESSO DOS DESLOCAMENTOS *
% Aplicativo com implementa��o do Processo dos Deslocamentos para a
% resolu��o de P�rticos Planos.
% Universidade Federal de Alagoas - UFAL
% Centro de Tecnologia - CTEC
% Programa de P�s-Gradua��o em Engenharia Civil - PPGEC
% -------------------------------------------------------------------------
% ALUNO: Ro�sio Oliveira Santos.
% Disciplina: Mec�nica Computacional de Estruturas
% Professor: Jo�o Carlos Cordeiro Barbirato(jccb@lccv.ufal.br)

%Entradas de dados:
% Quantidades de n�s e elementos:
% Nnos=input('N�mero de n�s: ');
% Nelem=input('N�mero de elementos: ');
% 
% %Coordenadas dos n�s
% for i=1:Nnos
%     disp('N�');
%     disp(i);
%     Cx(i)=input('Coordenada x: ');
%     Cy(i)=input('Coordenada y: ');
% end

% %Condi��es de contorno:
% disp(' Condi��es de contorno ');
% disp('Para a vincula��o dos n�s, utilize a seguinte simbologia:');
% disp('C�digo: 1 para deslocamentos prescritos e 0 para deslocamentos desconhecidos.');
% for i=1:Nnos
%     disp('Para o N�');
%     disp(i);
%     Cod(3*i-2)=input('Dire��o x - Vincula��o: ');
%     Cod(3*i-1)=input('Dire��o y - Vincula��o: ');
%     Cod(3*i)=input('rota��o - Vincula��o: ');
%     disp('A��es aplicadas:');
%     Fn(3*i-2)=input('Dire��o horizontal: ');
%     Fn(3*i-1)=input('Dire��o vertical: ');
%     Fn(3*i)=input('Momento: ');
% end

 %Conectividade dos elementos:
% for i=1:Nelem
%     disp(' Conectividade dos elementos ');
%     disp('Para o elemento');
%     disp(i);
%     conect(i,1)=input('N� j: ');
%     conect(i,2)=input('N� k: ');
%     prop(i,1)=input('Rigidez axial: (EA) ');
%     prop(i,2)=input('Rigidez � flex�o: (EI) ');
%     qforcas(i,1)=input('qx: ');
%     qforcas(i,2)=input('qy:  ');
% end

%entrada de dados (Exemplo):
Nnos = 5;
Nelem = 4;
 
%coordenadas
Cx = [0.0 0.0 2.0 3.0 3.0];
Cy = [0.0 1.0 1.0 2.0 0.0];
Cod = [1 1 1 0 0 0 1 1 0 0 0 0 0 1 0];
Fn = [0 0 0 0 0 0 0 0 0 3 -2 1 0 0 0];
conect = [1 2; 2 3; 3 4; 4 5]; 
prop = [10^6 1; 10^6 1; 10^6 1; 10^6 1];
qforcas = [0 0; 0 -10; 0 0; -5 0];

% MONTAGEM DA MATRIZ DE RIGIDEZ GLOBAL

R=zeros(3*Nnos);

for i=1:Nelem
    %Captura de informa��es dos n�s do elemento:
    %conectividade do elemento.
    j=conect(i,1);
    k=conect(i,2);
    %propriedades do elemento.
    ea=prop(i,1);
    ei=prop(i,2);
    %Comprimento do elemento:
    xj=Cx(j);
    xk=Cx(k);
    yj=Cy(j);
    yk=Cy(k);
    L=sqrt((xk-xj)^2+(yk-yj)^2);
    %Matriz de rigidez do elemento(i)tomada na coordenada local �:   
    rei(:,:)=[ea/L    0           0       ea/L    0           0
              0   12*ei/L^3   -6*ei/L^2   0   12*ei/L^3   6*ei/L^2
              0   -6*ei/L^2      4*ei/L   0   -6*ei/L^2    -2*ei/L
              ea/L    0           0       ea/L    0           0
              0   12*ei/L^3   -6*ei/L^2   0   12*ei/L^3   6*ei/L^2
              0   6*ei/L^2      -2*ei/L   0    6*ei/L^2     4*ei/L];
    
    % Sendo a Matriz de incid�ncia cinem�tica do elemento (i) dada por:
    cs=(xk-xj)/L;
    sn=(yk-yj)/L;
    
    betai(:,:)=[-cs -sn   0   0   0   0
                -sn  cs   0   0   0   0
                 0    0  -1   0   0   0
                 0    0   0   cs  sn  0
                 0    0   0   sn -cs  0
                 0    0   0   0   0   1];
             
    % C�lculo da matriz de rigidez do elemento (i) nas coordenadas globais:
    rgi(:,:) = betai(:,:)'*rei(:,:)*betai(:,:);
    
     %Montagem do vetor posi��o na matriz global
    j1=3*j-2;
    j2=3*j-1;
    j3=3*j;
    k1=3*k-2;
    k2=3*k-1;
    k3=3*k;
    vposi=[j1 j2 j3 k1 k2 k3];
         
       for lin=1:6
        for col=1:6
            R(vposi(lin),vposi(col)) = R(vposi(lin),vposi(col)) + rgi(lin,col);
        end
    end
end      

%MONTAGEM DO VETOR DE A��ES GLOBAIS

%Iniciando com o vetor de for�as zerado.
for i=1:(3*Nnos)
    Fne(i)=0;
end

% Intensidades nas coordenadas globais, distribu�das ao longo da barra
for elem=1:Nelem
   
    %Conectividade do elemento:
    j=conect(elem,1);
    k=conect(elem,2);
    
    %Comprimento do elemento:
    xj=Cx(j);
    xk=Cx(k);
    yj=Cy(j);
    yk=Cy(k);
    L=sqrt((xk-xj)^2+(yk-yj)^2);
    
    %Cargas
    qx=qforcas(elem,1);
    qy=qforcas(elem,2);
    
    %Matriz de incid�ncia cinem�tica do elemento (i):
    cs=(xk-xj)/L;
    sn=(yk-yj)/L;
    
    betai(:,:)=[-cs -sn   0   0   0   0
                -sn  cs   0   0   0   0
                 0    0  -1   0   0   0
                 0    0   0   cs  sn  0
                 0    0   0   sn -cs  0
                 0    0   0   0   0   1];
               
    % Chegamos ent�o as Cargas nas coordenadas globais:
    qxe = qx*cs + qy*sn;
    qye = qx*sn - qy*cs;
    
    Pne_ei = [qxe*L/2 qye*L/2 -qye*L^2/12 -qxe*L/2 -qye*L/2 -qye*L^2/12]';
    
    %Cargas nas coordenadas globais
    Pne_gi = betai(:,:)'*Pne_ei;
    
    %A��es no vetor das for�as nodais equivalentes
    j1=3*j-2;
    j2=3*j-1;
    j3=3*j;
    k1=3*k-2;
    k2=3*k-1;
    k3=3*k;
    
    vposi=[j1 j2 j3 k1 k2 k3];
    
    for lin=1:6
        Fne(vposi(lin)) = Fne(vposi(lin)) + Pne_gi(lin);
    end
end

%Completando o vetor das a��es globais
F=Fn-Fne;

% Por fim, completando o vetor das a��es globais:
F=Fn-Fne;

%C�lculo dos deslocamentos globais
%Resolu��o do sistema de equa��es lineares
% Para essa resolu��o, utilizou-se o m�todo do penalti

for i=1:3*Nnos
    if Cod(i)==1
        R(i,i)=R(i,i)*10^6;
%         F(i)=F(i)*10^6;
    end
end

%C�lculo dos deslocamentos globais
D = R\F';

%C�lculo das a��es nas extremidades das barras:

for elem=1:Nelem
        %Conectividade do elemento:
        j=conect(elem,1);
        k=conect(elem,2);
       
        % Propriedades:
        ea=prop(elem,1);
        ei=prop(elem,2);

        %Comprimento do elemento:
        xj=Cx(j);
        xk=Cx(k);
        yj=Cy(j);
        yk=Cy(k);
        L=((xk-xj)^2+(yk-yj)^2)^(1/2);

        %Cargas:
        qx=qforcas(elem,1);
        qy=qforcas(elem,2);
        
        %  Matriz de rigidez do elemento (i):   
    rei(:,:)=[ea/L    0           0       ea/L    0           0
              0   12*ei/L^3   -6*ei/L^2   0   12*ei/L^3   6*ei/L^2
              0   -6*ei/L^2      4*ei/L   0   -6*ei/L^2    -2*ei/L
              ea/L    0           0       ea/L    0           0
              0   12*ei/L^3   -6*ei/L^2   0   12*ei/L^3   6*ei/L^2
              0   6*ei/L^2      -2*ei/L   0    6*ei/L^2     4*ei/L];

        % Matriz de incid�ncia cinem�tica do elemento (i):
        cs=(xk-xj)/L;
        sn=(yk-yj)/L;

        betai(:,:)=[-cs -sn   0   0   0   0
                    -sn  cs   0   0   0   0
                     0    0  -1   0   0   0
                     0    0   0   cs  sn  0
                     0    0   0   sn -cs  0
                     0    0   0   0   0   1];

                 %Cargas nas coordenadas do elemento
        qxe = qx*cs + qy*sn;
        qye = qx*sn - qy*cs;

        Pne_ei = [qxe*L/2 qye*L/2 -qye*L^2/12 -qxe*L/2 -qye*L/2 -qye*L^2/12]';

        %A��es no vetor das for�as nodais equivalentes
        j1=3*j-2;
        j2=3*j-1;
        j3=3*j;
        k1=3*k-2;
        k2=3*k-1;
        k3=3*k;

        vposi=[j1 j2 j3 k1 k2 k3];
        % A Incid�ncia dos deslocamentos globais nos locais ser� dada por:
        for i=1:6
            Dg(i)=D(vposi(i));
        end
        
        %Transforma��o para as coordenadas locais
        De=betai*Dg';
        
        %Esfor�os solicitantes no elemento
        P_ei = Pne_ei + rei*De;
        
        El(elem).De = De;
        El(elem).P_ei = P_ei;
end
    
 %C�LCULO DAS REA��ES DE APOIO:
 
 RA=zeros(3*Nnos,1);
    
    for elem=1:Nelem
        %Conectividade do elemento:
        j=conect(elem,1);
        k=conect(elem,2);
        
        % Propriedades do elemento:
        ea=prop(elem,1);
        ei=prop(elem,2);

        %Comprimento do elemento:
        xj=Cx(j);
        xk=Cx(k);
        yj=Cy(j);
        yk=Cy(k);
        L=((xk-xj)^2+(yk-yj)^2)^(1/2);

        %Cargas:
        qx=qforcas(elem,1);
        qy=qforcas(elem,2);
        
        
    %Matriz de rigidez do elemento (i):   
    rei(:,:)=[ea/L    0           0       ea/L    0           0
              0   12*ei/L^3   -6*ei/L^2   0   12*ei/L^3   6*ei/L^2
              0   -6*ei/L^2      4*ei/L   0   -6*ei/L^2    -2*ei/L
              ea/L    0           0       ea/L    0           0
              0   12*ei/L^3   -6*ei/L^2   0   12*ei/L^3   6*ei/L^2
              0   6*ei/L^2      -2*ei/L   0    6*ei/L^2     4*ei/L];

        %Matriz de incid�ncia cinem�tica do elemento (i):
        cs=(xk-xj)/L;
        sn=(yk-yj)/L;

        betai(:,:)=[-cs -sn   0   0   0   0
                    -sn  cs   0   0   0   0
                     0    0  -1   0   0   0
                     0    0   0   cs  sn  0
                     0    0   0   sn -cs  0
                     0    0   0   0   0   1];

        %Cargas nas coordenadas do elemento:
        qxe = qx*cs + qy*sn;
        qye = qx*sn - qy*cs;

        Pne_ei = [qxe*L/2 qye*L/2 -qye*L^2/12 -qxe*L/2 -qye*L/2 -qye*L^2/12]';

        %A��es no vetor das for�as nodais equivalentes:
        j1=3*j-2;
        j2=3*j-1;
        j3=3*j;
        k1=3*k-2;
        k2=3*k-1;
        k3=3*k;

        vposi=[j1 j2 j3 k1 k2 k3];
        
        for i=1:6
           Pne_ei(i,1) = El(elem).P_ei(i);
        end
        
        %Cargas nas coordenadas globais:
        Pne_gi = betai*Pne_ei;

        %Rea��es nos apoios:
        for lin=1:6
            if Cod(vposi(lin))==1
                RA(vposi(lin)) = RA(vposi(lin)) + Pne_gi(lin);
            end
        end
    end
for i=1:3*Nnos
    if Cod(i)==1 && Fn(i)~=0
        RA(i) = RA(i)-Fn(i);
    end
end

% Os resultados ser�o mostrados em um arquivo no formato txt.

% SA�DA DE DADOS:

%ARQUIVO DE SA�DA DE DADOS
fp=fopen('resultado_portico.txt','wt');

fprintf(fp,'\n---                 UFAL/CTEC/PPGEC                ---');
fprintf(fp,'\n--- AN�LISE DE PORTICO PLANO - SA�DA DE DADOS   ---');
fprintf(fp,'\n---     MEC�NICA COMPUTACIONAL DE ESTRUTURAS    ---');
fprintf(fp,'\n---------------------------------------------------');
fprintf(fp,'\n---        Ro�sio Oliveira Santos      ---');


%Dados de Entrada
fprintf(fp,'*Dados de Entrada do problema:\n');
% informa��es dos N�s
fprintf(fp,'\nInforma��es dos N�s:\n');
fprintf(fp,'N�mero de n�s: %d\n',Nnos);
for i=1:Nnos
    fprintf(fp,'N�: %d\t(x,y)=(%f,%f)\tPresc: (%d,%d,%d)\tFor�a Ap: (%f,%f,%f)\n',i,Cx(i),Cy(i),Cod(3*i-2),Cod(3*i-1),Cod(3*i),Fn(3*i-2),Fn(3*i-1),Fn(3*i));        
end

%Elementos

fprintf(fp,'\nInforma��es dos Elementos:\n');
fprintf(fp,'N�meros de elemos: %d\n',Nelem);
for i=1:Nelem
    fprintf(fp,'Elem: %d\tN�J: %d\tN�K: %d\tEA: %f\tEI: %f\tqx: %f\tqy: %f\n',i,conect(i,1),conect(i,2),prop(i,1),prop(i,2),qforcas(i,1),qforcas(i,2));
end

%Dados de Sa�da  
fprintf(fp,'\n*Resultados obtidos:\n');
%Deslocamentos
fprintf(fp,'\nDeslocamentos de cada coordenada local dos elementos:\n');
for i=1:Nelem
    for j=1:6
        fprintf(fp,'%d\t%.3e\n',i,El(i).De(j));
    end
end

%Esfor��es internos solicitantes
fprintf(fp,'\nEsfor�os Internos Solicitantes:\n');
for i=1:Nelem
    for j=1:6
        fprintf(fp,'%d\t%.3e\n',i,El(i).P_ei(j));
    end
end

%Rea��es de apoio
fprintf(fp,'\nRea��es de Apoio:\n');
for i=1:Nnos
    for j=1:3
        fprintf(fp,'%d\t%.3e\n',i,RA(3*i+j-3));
    end
end

fclose(fp);
        
        
        
        
        
        
        
        
        
        
        
        

