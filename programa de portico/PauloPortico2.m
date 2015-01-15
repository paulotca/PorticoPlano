%Este programa é um software livre; você pode redistribuí-lo e/ou 
%
%   modificá-lo dentro dos termos da Licença Pública Geral GNU como 
%
%    publicada pela Fundação do Software Livre (FSF); na versão 2 da 
%
%    Licença, ou (na sua opinião) qualquer versão.
%
%
%
%   Este programa é distribuído na esperança de que possa ser útil, 
%
%    mas SEM NENHUMA GARANTIA; sem uma garantia implícita de ADEQUAÇÃO a qualquer
%
%    MERCADO ou APLICAÇÃO EM PARTICULAR. Veja a
%
%    Licença Pública Geral GNU para maiores detalhes.
%
%
%
%    Você deve ter recebido uma cópia da Licença Pública Geral GNU
%
%    junto com este programa, se não, escreva para a Fundação do Software
%
%    Livre(FSF) Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
clear;
clc;

%Algoritmo para análise de pórtico plano - Processo dos deslocamentos
%Paulo de Tarso Costa de Almeida
% %Quantidade de nós e elementos
% NumNos = input('Digite o número de nós da estrutura: ');
% NumEle = input('Digite o número de elementos: ');
% 
% %Capturando as coordenadas (x,y) dos nós
% fprintf('\n\nCoordenadas dos nós:\n');
% for i=1:NumNos
%     fprintf('\nNó %d:\n',i);
%     x(i) = input('  x = ');
%     y(i) = input('  y = ');
% end
% 
% fprintf('\n\nCondições de Contorno (Deslocamentos e Forças nodais)\n');
% fprintf('Obs.: +Para deslocamentos:\n');
% fprintf('     1 deslocamentos prescritos\n');
% fprintf('     0 deslocamentos desconhecidos\n');
% fprintf('      +Para força indicar o valor da força.');
% %Código: 1 deslocamentos prescritos e 0 deslocamentos desconhecidos
% for i=1:NumNos
%     fprintf('\nNó %d:\n',i);
%     CondContorno(3*i-2) = input('  x - Condição: ');
%     CondContorno(3*i-1) = input('  y - Condição: ');
%     CondContorno(3*i) = input('  rotação - Condição: ');
%     fprintf('\n  Força Nodal:\n');
%     ForLocais(3*i-2) = input('   Fx = ');
%     ForLocais(3*i-1) = input('   Fy = ');
%     ForLocais(3*i) = input('   Mz = ');
% end
% 
% fprintf('\n\nConecEleividade dos Elementos\n');
% 
% for i=1:NumEle
%     fprintf('\NumEleento %d:\n',i);
%     ConecEle(i,1) = input('  Nó j: ');
%     ConecEle(i,2) = input('  Nó k: ');
%     fprintf('\n  PropGeometricasriedades:\n');
%     PropGeometricas(i,1) = input('   EA = ');
%     PropGeometricas(i,2) = input('   EI = ');    
%     fprintf('\n  Forças Distribuídas:\n');
%     CarregaDistribu(i,1) = input('   qx = ');
%     CarregaDistribu(i,2) = input('   qy = ');    
% end

%entrada;

% %Exemplo1
 NumNos = 5;
 NumEle = 4;
% 
% %coordenadas
% x = [0.0 0.0 2.0 4.0 2.0];
 %y = [0.0 2.0 2.0 0.0 0.0];
% CondContorno = [1 1 1 1 1 0 0 0 0 1 1 1 1 1 1];
% ForLocais = [0 0 0 0 0 0 -5 -10 0 0 0 0 0 0 0];
 %ConecEle = [1 2; 2 3; 3 4; 3 5]; 
 %PropGeometricas = [10^6 1; 10^6 1; 10^6 1; 10^6 1];
 %CarregaDistribu = [0.5 0; 0 -1; 0 -1; -1 0.5];

%ExemploPortico
NumNos = 4;
NumEle =3;

%coordenadas
x = [0.0 0.0 8.0 12.0 ];
y = [0.0 4.0 4.0 0.0 ];
%Condições de Contorno (x1, y1, rot(z1))
CondContorno = [1 1 0 0 0 0 0 0 0 0 1 0];
ForLocais = [4 0 0 0 0 0 0 -2 10 0 0 0];
%ConecEleividade do elemento(nój nók)
ConecEle = [1 2; 2 3; 3 4]; 
%PropGeometricasriedades Geométricas (EA EI)
PropGeometricas = [10^6 1; 10^6 1; 10^6 1];
%Carregamento Distribuído (qx qy)
CarregaDistribu = [0 0; 0 -10; 0 0];

%Matriz de Rigidez Global
R = zeros(3*NumNos);
for i=1:NumEle
    %ConecEleividade do elemento i
    j = ConecEle(i,1);
    k = ConecEle(i,2);
    %PropGeometricasriedades do elemento i
    EA = PropGeometricas(i,1);
    EI = PropGeometricas(i,2);
    %comprimento do elemento i
    xj = x(j);
    xk = x(k);
    yj = y(j);
    yk = y(k);
    L  = sqrt((xk-xj)*(xk-xj)+(yk-yj)*(yk-yj));
    %Matriz de rigidez do elemento i
    r_ei = zeros(6);
    r_ei(1,1) = EA/L;
    r_ei(1,4) = EA/L;
    r_ei(2,2) = (12*EI)/(L^3);
    r_ei(2,3) = (-6*EI)/(L^2);
    r_ei(2,5) = (12*EI)/(L^3);
    r_ei(2,6) = (6*EI)/(L^2);
    r_ei(3,2) = (-6*EI)/(L^2);
    r_ei(3,3) = (4*EI)/(L);
    r_ei(3,5) = (-6*EI)/(L^2);
    r_ei(3,6) = (-2*EI)/(L);
    r_ei(4,1) = EA/L;
    r_ei(4,4) = EA/L;
    r_ei(5,2) = (12*EI)/(L^3);
    r_ei(5,3) = (-6*EI)/(L^2);
    r_ei(5,5) = (12*EI)/(L^3);
    r_ei(5,6) = (6*EI)/(L^2);
    r_ei(6,2) = (6*EI)/(L^2);
    r_ei(6,3) = (-2*EI)/(L);
    r_ei(6,5) = (6*EI)/(L^2);
    r_ei(6,6) = (4*EI)/(L);
    
    %matriz de incidência cinemática
    cs = (xk-xj)/L;
    sn = (yk-yj)/L;
    
    beta_ei = zeros(6);
    beta_ei(1,1) = -cs;
    beta_ei(1,2) = -sn;
    beta_ei(2,1) = -sn;
    beta_ei(2,2) = cs;
    beta_ei(3,3) = -1;
    beta_ei(4,4) = cs;
    beta_ei(4,5) = sn;
    beta_ei(5,4) = sn;
    beta_ei(5,5) = -cs;
    beta_ei(6,6) = 1;
    
    %Montagem de r_e nas coordenadas globais
    r_gi = beta_ei'*r_ei*beta_ei;
    
    %Montagem do vetor posição na matriz global
    J1 = 3*j - 2;
    J2 = 3*j - 1;
    J3 = 3*j;
    
    K1 = 3*k - 2;
    K2 = 3*k - 1;
    K3 = 3*k;    
    
    v_posi = [J1 J2 J3 K1 K2 K3];
    
    for lin=1:6
        for col=1:6
            R(v_posi(lin),v_posi(col)) = R(v_posi(lin),v_posi(col)) + r_gi(lin,col);
        end
    end
end

%Montagem do vetor de ações globais
%Para intensidades nas coordenadas globais, distribuídas ao longo da barra
for i=1:3*NumNos
    ForLocaise(i) = 0;
end
for i=1:NumEle
    %ConecEleividade do elemento i
    j = ConecEle(i,1);
    k = ConecEle(i,2);
    
    %comprimento do elemento i
    xj = x(j);
    xk = x(k);
    yj = y(j);
    yk = y(k);
    L  = sqrt((xk-xj)*(xk-xj)+(yk-yj)*(yk-yj));
    
    %cargas
    qx = CarregaDistribu(i,1);
    qy = CarregaDistribu(i,2);

    %matriz de incidência cinemática
    cs = (xk-xj)/L;
    sn = (yk-yj)/L;
    
    beta_ei = zeros(6);
    beta_ei(1,1) = -cs;
    beta_ei(1,2) = -sn;
    beta_ei(2,1) = -sn;
    beta_ei(2,2) = cs;
    beta_ei(3,3) = -1;
    beta_ei(4,4) = cs;
    beta_ei(4,5) = sn;
    beta_ei(5,4) = sn;
    beta_ei(5,5) = -cs;
    beta_ei(6,6) = 1;
    
    %cargas nas coordenadas do elemento
    qxe = qx*cs + qy*sn;
    qye = qx*sn - qy*cs;
    
    Pne_ei(1,1) = qxe*L/2;
    Pne_ei(2,1) = qye*L/2;
    Pne_ei(3,1) =-qye*(L*L)/12;
    Pne_ei(4,1) =-qxe*L/2;
    Pne_ei(5,1) =-qye*L/2;
    Pne_ei(6,1) =-qye*(L*L)/12;
    
    %cargas nas coordenadas globais
    Pne_gi = beta_ei'*Pne_ei;
    
    %Ações no vetor de forças nodais equivalentes
    J1 = 3*j - 2;
    J2 = 3*j - 1;
    J3 = 3*j;
    
    K1 = 3*k - 2;
    K2 = 3*k - 1;
    K3 = 3*k; 
    
    v_posi = [J1 J2 J3 K1 K2 K3];
    
    for lin=1:6
        ForLocaise(v_posi(lin)) = ForLocaise(v_posi(lin)) + Pne_gi(lin);
    end
end

%Completando o vetor das ações globais
F = ForLocais - ForLocaise;

%Cálculo dos deslocamentos globais
%Resolução do sistema de equações lineares
for i=1:3*NumNos
    if CondContorno(i) == 1
        R(i,i) = R(i,i)*(10^6);
        %F(i) = F(i)*(10^6);
    end
end

%Cálculo dos deslocamentos globais
D = R\F';

%Cálculo das ações nas extremidades das barras
for i=1:NumEle
    %ConecEleividade do elemento i
    j = ConecEle(i,1);
    k = ConecEle(i,2);
    
    %comprimento do elemento i
    xj = x(j);
    xk = x(k);
    yj = y(j);
    yk = y(k);
    L  = sqrt((xk-xj)*(xk-xj)+(yk-yj)*(yk-yj));
    
    %cargas
    qx = CarregaDistribu(i,1);
    qy = CarregaDistribu(i,2);

    %matriz de incidência cinemática
    cs = (xk-xj)/L;
    sn = (yk-yj)/L;
    
    r_ei = zeros(6);
    r_ei(1,1) = EA/L;
    r_ei(1,4) = EA/L;
    r_ei(2,2) = (12*EI)/(L^3);
    r_ei(2,3) = (-6*EI)/(L^2);
    r_ei(2,5) = (12*EI)/(L^3);
    r_ei(2,6) = (6*EI)/(L^2);
    r_ei(3,2) = (-6*EI)/(L^2);
    r_ei(3,3) = (4*EI)/(L);
    r_ei(3,5) = (-6*EI)/(L^2);
    r_ei(3,6) = (-2*EI)/(L);
    r_ei(4,1) = EA/L;
    r_ei(4,4) = EA/L;
    r_ei(5,2) = (12*EI)/(L^3);
    r_ei(5,3) = (-6*EI)/(L^2);
    r_ei(5,5) = (12*EI)/(L^3);
    r_ei(5,6) = (6*EI)/(L^2);
    r_ei(6,2) = (6*EI)/(L^2);
    r_ei(6,3) = (-2*EI)/(L);
    r_ei(6,5) = (6*EI)/(L^2);
    r_ei(6,6) = (4*EI)/(L);
    
    beta_ei = zeros(6);
    beta_ei(1,1) = -cs;
    beta_ei(1,2) = -sn;
    beta_ei(2,1) = -sn;
    beta_ei(2,2) = cs;
    beta_ei(3,3) = -1;
    beta_ei(4,4) = cs;
    beta_ei(4,5) = sn;
    beta_ei(5,4) = sn;
    beta_ei(5,5) = -cs;
    beta_ei(6,6) = 1;
    
    %cargas nas coordenadas do elemento
    qxe = qx*cs + qy*sn;
    qye = qx*sn - qy*cs;
    
    Pne_ei(1,1) = qxe*L/2;
    Pne_ei(2,1) = qye*L/2;
    Pne_ei(3,1) =-qye*(L*L)/12;
    Pne_ei(4,1) =-qxe*L/2;
    Pne_ei(5,1) =-qye*L/2;
    Pne_ei(6,1) =-qye*(L*L)/12;
    
    %Ações no vetor das forças nodais equivalentes
    J1 = 3*j - 2;
    J2 = 3*j - 1;
    J3 = 3*j;
    
    K1 = 3*k - 2;
    K2 = 3*k - 1;
    K3 = 3*k; 
    
    v_posi = [J1 J2 J3 K1 K2 K3];
    
    %incidência dos deslocamentos globais nos locais
    for l=1:6
        D_g(l,1) = D(v_posi(l));       
    end
    
    %transformação para as coordenadas locais
    D_e = beta_ei*D_g;
        
    %esforços solicitantes no elemento
    P_ei = Pne_ei + r_ei*D_e;
    
    elem(i).D_e = D_e;
    elem(i).P_ei = P_ei;
end

%Reações de apoio
Reacoes=zeros(3*NumNos,1);

for i=1:NumEle
        
    %ConecEleividade do elemento i
    j = ConecEle(i,1);
    k = ConecEle(i,2);
    
    %comprimento do elemento i
    xj = x(j);
    xk = x(k);
    yj = y(j);
    yk = y(k);
    L  = sqrt((xk-xj)*(xk-xj)+(yk-yj)*(yk-yj));
    
    %matriz de incidência cinemática
    cs = (xk-xj)/L;
    sn = (yk-yj)/L;
    
    beta_ei = zeros(6);
    beta_ei(1,1) = -cs;
    beta_ei(1,2) = -sn;
    beta_ei(2,1) = -sn;
    beta_ei(2,2) = cs;
    beta_ei(3,3) = -1;
    beta_ei(4,4) = cs;
    beta_ei(4,5) = sn;
    beta_ei(5,4) = sn;
    beta_ei(5,5) = -cs;
    beta_ei(6,6) = 1;
    
    for l=1:6
        Pne_ei(l,1) = elem(i).P_ei(l);
    end
    
    %cargas nas coordenadas globais
    Pne_gi = beta_ei*Pne_ei;
    
    %Ações no vetor de forças nodais equivalentes
    J1 = 3*j - 2;
    J2 = 3*j - 1;
    J3 = 3*j;
    
    K1 = 3*k - 2;
    K2 = 3*k - 1;
    K3 = 3*k; 
    
    v_posi = [J1 J2 J3 K1 K2 K3];
    
    for lin=1:6
        if CondContorno(v_posi(lin))==1
            Reacoes(v_posi(lin)) = Reacoes(v_posi(lin)) + Pne_gi(lin);
        end
    end
    
end
for i=1:3*NumNos
    if CondContorno(i)==1 && ForLocais(i)~=0
        Reacoes(i) = Reacoes(i)-ForLocais(i);
    end
end

%Arquivo de saída
fp=fopen('SaidaResultadoPorticos.txt','wt');


fprintf(fp,'\n---                UFAL/CTEC/PPGEC            ---');
fprintf(fp,'\n---       PORTICO PLANO RESULTADOS OBTIDOS    ---');
fprintf(fp,'\n---     MECÂNICA COMPUTACIONAL DE ESTRUTURAS  ---');
fprintf(fp,'\n-------------------------------------------------');
fprintf(fp,'\n---       Paulo de Tarso Costa de Almeida     ---');

%Matriz de rigidez
fprintf(fp,'Matriz de Rigidez\n');

for i=1:3*NumNos
    for j=1:3*NumNos
        fprintf(fp,'%.3e\t',R(i,j));
    end
    fprintf(fp,'\n');
end

%Vetor de forcas
fprintf(fp,'\nVetor de Forças F\n');

for i=1:NumNos
    for j=1:3
        fprintf(fp,'%d\t%.3e\n',i,F(3*i+j-3));
    end
end

%Deslocamentos
fprintf(fp,'\nDeslocamentos\n');
%fprintf(fp,'%d\n',3*NumNos);

for i=1:NumEle
    for j=1:6
        fprintf(fp,'%d\t%.3e\n',i,elem(i).D_e(j));
    end
end

%Forcas nos elementos
fprintf(fp,'\nEsforços Internos - Pe\n');
%fprintf(fp,'%d\n',3*NumNos);

for i=1:NumEle
    for j=1:6
        fprintf(fp,'%d\t%.3e\n',i,elem(i).P_ei(j));
    end
end

%Reacoes de apoio
fprintf(fp,'\nReações de apoio\n');
%fprintf(fp,'%d\n',3*NumNos);

for i=1:NumNos
    for j=1:3
        fprintf(fp,'%d\t%.3e\n',i,Reacoes(3*i+j-3));
    end
end

fclose(fp);
