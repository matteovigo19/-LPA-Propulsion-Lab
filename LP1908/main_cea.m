CEA_output = CEA( ...
    'problem','o/f',3, ...
    'rocket','frozen','nfz',1, ...
    'p,bar',10, ...
    'reactants', ...
        'oxid','O2','wt%',100,'t,k',300, ...
        'fuel','HTPB','wt%',100,'t,k',300, ...
        'h,kj/mol',-58, ...
        'C',7.075,'H',10.650,'O',0.223,'N',0.063, ...
    'output', ...
        'massf','short','transport', ...
    'end','screen');

% A noi interessa:
% - T flame
% - Molar mass
% - specific heat ratio

Tf = CEA_output.output.froz.temperature(1);    % K
Mm = CEA_output.output.froz.mw(1);             % g/mol
cp = CEA_output.output.froz.cp_tran.froz(1);   % kJ/kg K

cp = cp * 1e3;    % J/kg K
R  = 8314.5 / Mm; % J/kg K
cv = cp - R;      % J/kg K
k  = cp / cv;     % -

%% Let's create a table, so we do not have to run the CEA all the time

p_vector   = [1, 1.01325, 2:1:100];     % bar
o_f_vector = [1:0.1:3, 3.5:0.5:20];

T_matrix = zeros(length(p_vector), length(o_f_vector)); % K
R_matrix = zeros(length(p_vector), length(o_f_vector)); % J/kg K
k_matrix = zeros(length(p_vector), length(o_f_vector)); % -

for i = 1:length(p_vector)
    for j = 1:length(o_f_vector)
        
        p   = p_vector(i);
        o_f = o_f_vector(j);
        
        CEA_output = CEA( ...
            'problem','o/f',o_f, ...
            'rocket','frozen','nfz',1, ...
            'p,bar',p, ...
            'reactants', ...
                'oxid','O2','wt%',100,'t,k',300, ...
                'fuel','HTPB','wt%',100,'t,k',300, ...
                'h,kj/mol',-58, ...
                'C',7.075,'H',10.650,'O',0.223,'N',0.063, ...
            'output', ...
                'massf','short','transport', ...
            'end', 'screen');
        % Estrazione proprietà
        Tf = CEA_output.output.froz.temperature(1);    % K
        Mm = CEA_output.output.froz.mw(1);             % g/mol
        cp = CEA_output.output.froz.cp_tran.froz(1);   % kJ/kg K
        
        cp = cp * 1e3;       % J/kg K
        R  = 8314.5 / Mm;    % J/kg K
        cv = cp - R;         % J/kg K
        k  = cp / cv;        % -
        
        % Salvataggio matrici
        T_matrix(i,j) = Tf;
        R_matrix(i,j) = R;
        k_matrix(i,j) = k;
        
    end
end

%% Plot
figure;
subplot(1,3,1)
surf(o_f_vector, p_vector, T_matrix);
ylabel("pressione, p, bar")
xlabel("oxider to fuel ratio, o/f, -")
zlabel("temperatura, T, k")
subplot(1,3,2)
surf(o_f_vector, p_vector, R_matrix);
ylabel("pressione, p, bar")
xlabel("oxider to fuel ratio, o/f, -")
zlabel("costante del gas, R, J/kg K")
subplot(1,3,3)
surf(o_f_vector, p_vector, k_matrix);
ylabel("pressione, p, bar")
xlabel("oxider to fuel ratio, o/f, -")
zlabel("specific heat ratio, k, -")

savedata(o_f_vector, p_vector, T_matrix,'CEA_T_of_p.txt');
savedata(o_f_vector, p_vector, R_matrix,'CEA_R_of_p.txt');
savedata(o_f_vector, p_vector, k_matrix,'CEA_k_of_p.txt');
 
%% da dati discreti → funzione continua
close all;clear;clc;

temp = readmatrix("CEA_T_of_p.txt");
o_f_vector = temp(1,2:end);
p_vector = temp(2:end,1);
T_matrix = temp(2:end,2:end);
temp = readmatrix("CEA_R_of_p.txt");
R_matrix = temp(2:end,2:end);
temp = readmatrix("CEA_k_of_p.txt");
k_matrix = temp(2:end,2:end);

T_fun_of_p = griddedInterpolant({o_f_vector, p_vector},T_matrix',"linear","none");
R_fun_of_p = griddedInterpolant({o_f_vector, p_vector},R_matrix',"linear","none");
k_fun_of_p = griddedInterpolant({o_f_vector, p_vector},k_matrix',"linear","none");

%% derivate parziali
RT_matrix = R_matrix.*T_matrix;
[dRTdOF_matrix, dRTdp_matrix] = gradient(RT_matrix, o_f_vector, p_vector); %il bro scrive T_matrix ma non sono convinto, per me dovrebbe essere RT_matrix
% [dTdOF_matrix, dTdp_matrix] = gradient(T_matrix, o_f_vector, p_vector);
% [dRdOF_matrix, dRdp_matrix] = gradient(R_matrix, o_f_vector, p_vector);
% [dkdOF_matrix, dkdp_matrix] = gradient(k_matrix, o_f_vector, p_vector);

dRTdp_fun_OF_p = griddedInterpolant({o_f_vector, p_vector}, dRTdp_matrix',"linear","none");
dRTdOF_fun_OF_p = griddedInterpolant({o_f_vector, p_vector}, dRTdOF_matrix',"linear","none");
% dTdp_fun_of_p = griddedInterpolant({o_f_vector, p_vector}, dTdp_matrix',"linear","none");
% dTdOF_fun_of_p = griddedInterpolant({o_f_vector, p_vector}, dTdOF_matrix',"linear","none");
% dRdp_fun_of_p = griddedInterpolant({o_f_vector, p_vector}, dRdp_matrix',"linear","none");
% dRdOF_fun_of_p = griddedInterpolant({o_f_vector, p_vector}, dRdOF_matrix',"linear","none");
% dkdp_fun_of_p = griddedInterpolant({o_f_vector, p_vector}, dkdp_matrix',"linear","none");
% dkdOF_fun_of_p = griddedInterpolant({o_f_vector, p_vector}, dkdOF_matrix',"linear","none");

%% save in a .mat file for future use
save('CEA_functions.mat', ...
    'T_fun_of_p', 'R_fun_of_p', 'k_fun_of_p', ...
    'dRTdp_fun_OF_p','dRTdOF_fun_OF_p', ...
    '-mat');

%% ------------------------------------
function savedata(X,Y,Z, filename)
% 1. Create the header row (0 followed by your X vector)
header = [0, X];

% 2. Attach the Y vector as the first column to the Z matrix
body = [Y',Z];

% 3. Combine them
final_output = [header; body];

% 4. Save as a tab-delimited text file
writematrix(final_output, filename, 'Delimiter', 'tab');
end