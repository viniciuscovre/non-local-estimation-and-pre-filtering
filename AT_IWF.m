function filtrada = AT_IWF(img_antes, pre)

jan = 9;
dominio = 'ansc';

sigma2v = 1;

%cria a máscara do tamanho [jan jan] para um filtro da média ('average')
mean_filt = fspecial('average', [1 jan]);
mf = pre;
mg = imfilter(img_antes, mean_filt); %media da imagem ruidosa

% M = max(img_antes(:));
% m = min(img_antes(:));
% %img_antes = ((img_antes-m)/(M-m))*255;
% 
% M_ = max(pre(:));
% m_ = min(pre(:));
% % pre = ((pre-m_)/(M_-m_))*255;
% pre = ((pre-m_)/(M_-m_))*(M-m)+m;
% 
% M = max(M,M_);
% m = min(m,m_);

[l,c]=size(img_antes);
pad = floor(jan/2);

acumulador = zeros(l,c);

i=0;
for j = -pad : pad
    deslocada = circshift(pre, [i j]);
    diferenca = (deslocada - mf).^2;
    acumulador = acumulador + diferenca;
end

vf = acumulador/((jan^2)-1); % n = jan*jan (fórmula da variância)

img_antes = padarray(img_antes, [pad pad], 'symmetric');
vf = padarray(vf, [pad pad], 'symmetric');
mf = padarray(mf, [pad pad], 'symmetric');

ro = 0.95; % roV = roH

% Cálculo dos pesos para Rgg
pesos_rgg = zeros(jan^2, jan^2);
Ai=0; %Acumulador de linhas
for I = 1 : jan
    for J = 1 : jan
        Ai = Ai + 1;
        Aj = 0;%Acumulador de colunas
        for i = 1 : jan
            for j = 1 : jan
                Aj = Aj+1;
                pesos_rgg(Ai,Aj) = ro^sqrt((I-i)^2 + (J-j)^2);
            end
        end
    end
end

pesos_rff = pesos_rgg(ceil(jan^2/2),:)'; %linha central de Rgg equivale aos pesos de Rff

for i = pad+1 : l+pad
    for j = pad+1 : c+pad
        
        patch = img_antes(i-pad : i+pad, j-pad : j+pad) - mf(i-pad : i+pad, j-pad : j+pad); % Define o patch
        
        % Cálculo de Rgg
        Rgg = vf(i,j)*pesos_rgg;
        diag_princ = vf(i,j) * ones(jan^2, 1) + sigma2v; % sigma2v = 1 no domínio de Anscombe
        Rgg = Rgg - diag(diag(Rgg)) + diag(diag_princ); % seta a diagonal principal
        
        % Cálculo de Rff
        Rff = vf(i,j) * pesos_rff;
        
        % Cálculo de 'a'
        a = Rgg\Rff;
        
        filtrada(i-pad, j-pad) = mf(i,j) + sum(patch(:) .* a(:));
    end
end
filtrada(vf == 0) = pre(vf == 0);
% filtrada = filtrada * (M-m)/255 + m;
filtrada = noise_transform(filtrada, [dominio '_inverse']);
end