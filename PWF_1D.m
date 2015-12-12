function filtrada = PWF_1D(img_antes, pre)
% filtro de wiener generalizado (para projeções no domínio de anscombe)

jan = 3;

%pre = wiener(img_antes);
alfa = 1; %padrão fora do domínio de anscombe: 0.85

%img_antes = noise_transform(img_antes,'ansc');
%imgF = noise_transform(imgF,'ansc');

% M = max(img_antes(:));
% m = min(img_antes(:));
% img_antes = ((img_antes-m)/(M-m))*255;
%
% M_ = max(pre(:));
% m_ = min(pre(:));
% pre = ((pre-m_)/(M_-m_))*255;
%
% M = max(M,M_);
% m = min(m,m_);

[l,c]=size(img_antes);
d = floor(jan/2); %d = deslocamento

%cria a máscara do tamanho [jan jan] para um filtro da média ('average')
mean_filt = fspecial('average', [1 jan]);

mf = pre;
%mg = imfilter(img_antes, mean_filt); %media da imagem ruidosa
for k=1 : l
    mg(k,:) = medida_sinal(img_antes(k,:),'media',jan);
end
vr = 1; %variancia do ruido == 1 para trans ansc
ac_vf = zeros(l,c); %ac_vf = acumulador para a variância da imagem filtrada
ac_gKL = zeros(l,c);

i = 0;
for j = -d : d
    deslocadaF = circshift(pre, [i j]);
    dif_vf = (deslocadaF - mf).^2;
    ac_vf = ac_vf + dif_vf;
    
    deslocadaG = circshift(img_antes, [i j]);
    if (i==0 && j==0)
        continue;
    end
    dif_gKL = (deslocadaG - mg);
    ac_gKL = ac_gKL + dif_gKL;
end
% vf = ac_vf/((jan*jan)-1); % n = jan*jan (fórmula da variância)
clear i
for i=1 : l
    vf = medida_sinal(mf(i,:),'variancia',jan); % --> sigma^2_f (i,j)
    %Equação 5.13 da tese do Denis:
    filtrada(i,:) = mf(i,:) + (vf./(vf+vr)).*((alfa*(img_antes(i,:)-mf(i,:))) + ((1-alfa)*ac_gKL(i,:)));
    % filtrada = ((filtrada-min(filtrada(:)))/(max(filtrada(:))-min(filtrada(:))))*255;
    % % filtrada = filtrada * (M-m)/255 + m;
end
filtrada = noise_transform(filtrada,'ansc_inverse');
end