% Script para a aplicação de uma filtragem com métodos contextuais e 
% estimação estatística por meios não-locais.
%
% ARTIGO WVC
%
% Abordagem: usar AT-NLM e AT-BM3D pare estimar a média e variância da 
% imagem. Sem sair do domínio de anscombe, aplicar filtros contextuais 
% (PWF, GWF, IWF e SWF) tendo como entrada as estimações estatísticas feita
% pelos métodos não-locais.

cd ~/Documents/Vinicius/dados_de_projecao/embrapa

l = dir('*3s.dat');

f1 = {'AT_NLM', 'AT_BM3D'}; % primeira pré-filtragem
f2 = {'PWF_1D', 'AT_IWF', 'AT_SWF', 'AT_GWF'}; % segunda pré-filtragem

PSNR = struct('name','', 'NoisyImage',.0, 'AT_NLM_D_PWF_1D',.0, 'AT_NLM_D_AT_IWF',.0,...
    'AT_NLM_D_AT_SWF',.0, 'AT_NLM_D_AT_GWF',.0, 'AT_BM3D_D_PWF_1D',.0,...
    'AT_BM3D_D_AT_IWF',.0, 'AT_BM3D_D_AT_SWF',.0, 'AT_BM3D_D_AT_GWF',.0);

SSIM = PSNR;

for k=1 : length(l)
    tomog = open_file_proj(l(k).name);
    tomog_ansc = noise_transform(tomog, 'ansc');
    reconstruida = retroprojecao(tomog);
    reconstruida = im2double(reconstruida);
    
    aux = strsplit(l(k).name, '.');
    PSNR(k).name = aux{1};
    SSIM(k).name = aux{1};
    
    img_ideal = open_file_proj(strrep(l(k).name,'_3s', '_20s'));
    %ao invés de abrir os dados _3s como sendo ideal, abra a _20s
    img_ideal = retroprojecao(img_ideal);
    
    PSNR(k).('NoisyImage') = psnr(reconstruida, img_ideal);
    SSIM(k).('NoisyImage') = ssim(reconstruida, img_ideal);
    
    cd ~/Documents/Vinicius/dupla_pre/imagens
    imwrite(img_ideal, [aux{1} '_ORIGINAL' '.png']);
    cd ~/Documents/Vinicius/dados_de_projecao/embrapa
    
    [lin,col] = size(reconstruida);
    
    for i=1 : length(f1)
        
        % Primeira pré-filtragem
        img_f1 = eval([f1{i} '(tomog_ansc)']);
        
        for j=1 : length(f2)
            
            % Segunda pré-filtragem
            img_f2 = eval([f2{j} '(tomog_ansc, img_f1)']);
            
            img_f2 = retroprojecao(img_f2);
            
            string = strcat(f1{i}, '_D_', f2{j});
            
            PSNR(k).(string) = psnr(img_ideal, img_f2);
            SSIM(k).(string) = ssim(img_ideal, img_f2);
            
            cd ~/Documents/Vinicius/dupla_pre/imagens
            imwrite(img_f2, [aux{1} '_' f1{i} '_D_' f2{j} '.png']);
            cd ~/Documents/Vinicius/dados_de_projecao/embrapa
            
            clear img_f2 string
        end
        
        clear img_f1
    end
    
    clear i j lin col aux reconstruida img_ideal tomog tomog_ansc
end

clear k f1 f2 l
