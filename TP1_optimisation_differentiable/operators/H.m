function blur = H(x,type,varargin)
%% H     op�rateur de flou
%     BLUR = H(X,TYPE,N)  floute X par un flou de type TYPE et de taille 
%     N impaire
%     . TYPE = 'gaussian','uniform'
%
%     Autrement dit, H(X,TYPE,N) calcule le produit matriciel :   
%           H*x 
%
%     Compatible :
%     . 1D si X est un signal,
%     . 2D si X est une image non vectoris�e.
% 
% 
%     Ex 1 (2D):
%       im = im2double(imread('cameraman.tif'));
%       
%       ker = 'gaussian';  
%       im_blur = H(im,ker);
% 
%       figure(1); clf;
%       subplot(121); imshow(im,[]);      title('image originale');
%       subplot(122); imshow(im_blur,[]); title('image flout�e');
%%

siz = size(x);

dim = (min(size(x)) > 1) + 1;

% check parameters
if dim == 1, siz = sort(siz,'descend'); end

if nargin ~= 3
    error('Nombre d''entr�es incorrect, la fonction attend 3 arguments.'); 
else 
    if max(strcmp(type,{'Arbitrary','arbitrary'})) == 1
        N = size(varargin{1});
    else
        N = varargin{1};
        if numel(N) > 1
            error('La taille du masque doit �tre un entier avec les options ''gaussian'' et ''uniform''.');
        end
        
        if dim == 1, N = [N,1]; end
        if dim == 2, N = [N,N]; end
    end
end

if min(mod(N,2)) == 0, error('Les dimensions du masque doivent �tre impaires.'); end

if sum(N > siz) > 0, error('Le masque est plus grand que le signal/l''image.'); end


switch type
    case 'gaussian', ker = fspecial('gaussian',N,N(1)/6);  
    case 'uniform',  ker = ones(N);
    case 'none',     ker = 1;
end



switch dim
    case 1
        if isrow(x), ker = ker'; end
        blur = real(ifft(psf2otf(ker,size(x)).*fft(x)));
        
    case 2
        blur = real(ifft2(psf2otf(ker,size(x)).*fft2(x)));
end

end