function im_new_RGB=visibility_enhancement(im1)

    im_new_RGB=zeros(size(im1));
    for xx=1:3
        im=im1(:,:,xx);
        siz_im=size(im);
        [count,binLocations]=imhist(im);

        [C,U,LUT,H]=FCM21(im,3,2);
        th=0.25
        x1=0;

        for i=1:255
            if U(i,2)>th
                x1=i;break
            end
        end

        for i=1:255
            if U(i,1)<th
                x2=i;break
            end
        end


        for i=1:255
            if U(i,3)>th
                x3=i;break
            end
        end


        for i=x1+1:255
            if U(i,2)<th
                x4=i;break
            end
        end
        x1=double(x1+min(min(im1(:,:,xx))));
        x2=double(x2+min(min(im1(:,:,xx))));
        x3=double(x3+min(min(im1(:,:,xx))));
        x4=double(x4+min(min(im1(:,:,xx))));
        size_clus=size(U);
        U1=zeros(256,3);
        xxxx(1:256)=1:256;
        U1(uint16(min(min(im(:,:,:))))+1:uint16(max(max(im(:,:,:))))+1,1)=U(:,1);
        U1(uint16(min(min(im(:,:,:))))+1:uint16(max(max(im(:,:,:))))+1,2)=U(:,2);
        U1(uint16(min(min(im(:,:,:))))+1:uint16(max(max(im(:,:,:))))+1,3)=U(:,3);

        for i=1:256
            if(i>C(1))&&(U1(i,1)>0)
                U1(i,1)=1+(1-U1(i,1));
            end
            if(i>C(2))&&(U1(i,2)>0)
                U1(i,2)=1+(1-U1(i,2));
            end
            if(i>C(3))&&(U1(i,3)>0)
                U1(i,3)=1+(1-U1(i,3));
            end
        end

        U2=zeros(256,3);
        xxxx(1:256)=1:256;
        im_new1=zeros(siz_im);
        im_new2=zeros(siz_im);
        im_new3=zeros(siz_im);

        for i=1:siz_im(1)
            for j=1:siz_im(2)
                im_new1(i,j)=im(i,j).*U1(im(i,j)+1,1);
                im_new2(i,j)=im(i,j).*U1(im(i,j)+1,2);
                im_new3(i,j)=im(i,j).*U1(im(i,j)+1,3);
            end
        end
        im_new=(((im_new1)/3)+(im_new2)/3+(im_new3)/3);
        im_new=uint8(im_new);
        im_new_RGB(:,:,xx)=im_new;
    end
    im_new_RGB=uint8(im_new_RGB);

end


function[C,U,LUT,H]=FCM21(im,c,q,opt)

    if nargin<2||isempty(c),c=2;end
    if nargin<3||isempty(q),q=2;end
    if nargin<4||isempty(opt),opt=true;end
    if nargin<1||isempty(im)
        error('error: less input arguments')
    end
    msg='change variable to specify centroids.';
    if~isnumeric(c)||~isvector(c)
        error(msg)
    end
    if numel(c)==1&&(~isnumeric(c)||round(c)~=c||c<2)
        error(msg)
    end
    if~isnumeric(q)||numel(q)~=1||q<1.1
        error('Error')
    end
    if~islogical(opt)||numel(opt)>1
        error('Error')
    end
    if isempty(strfind(class(im),'int'))
        error('Please specify the image in integer format like uint8, int16 etc.')
    end
    Imin=double(min(im(:)));
    Imax=double(max(im(:)));
    I=(Imin:Imax)';
    if numel(c)>1
        C=c;
        opt=true;
    else
        if opt
            dI=(Imax-Imin)/c;
            C=Imin+dI/2:dI:Imax;
        else
            [C,~,H]=FCM2(im,c);
        end
    end
    if opt
        H=hist(double(im(:)),I);
        H=H(:);
    end
    % clearim
    dC=Inf;
    while dC>1E-5
        C0=C;
        D=abs(bsxfun(@minus,I,C));
        D=D.^(2/(q-1))+eps;
        U=bsxfun(@times,D,sum(1./D,2));
        U=1./(U+eps);
        UH=bsxfun(@times,U.^q,H);
        C=sum(bsxfun(@times,UH,I),1)./sum(UH,1);
        C=sort(C,'ascend');
        dC=max(abs(C-C0));

    end
    [~,LUT]=max(U,[],2);
end

function[C,LUT,H]=FCM2(im,c)

    if nargin<2||isempty(c),c=2;end
    if nargin<1||isempty(im)
        error('error: less input arguments')
    end
    msg='change variable to specify centroids.';
    if~isnumeric(c)||~isvector(c)
        error(msg)
    end
    if numel(c)==1&&(~isnumeric(c)||round(c)~=c||c<2)
        error(msg)
    end
    if isempty(strfind(class(im),'int'))
        error('Please specify the image in integer format like uint8, int16 etc.')
    end
    Imin=double(min(im(:)));
    Imax=double(max(im(:)));
    I=(Imin:Imax)';
    H=hist(double(im(:)),I);
    H=H(:);
    if numel(c)>1
        C=c;
        c=numel(c);
    else
        dI=(Imax-Imin)/c;
        C=Imin+dI/2:dI:Imax;
    end
    IH=I.*H;
    dC=Inf;
    while dC>1E-8
        C0=C;
        D=abs(bsxfun(@minus,I,C));
        [~,LUT]=min(D,[],2);
        for j=1:c
            C(j)=sum(IH(LUT==j))/sum(H(LUT==j));
        end
        C=sort(C,'ascend');
        dC=max(abs(C-C0));

    end
end


