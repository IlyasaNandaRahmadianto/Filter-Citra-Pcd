citra=imread ('kuda.jpg'); //Import Gambar Dan Sesuaikan

R=rgb2gray(citra);    

filter = input('Pilih 1 Filter " batas / rata / median / konvolusi / lolos-rendah / lolos-tinggi / hboost / emboss " = ','s');
switch filter
    case 'batas'    %Ilyasa_07316
        
        [m,n] = size(R);
        G=zeros(m,n);
        f=double (R);
        
        for x =2 : m-1
            for y=2 : n-1

            minpiksel= min([f(x-1,y-1) f(x-1,y) f(x-1,y+1)  ...
            f(x,y-1) f(x,y+1)                               ...
            f(x+1,y-1) f(x+1,y) f(x+1,y+1)]);
            makspiksel=max ([f(x-1,y-1) f(x-1,y) f(x-1,y+1) ...
            f(x,y-1) f(x,y+1)                               ...
            f(x+1,y-1) f(x+1,y) f(x+1,y+1)]);
                if f(x,y) < minpiksel
                    G(x,y) = minpiksel;
                else
                    if f(x,y) > makspiksel
                        G(x,y) = makspiksel;
                    else
                        G(x,y) = f(x,y);
                    end
                end
            end
        end

        G=uint8(G);
        subplot(1,2,1); imshow(R); title('Citra Asli');
        subplot(1,2,2); imshow(G); title('Filter Batas');
        
        %Hitung MSE, RMSE, PNSR
        [baris,kolom]=size(R);
        selisih=zeros(baris,kolom);
        
        for i=1:baris
            for j=1:kolom
                selisih(i,j)=G(i,j)-R(i,j);
               

            end
        end
        
        %hitung kuadrat
        kuadrat=zeros(baris,kolom);
        for i=1:baris
            for j=1:kolom
                kuadrat(i,j)=selisih(i,j)*selisih(i,j);
                
            end
        end

        kuadrat;
        totalkuwadrat=0;
        
        for i=1:baris
            for j=1:kolom
                totalkuwadrat=totalkuwadrat+kuadrat(i,j);
                
            end
        end

        totalkuwadrat;
        
        MSE=totalkuwadrat/(baris*kolom);
        
        % Hitung RMSE
        RMSE = sqrt(MSE);
        
        %Hitung PSNR
        PSNR = 10*log((255^2)/MSE);
        
        fprintf('\nFilter Batas\n')
        fprintf('MSE  : %12.4f;\n',MSE)
        fprintf('RMSE : %12.4f;\n',RMSE)
        fprintf('PSNR : %12.4f;\n',PSNR)
        
    case 'rata'     %Ilyasa_07316
        
        [m,n] = size(R);
        G=zeros(m,n);
        f=double (R);
        
        for x =2 : m-1
            for y=2 : n-1

            jum=f(x-1,y-1)+ f(x-1,y) +f(x-1,y+1)+   ...
            f(x,y-1) +f(x,y)+f(x,y+1)+              ...
            f(x+1,y-1)+f(x+1,y)+f(x+1,y+1);
            G(x,y) = uint8(1/9 * jum);
            end
        end

            G=uint8(G);
            subplot(1,2,1); imshow(R); title('Citra Asli');
            subplot(1,2,2); imshow(G); title('Filter Rata-Rata');
            
            %Hitung MSE, RMSE, PNSR
            [baris,kolom]=size(R);
            selisih=zeros(baris,kolom);

            for i=1:baris
                for j=1:kolom
                    selisih(i,j)=G(i,j)-R(i,j);


                end
            end

            %hitung kuadrat
            kuadrat=zeros(baris,kolom);
            for i=1:baris
                for j=1:kolom
                    kuadrat(i,j)=selisih(i,j)*selisih(i,j);

                end
            end

            kuadrat;
            totalkuwadrat=0;

            for i=1:baris
                for j=1:kolom
                    totalkuwadrat=totalkuwadrat+kuadrat(i,j);

                end
            end

            totalkuwadrat;

            MSE=totalkuwadrat/(baris*kolom);

            % Hitung RMSE
            RMSE = sqrt(MSE);

            %Hitung PSNR
            PSNR = 10*log((255^2)/MSE);

            fprintf('\nFilter Rata-Rata\n')
            fprintf('MSE  : %12.4f;\n',MSE)
            fprintf('RMSE : %12.4f;\n',RMSE)
            fprintf('PSNR : %12.4f;\n',PSNR)
            
     case 'median'      %Ilyasa_07316
         
        [m,n] = size(R);
        G=zeros(m,n);
        f=double (R);
         
         for x =2 : m-1
            for y=2 : n-1

            data=[f(x-1,y-1) f(x-1,y) f(x-1,y+1)    ...
            f(x,y-1) f(x,y) f(x,y+1)                ...
            f(x+1,y-1) f(x+1,y) f(x+1,y+1)];
            % Urutkan
            for i=1 : 8
                for j=i+1 : 9
                    if data(i) > data(j)
                    tmp = data(i);
                    data(i) = data(j);
                    data(j) = tmp;
                    end
                end
            end
            % Ambil nilai median
            G(x,y) = data(5);
            end
        end
            G=uint8(G);
            subplot(1,2,1); imshow(R); title('Citra Asli');
            subplot(1,2,2); imshow(G); title('Filter Median');
            
            %Hitung MSE, RMSE, PNSR
            [baris,kolom]=size(R);
            selisih=zeros(baris,kolom);

            for i=1:baris
                for j=1:kolom
                    selisih(i,j)=G(i,j)-R(i,j);


                end
            end

            %hitung kuadrat
            kuadrat=zeros(baris,kolom);
            for i=1:baris
                for j=1:kolom
                    kuadrat(i,j)=selisih(i,j)*selisih(i,j);

                end
            end

            kuadrat;
            totalkuwadrat=0;

            for i=1:baris
                for j=1:kolom
                    totalkuwadrat=totalkuwadrat+kuadrat(i,j);

                end
            end

            totalkuwadrat;

            MSE=totalkuwadrat/(baris*kolom);

            % Hitung RMSE
            RMSE = sqrt(MSE);

            %Hitung PSNR
            PSNR = 10*log((255^2)/MSE);

            fprintf('\nFilter Median\n')
            fprintf('MSE  : %12.4f;\n',MSE)
            fprintf('RMSE : %12.4f;\n',RMSE)
            fprintf('PSNR : %12.4f;\n',PSNR)
            
     case 'konvolusi'       %Ilyasa_07316
         
        [m,n] = size(R);
        G=zeros(m,n);
        f=double (R);
         
        h1= -1; %h(x-1,y-1)
        h2= 1; %h(x-1,y)
        h3=1; %h(x-1,y+1)
        h4=1; %h(x,y-1);
        h5=-2; %h(x,y);
        h6=1; %h(x,y+1);
        h7=-1; %h(x+1,y-1);
        h8=-1; %h(x+1,y);
        h9=1; %h(x+1,y+1);

        for x =2 : m-1
            for y=2 : n-1
            G(x,y)=h1*f(x+1,y+1)+ h2*f(x+1,y)+ h3*f(x+1,y-1)+   ...
            h4*f(x,y+1)+ h5*f(x,y)+ h6*f(x,y-1)+                ...
            h7* f(x-1,y+1)+ h8*f(x-1,y) + h9*f(x-1,y-1);
            %G(x,y) = uint8(G(x,y));
            end
        end
        G=uint8(G);
        subplot(1,2,1); imshow(R); title('Citra Asli');
        subplot(1,2,2); imshow(G); title('Filter Konvolusi');
        
        %Hitung MSE, RMSE, PNSR
        [baris,kolom]=size(R);
            selisih=zeros(baris,kolom);

            for i=1:baris
                for j=1:kolom
                    selisih(i,j)=G(i,j)-R(i,j);


                end
            end

            %hitung kuadrat
            kuadrat=zeros(baris,kolom);
            for i=1:baris
                for j=1:kolom
                    kuadrat(i,j)=selisih(i,j)*selisih(i,j);

                end
            end

            kuadrat;
            totalkuwadrat=0;

            for i=1:baris
                for j=1:kolom
                    totalkuwadrat=totalkuwadrat+kuadrat(i,j);

                end
            end

            totalkuwadrat;

            MSE=totalkuwadrat/(baris*kolom);

            % Hitung RMSE
            RMSE = sqrt(MSE);

            %Hitung PSNR
            PSNR = 10*log((255^2)/MSE);

            fprintf('\nFilter Konvolusi\n')
            fprintf('MSE  : %12.4f;\n',MSE)
            fprintf('RMSE : %12.4f;\n',RMSE)
            fprintf('PSNR : %12.4f;\n',PSNR)
        
     case 'lolos-rendah'     %Ilyasa_07316
         
        [m,n] = size(R);
        G=zeros(m,n);
        f=double (R);
         
        h1= 0; %h(x-1,y-1)
        h2= 1; %h(x-1,y)
        h3=0; %h(x-1,y+1)
        h4=1; %h(x,y-1);
        h5=2; %h(x,y);
        h6=1; %h(x,y+1);
        h7=0; %h(x+1,y-1);
        h8=1; %h(x+1,y);
        h9=0; %h(x+1,y+1);
        
        for x =2 : m-1
            for y=2 : n-1

                G(x,y)=1/6*(h1*f(x+1,y+1)+ h2*f(x+1,y)+ h3*f(x+1,y-1)+  ...
                h4*f(x,y+1)+ h5*f(x,y)+ h6*f(x,y-1)+                    ...
                h7* f(x-1,y+1)+ h8*f(x-1,y) + h9*f(x-1,y-1));
                G(x,y) = uint8(G(x,y));
            end
        end

        G=uint8(G);
        subplot(1,2,1); imshow(R); title('Citra Asli');
        subplot(1,2,2); imshow(G); title('Filter Lolos Rendah');
        
        %Hitung MSE, RMSE, PNSR
        [baris,kolom]=size(R);
            selisih=zeros(baris,kolom);

            for i=1:baris
                for j=1:kolom
                    selisih(i,j)=G(i,j)-R(i,j);


                end
            end

            %hitung kuadrat
            kuadrat=zeros(baris,kolom);
            for i=1:baris
                for j=1:kolom
                    kuadrat(i,j)=selisih(i,j)*selisih(i,j);

                end
            end

            kuadrat;
            totalkuwadrat=0;

            for i=1:baris
                for j=1:kolom
                    totalkuwadrat=totalkuwadrat+kuadrat(i,j);

                end
            end

            totalkuwadrat;

            MSE=totalkuwadrat/(baris*kolom);

            % Hitung RMSE
            RMSE = sqrt(MSE);

            %Hitung PSNR
            PSNR = 10*log((255^2)/MSE);

            fprintf('\nFilter Lolos Rendah\n')
            fprintf('MSE  : %12.4f;\n',MSE)
            fprintf('RMSE : %12.4f;\n',RMSE)
            fprintf('PSNR : %12.4f;\n',PSNR)
        
     case 'lolos-tinggi'     %Ilyasa_07316
         
        [m,n] = size(R);
        G=zeros(m,n);
        f=double (R);
         
        h1= 0; %h(x-1,y-1)
        h2= -1; %h(x-1,y)
        h3=0; %h(x-1,y+1)
        h4=-1; %h(x,y-1);
        h5=4; %h(x,y);
        h6=-1; %h(x,y+1);
        h7=0; %h(x+1,y-1);
        h8=-1; %h(x+1,y);
        h9=0; %h(x+1,y+1);
        
        for x =2 : m-1
            for y=2 : n-1

                G(x,y)=h1*f(x+1,y+1)+ h2*f(x+1,y)+ h3*f(x+1,y-1)+   ...
                h4*f(x,y+1)+ h5*f(x,y)+ h6*f(x,y-1)+                ...
                h7* f(x-1,y+1)+ h8*f(x-1,y) + h9*f(x-1,y-1);
                G(x,y) = uint8(G(x,y));
            end
        end

        G=uint8(G);
        subplot(1,2,1); imshow(R); title('Citra Asli');
        subplot(1,2,2); imshow(G); title('Filter Lolos Tinggi');
        
        %Hitung MSE, RMSE, PNSR
        [baris,kolom]=size(R);
            selisih=zeros(baris,kolom);

            for i=1:baris
                for j=1:kolom
                    selisih(i,j)=G(i,j)-R(i,j);


                end
            end

            %hitung kuadrat
            kuadrat=zeros(baris,kolom);
            for i=1:baris
                for j=1:kolom
                    kuadrat(i,j)=selisih(i,j)*selisih(i,j);

                end
            end

            kuadrat;
            totalkuwadrat=0;

            for i=1:baris
                for j=1:kolom
                    totalkuwadrat=totalkuwadrat+kuadrat(i,j);

                end
            end

            totalkuwadrat;

            MSE=totalkuwadrat/(baris*kolom);

            % Hitung RMSE
            RMSE = sqrt(MSE);

            %Hitung PSNR
            PSNR = 10*log((255^2)/MSE);

            fprintf('\nFilter Lolos Tinggi\n')
            fprintf('MSE  : %12.4f;\n',MSE)
            fprintf('RMSE : %12.4f;\n',RMSE)
            fprintf('PSNR : %12.4f;\n',PSNR)
        
     case 'hboost'      %Ilyasa_07316
         
        [m,n] = size(R);
        G=zeros(m,n);
        f=double (R);
         
        h1= -1; %h(x-1,y-1)
        h2= -1; %h(x-1,y)
        h3=-1; %h(x-1,y+1)
        h4=-1; %h(x,y-1);
        h5=11; %h(x,y);
        h6=-1; %h(x,y+1);
        h7=-1; %h(x+1,y-1);
        h8=-1; %h(x+1,y);
        h9=-1; %h(x+1,y+1);
        
        for x =2 : m-1
            for y=2 : n-1

                G(x,y)=h1*f(x+1,y+1)+ h2*f(x+1,y)+ h3*f(x+1,y-1)+   ...
                h4*f(x,y+1)+ h5*f(x,y)+ h6*f(x,y-1)+                ...
                h7* f(x-1,y+1)+ h8*f(x-1,y) + h9*f(x-1,y-1);
                G(x,y) = uint8(G(x,y));
            end
        end

        G=uint8(G);
        subplot(1,2,1); imshow(R); title('Citra Asli');
        subplot(1,2,2); imshow(G); title('Filter High Boost');
        
        %Hitung MSE, RMSE, PNSR
        [baris,kolom]=size(R);
            selisih=zeros(baris,kolom);

            for i=1:baris
                for j=1:kolom
                    selisih(i,j)=G(i,j)-R(i,j);


                end
            end

            %hitung kuadrat
            kuadrat=zeros(baris,kolom);
            for i=1:baris
                for j=1:kolom
                    kuadrat(i,j)=selisih(i,j)*selisih(i,j);

                end
            end

            kuadrat;
            totalkuwadrat=0;

            for i=1:baris
                for j=1:kolom
                    totalkuwadrat=totalkuwadrat+kuadrat(i,j);

                end
            end

            totalkuwadrat;

            MSE=totalkuwadrat/(baris*kolom);

            % Hitung RMSE
            RMSE = sqrt(MSE);

            %Hitung PSNR
            PSNR = 10*log((255^2)/MSE);

            fprintf('\nFilter High Boost\n')
            fprintf('MSE  : %12.4f;\n',MSE)
            fprintf('RMSE : %12.4f;\n',RMSE)
            fprintf('PSNR : %12.4f;\n',PSNR)
        
     case 'emboss'      %Ilyasa_07316
         
        [m,n] = size(R);
        G=zeros(m,n);
        f=double (R);
         
        h1= -1; %h(x-1,y-1)
        h2= 0; %h(x-1,y)
        h3= 0; %h(x-1,y+1)
        h4= 0; %h(x,y-1);
        h5= 0; %h(x,y);
        h6= 0; %h(x,y+1);
        h7= 0; %h(x+1,y-1);
        h8= 0; %h(x+1,y);
        h9= 1; %h(x+1,y+1);
        
        for x =2 : m-1
            for y=2 : n-1

                G(x,y)=h1*f(x+1,y+1)+ h2*f(x+1,y)+ h3*f(x+1,y-1)+   ...
                h4*f(x,y+1)+ h5*f(x,y)+ h6*f(x,y-1)+                ...
                h7* f(x-1,y+1)+ h8*f(x-1,y) + h9*f(x-1,y-1);
                G(x,y) = uint8(G(x,y));
            end
        end

        G=uint8(G);
        subplot(1,2,1); imshow(R); title('Citra Asli');
        subplot(1,2,2); imshow(G); title('Filter Emboss');
        
        %Hitung MSE, RMSE, PNSR
        [baris,kolom]=size(R);
            selisih=zeros(baris,kolom);

            for i=1:baris
                for j=1:kolom
                    selisih(i,j)=G(i,j)-R(i,j);


                end
            end

            %hitung kuadrat
            kuadrat=zeros(baris,kolom);
            for i=1:baris
                for j=1:kolom
                    kuadrat(i,j)=selisih(i,j)*selisih(i,j);

                end
            end

            kuadrat;
            totalkuwadrat=0;

            for i=1:baris
                for j=1:kolom
                    totalkuwadrat=totalkuwadrat+kuadrat(i,j);

                end
            end

            totalkuwadrat;

            MSE=totalkuwadrat/(baris*kolom);

            % Hitung RMSE
            RMSE = sqrt(MSE);

            %Hitung PSNR
            PSNR = 10*log((255^2)/MSE);

            fprintf('\nFilter Emboss\n')
            fprintf('MSE  : %12.4f;\n',MSE)
            fprintf('RMSE : %12.4f;\n',RMSE)
            fprintf('PSNR : %12.4f;\n',PSNR)
        
     %case 'konvolusi'
end