BEGIN {imin=100.0; imax=-100.0; jmin=100.0; jmax=-100.0; imean = 0.0; jmean = 0.0;}
/EIGENVALUES/ { Ix=$3;
                Iy=$4;
                Iz=$5;
                Tr = Ix+Iy+Iz;
                I2=Ix*Iy+Iy*Iz+Iz*Ix;
                L3=Tr/3.;
                Jx=(Ix-L3);
                Jy=(Iy-L3);
                Jz=(Iz-L3);
                ii=1-3.*I2/Tr/Tr;
                jj=27*Jx*Jy*Jz/Tr**3;
                if (ii > imax) imax=ii;
                if (jj > jmax) jmax=jj;
                if (ii < imin) imin=ii;
                if (jj < jmin) jmin=jj;
                n=n+1;
                imean = imean+ii;
                jmean = jmean+jj;
                imean2 = imean2+ii*ii;
                jmean2 = jmean2+jj*jj;
                printf "%8.3f %8.3f\n", ii,jj;
              }
END { printf "# I = %8.3f %8.3f +/- %8.3f %8.3f J = %8.3f %8.3f +/- %8.3f %8.3f\n", imin, imean/n, sqrt(imean2/n-imean*imean/n/n), imax, jmin, jmean/n, sqrt(jmean2/n-jmean*jmean/n/n), jmax; }
