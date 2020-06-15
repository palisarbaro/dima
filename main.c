#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include <stdlib.h>
int N = 2; 
double *fi; // угол поворота паралеллипипеда
double *width, *length, *height; 
double *x_c, *y_c, *z_c; // центр
double *result1_x, *result1_y, *result1_z;
double *result2_x, *result2_y, *result2_z;
int *power;
int *td;
double A[6], B[6], C[6], D[6];
#define new(n,t) (t *)calloc(n,sizeof(t))
double pi = 3.1415926535897932;
double mmm = 0.0000000000001;
#define eq(x,y) (abs(x-y)<mmm)
struct obj{
    uint16_t dist:9;
    uint16_t size:4;
    uint16_t pow:3;
};
struct codegram{
    uint16_t td;
    uint16_t az;
    uint16_t um;
    struct obj objs[5];
};
int M = 10;
struct codegram* codegrams;
int cd_iter;
void alloc(){
    fi = new(N,double);
    width = new(N,double);
    length = new(N,double);
    height = new(N,double);
    x_c = new(N,double);
    y_c = new(N,double);
    z_c = new(N,double);
    result1_x = new(N,double);
    result1_y = new(N,double);
    result1_z = new(N,double);
    result2_x = new(N,double);
    result2_y = new(N,double);
    result2_z = new(N,double);
    power = new(N,int);
    td = new(N,int);
    codegrams = new(M, struct codegram);
    memset(codegrams,0,sizeof(struct codegram)*M);
}
void transform_cube_to_parallelepiped(int i, double *x, double *y, double *z){
    for(int p=0;p<8;p++){
        double nx,ny,nz;
        nx = x_c[i] + x[p] * width[i]  * cos(fi[i]) - y[p] * width[i]  * sin(fi[i]);
        ny = y_c[i] + x[p] * length[i] * sin(fi[i]) + y[p] * length[i] * cos(fi[i]); 
        nz = z_c[i] + z[p] * height[i];
        x[p] = nx;
        y[p] = ny;
        z[p] = nz;
    }
}
double *lines[3]; 
void add(int k1, int k2, double m){
    for(int Q=0;Q<4;Q++){  
        lines[k1][Q] += m*lines[k2][Q];
    }
}
void mult(int k,double m){
    for(int Q=0;Q<4;Q++){  
        lines[k][Q]  *= m; 
    }

}
int solve_gauss(int n){
    if (eq(lines[0][0],0)) {
        add(0,1,1);
        if (eq(lines[0][0],0)) {
            add(0,2,1);
            if (eq(lines[0][0],0)) {
                printf("unlikely case \n");
                return 0;
            }
        }
    }
    
    mult(0,1./lines[0][0]);
    add(1,0,-lines[1][0]);
    add(2,0,-lines[2][0]);
    
    if (eq(lines[1][1],0)) {
        add(1,2,1);
        if (eq(lines[1][1],0)) {
            printf("unlikely case \n");
            return 0;
        }
    }
    mult(1,1./lines[1][1]);
    add(2,1,-lines[2][1]);
    if (eq(lines[2][2],0)) {
        printf("unlikely case \n");
        return 0;
    }
    mult(2,1./lines[2][2]);
    add(1,2,-lines[1][2]);
    add(0,2,-lines[0][2]);

    add(0,1,-lines[0][1]);
    /*
    for(int q = 0; q<3;q++){
        printf("%lf    %lf    %lf   %lf \n",lines[q][0],lines[q][1],lines[q][2],lines[q][3]);
    }
    printf("\n");
    */
    A[n] = lines[0][3];
    B[n] = lines[1][3];
    C[n] = lines[2][3];
    D[n] = -1; 
    /*
    printf(" x*%lf +  y*%lf  +  z*%lf  + %lf \n", A[n],B[n],C[n],D[n]);
        printf("\n");
*/
    return 1;
}
int check_point_matches_plain(int n,double x, double y, double z){
    return A[n]*x+B[n]*y+C[n]*z+D[n]>=0;
}
void get_planes(int i, double *x, double *y, double *z){
    
    int planes_points[] = {
        0,1,5,
        3,2,6,
        0,1,2,
        4,5,6,
        0,3,7,
        1,2,6
    };
    for(int n=0;n<6;n++){
        int p1 = planes_points[3*n];
        int p2 = planes_points[3*n+1];
        int p3 = planes_points[3*n+2];
        double line1[4] = {x[p1], y[p1], z[p1], 1};
        double line2[4] = {x[p2], y[p2], z[p2], 1};
        double line3[4] = {x[p3], y[p3], z[p3], 1};
        
        lines[0] = line1;
        lines[1] = line2;
        lines[2] = line3; 
        if(! solve_gauss(n)){
            double line1[4] = {x[p1]+0.01, y[p1]+0.01, z[p1]+0.01, 1};
            double line2[4] = {x[p2]+0.01, y[p2]+0.01, z[p2]+0.01, 1};
            double line3[4] = {x[p3]+0.01, y[p3]+0.01, z[p3]+0.01, 1};
            
            lines[0] = line1;
            lines[1] = line2;
            lines[2] = line3; 
            solve_gauss(n);
        }
    }
}

void find_intersections(int i, double azimut, double mesto){
    double a = cos(azimut)*sin(mesto);
    double b = sin(azimut)*sin(mesto);
    double c = cos(mesto); 
    double result_x[2];
    double result_y[2];
    double result_z[2];

    int s = 0;
    for (int n=0;n<6;n++){
        int t = -D[n]/(A[n]*a + B[n]*b + C[n]*c);
        double x = a*t;
        double y = b*t;
        double z = c*t; 
        int f = 1;
        for(int m=0;m<6;m++){
            if (m==n) continue;
            if(!check_point_matches_plain(m,x,y,z)){
                f = 0;
            }
        }
        if(f){ //the point on the parallelepiped
            result_x[s] = x;
            result_y[s] = y;
            result_z[s] = z;
            s++;
        }
    }
    if(s>0){
        if(s!=2) printf("Oops s=%d",s);
        printf("(%lf,%lf,%lf) \n", result_x[0],result_y[0],result_z[0]);
        printf("(%lf,%lf,%lf) \n", result_x[1],result_y[1],result_z[1]);
        result1_x[i] = result_x[0];
        result1_y[i] = result_y[0];
        result1_z[i] = result_z[0]; 
        result2_x[i] = result_x[1];
        result2_y[i] = result_y[1];
        result2_z[i] = result_z[1]; 
    }
    else{
        result1_x[i] = -1;
        result1_y[i] = -1;
        result1_z[i] = -1; 
        result2_x[i] = -1;
        result2_y[i] = -1;
        result2_z[i] = -1;
        printf("miss \n");
    }
}

void process(double azimut, double mesto) {
    for(int i=0;i<N;i++){
        double parallelepiped_x[] = {-1,  1,  1, -1, -1,  1, 1, -1};
        double parallelepiped_y[] = {-1, -1,  1,  1, -1, -1, 1,  1};
        double parallelepiped_z[] = {-1, -1, -1, -1,  1,  1, 1,  1};
        
        transform_cube_to_parallelepiped(i,parallelepiped_x,parallelepiped_y,parallelepiped_z);
        
        get_planes(i,parallelepiped_x,parallelepiped_y,parallelepiped_z);
        for(int n=0;n<6;n++){
            if(!check_point_matches_plain(n,x_c[i],y_c[i],z_c[i])){
                A[n] *= -1;
                B[n] *= -1;
                C[n] *= -1;
                D[n] *= -1;
            }
        }
        find_intersections(i,azimut,mesto);
        /*
        for(int i=0;i<8;i++){
            printf("%lf, %lf, %lf \n",parallelepiped_x[i],parallelepiped_y[i],parallelepiped_z[i]);
        }
        */
    }
}
int set_codegram(int dist,int size,int i,int k,int az,int um,int td){
    codegrams[cd_iter].az = az;
    codegrams[cd_iter].um = um;
    codegrams[cd_iter].td = td;
    codegrams[cd_iter].objs[k].dist = dist;
    codegrams[cd_iter].objs[k].pow = power[i];
    codegrams[cd_iter].objs[k].size = size;
    k++;
    if(k==5){
        k=0;
        cd_iter++;
    }
    return k;
}
void write(int dd){
    cd_iter = 0;
    for(int i=0;i<N;i++){
        if (result1_x[i]==-1 || result2_x[i]==-1){
            continue;
        }
        double len1 = sqrt( result1_x[i]*result1_x[i] + result1_y[i]*result1_y[i] + result1_z[i]*result1_z[i] );
        double len2 = sqrt( result2_x[i]*result2_x[i] + result2_y[i]*result2_y[i] + result2_z[i]*result2_z[i] );
        if(len1>len2){
            double tmp = len1;
            len1 = len2;
            len2 = tmp;
        }
        int az = 2;
        int um = 4;
        int dist = ceil(len1/dd);
        int size = floor(len2/dd) - dist + 1;
        if(dist*dd>len2){
            printf("QQ");
            continue;
        }
        int k = 0;
        while(size>15){
            k = set_codegram(dist,15,i,k,az,um,td[i]);
            size -= 15;
            dist+=15;
        }
        set_codegram(dist,size,i,k,az,um,td[i]);
        cd_iter++;
    }
}
void print_c(){
    for(int i=0;i<cd_iter;i++){
        printf("#%d \n", i);
        for(int k=0;k<5;k++){
              printf("dist: %d      ",codegrams[i].objs[k].dist);
        }
        printf("\n");
        for(int k=0;k<5;k++){
              printf("size: %d      ",codegrams[i].objs[k].size);
        }
        printf("\n");

    }
}

const int T = 5;
void init(){
    for(int i=0;i<N;i++) {
        int f,w,l,h,x,y,z,pow,t,vx,vy,vz;
        scanf("%d %d %d %d %d %d %d %d %d %d %d %d", &f, &pow, &t, &w, &l, &h, &x, &y, &z, &vx, &vy, &vz);
        fi[i] = f * 2*pi/360;
        width[i] = w;
        length[i] = l;
        height[i] = h;  
        x_c[i] = x + vx*T;
        y_c[i] = y + vy*T;
        z_c[i] = z + vz*T; 
        power[i] = pow;
        td[i] = t;
    }
}
void free_f(){
    free(fi);
    free(width);
    free(length);
    free(height);
    free(x_c);
    free(y_c);
    free(z_c);
    free(result1_x);
    free(result1_y);
    free(result1_z);
    free(result2_x);
    free(result2_y);
    free(result2_z);
    free(power);
    free(td);
}
int main(){
    alloc();
    init();
    for(int i=0;i<N;i++) {
        printf("%lf %lf %lf %lf %lf %lf %lf \n", fi[i], width[i], length[i], height[i], x_c[i], y_c[i], z_c[i]);
    }
    process(0.0,0.0);
    write(10);
    print_c();
    free_f();
    return 0;
}