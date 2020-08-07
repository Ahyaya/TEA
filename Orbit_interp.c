#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TTref_sec 43135.816

int interp_spline(double *interp_out, double *T, double *X, double *DX, int head, int tail, double t)
{
    double t_tn, t_tn1, tn_tn1;
    int pf=head;
    while(pf<tail)
    {
        if(((*(T+pf))-t)*((*(T+pf+1))-t)<1e-8)
        {
            t_tn=t-(*(T+pf));
            t_tn1=t-(*(T+pf+1));
            tn_tn1=(*(T+pf))-(*(T+pf+1));
            *interp_out=((*(X+pf))*((1-2*t_tn/tn_tn1)*t_tn1*t_tn1)+(*(X+pf+1))*((1+2*t_tn1/tn_tn1)*t_tn*t_tn)+(*(DX+pf))*(t_tn*t_tn1*t_tn1)+(*(DX+pf+1))*(t_tn1*t_tn*t_tn))/tn_tn1/tn_tn1;
            return 1;
        }
        pf++;
    }
    printf("no match unit within index\n");
    return 0;
}

struct UTC_time
{
    int year,month,day,hr,min;
    double sec;
};

char *monv[12]={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};

int isLeapYear(int year)
{
    if(year%100)
    {
        if(year%4){return 0;}else{return 1;}
    }else{
        if(year%400){return 0;}else{return 1;}
    }
    return 0;
}

double UTC2TTsec(struct UTC_time time)
{
    int pf_year=2000,pf_month=1,temp_day=0,pf=isLeapYear(time.year);
    int daylist[2][13]={0,31,28,31,30,31,30,31,31,30,31,30,31,0,31,29,31,30,31,30,31,31,30,31,30,31};
    double TTsec=time.sec+60.0*time.min+3600.0*time.hr+86400.0*(time.day-1)-TTref_sec;
    while(pf_month<time.month)
    {
        temp_day+=daylist[pf][pf_month++];
    }
    while(pf_year<time.year)
    {
        temp_day+=(365+isLeapYear(pf_year++));
    }
    TTsec+=86400.0*temp_day;

    return TTsec;
}

int usage_print()
{
    printf("Usage: Orbit_interp [option]\n\t[-h] or [--help] Display this information.\n");
    printf("\t[-i filename] Specify the input file to convert.\n");
    printf("\t[-o filename] Costomize output filename (optional).\n");
    printf("\t[-s number] or [--skip number] Skip the first number of lines (optional).Default value is 0.\n");
    printf("\t[-n number] Total number of lines intend to load (optional).Default value (also max value) is 1024.\n");
    return 0;
}

int main(int argc, char const *argv[])
{
    struct UTC_time c_time;
    
    double T[1024], X[1024], Y[1024], Z[1024], DX[1024], DY[1024], DZ[1024];
    double t, temp_X, temp_Y, temp_Z;
    FILE *inf, *outf;
    int pf,nameip,namepf,skip=0,total_line=1024,pf_EOF=0;

    if(argc>1)
    {
        pf=1;
        while(pf<argc)
        {
            
            if(!strcmp(argv[pf],"-h"))
            {
                usage_print();
                return 0;
            }else if(!strcmp(argv[pf],"-help")){
                usage_print();
                return 0;
            }else if(!strcmp(argv[pf],"--help")){
                usage_print();
                return 0;
            }else if(!strcmp(argv[pf],"-i")){
                pf++;nameip=pf;
            }else if(!strcmp(argv[pf],"-o")){
                pf++;namepf=pf;
            }else if(!strcmp(argv[pf],"-s")){
                skip=atoi(argv[++pf]);
            }else if(!strcmp(argv[pf],"--skip")){
                skip=atoi(argv[++pf]);
            }else if(!strcmp(argv[pf],"-n")){
                total_line=atoi(argv[++pf]);
            }else{
                printf("unknown argument\n\n");
                usage_print();
                return -1;
            }
            pf++;
        }
        
    }else{usage_print();return -1;}

    if(!nameip){usage_print();return -1;}
    if(NULL==(inf=fopen(argv[nameip],"r"))){printf("Unable to read file %s\n",argv[nameip]);return -1;}
    if(namepf){outf=fopen(argv[namepf],"w");}else{outf=fopen("orbit_ouput","w");}
    if(skip){
        for(pf=0;pf<skip;pf++)
        {
            fscanf(inf,"%*d,%*d,%*d@%*d:%*d:%*lf");
            fscanf(inf,"%*lf%*lf%*lf%*lf%*lf%*lf");
        }
    }
    for(pf=0;pf<total_line;pf++)
    {
        if(fscanf(inf,"%d,%d,%d@%d:%d:%lf",&c_time.year,&c_time.month,&c_time.day,&c_time.hr,&c_time.min,&c_time.sec)<5)
        {
            pf_EOF=pf;printf("EOF incident at ttl: %d\n",pf);break;
        }
        fscanf(inf,"%lf%lf%lf%lf%lf%lf",&X[pf],&Y[pf],&Z[pf],&DX[pf],&DY[pf],&DZ[pf]);
        T[pf]=UTC2TTsec(c_time);
    }
    if(pf_EOF) total_line=pf_EOF;

    for(t=T[0],pf=0;t<T[total_line-1];t+=30)
    {
        interp_spline(&temp_X,T,X,DX,0,total_line-1,t);
        interp_spline(&temp_Y,T,Y,DY,0,total_line-1,t);
        interp_spline(&temp_Z,T,Z,DZ,0,total_line-1,t);
        fprintf(outf,"%lf\t%lf\t%lf\t%lf\n",t,temp_X,temp_Y,temp_Z);
    }
    fclose(outf);
    fclose(inf);


    printf("Done!\n");
    return 0;

}
