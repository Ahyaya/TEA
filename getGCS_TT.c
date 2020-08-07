#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define pi 3.141592653589793
#define e 0.0167
#define T 31556925.51
#define au_a 149602022500.0
#define Earth_Rotation 7.2921159e-5
#define TTref_sec 43135.816
#define theta_X -1.791513047394897
#define Fi 0.4089888213840046
#define cosFi 0.9175234191066961
#define sinFi 0.3976817513926911
#define TTref_sec 43135.816

double dt,Exp_time,alpha_X0;

struct UTC_time
{
    int year,month,day,hr,min;
    double sec;
};

struct UTC_time d_time;

struct observ
{
    double t,theta,EOT;
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

struct observ ref={615553864.184,0,-263.4};
int get_alpha_X0()
{
    double td=ref.t-599572864.184, theta=ref.theta, alpha_theta, alpha_tx;
    while(td>86400){td-=86400;}
    alpha_theta=(1-td/43200)*pi;
    alpha_tx=acos(-cos(theta-theta_X)/sqrt(cos(theta-theta_X)*cos(theta-theta_X)+cos(Fi)*cos(Fi)*sin(theta-theta_X)*sin(theta-theta_X)));
    alpha_tx=sin(theta-theta_X)>0?alpha_tx:-alpha_tx;
    alpha_X0=alpha_theta+alpha_tx;
    return 0;
}

int ecs2gcs(double *lon, double *lat, double t, double X, double Y, double Z)
{
    double alpha_ER=Earth_Rotation*(t-ref.t), alpha_EOTdrift=Earth_Rotation*ref.EOT;
    double clockwise=cosFi*Y+sinFi*Z, alpha_RX=acos(X/sqrt(X*X+clockwise*clockwise));
    *lon=clockwise>0?(alpha_X0-alpha_ER-alpha_EOTdrift+alpha_RX):(alpha_X0-alpha_ER-alpha_EOTdrift-alpha_RX);
    *lat=asin((cosFi*Z-sinFi*Y)/sqrt(X*X+Y*Y+Z*Z));
    return 0;
}

int TTsec2UTC(double t)
{
    struct UTC_time time={2000,1,1,11,58,55.816};
    int pf=0,yearType[2]={365,366},flag[5]={0};
    int daylist[2][13]={0,31,28,31,30,31,30,31,31,30,31,30,31,0,31,29,31,30,31,30,31,31,30,31,30,31};
    double temp=86400.0*yearType[isLeapYear(time.year)];
    while(t>temp){time.year++;t-=temp;temp=86400.0*yearType[isLeapYear(time.year)];}
    pf=isLeapYear(time.year);temp=86400.0*daylist[pf][time.month];
    while(t>temp&&time.month<12){time.month++;t-=temp;temp=86400.0*daylist[pf][time.month];}
    temp=0;pf=daylist[pf][time.month];
    while(t>86400&&time.day<pf){time.day++;t-=86400;}
    while(t>3600&&time.hr<23){time.hr++;t-=3600;}
    while(t>60&&time.min<59){time.min++;t-=60;}
    time.sec+=t;
    while(flag[0]*flag[1]*flag[2]*flag[3]*flag[4]<1)
    {
        if(time.sec<60){flag[0]=1;}else{flag[0]=0;time.sec-=60;time.min++;}
        if(time.min<60){flag[1]=1;}else{flag[1]=0;time.min-=60;time.hr++;}
        if(time.hr<24){flag[2]=1;}else{flag[2]=0;time.hr-=24;time.day++;}
        if(time.day<=daylist[isLeapYear(time.year)][time.month]){flag[3]=1;}else{flag[3]=0;time.day=1;time.month++;}
        if(time.month<=12){flag[4]=1;}else{flag[4]=0;time.month=1;time.year++;}
    }
    d_time=time;
    return 0;
}

int message_print(double t, double alpha, double delta)
{
    while(alpha+pi<0){alpha+=2*pi;}
    char Hor[2]="WE",Ver[2]="SN";
    int degtran_1[3],degtran_2[3];
    int sg_1=(alpha>0)?(1):(0);
    int sg_2=(delta>0)?(1):(0);

    alpha=fabs(alpha);delta=fabs(delta);
    degtran_1[0]=(int) (alpha*180/pi);degtran_1[1]=(int) 60*(alpha*180/pi-degtran_1[0]);degtran_1[2]=(int) 60*(60*(alpha*180/pi-degtran_1[0])-degtran_1[1]);
    degtran_2[0]=(int) (delta*180/pi);degtran_2[1]=(int) 60*(delta*180/pi-degtran_2[0]);degtran_2[2]=(int) 60*(60*(delta*180/pi-degtran_2[0])-degtran_2[1]);
    printf("Terrestrial Time footprint: %.3lf\n\n",t);  
    printf("delta: %d°%d\'%d\" %c\n",degtran_2[0],degtran_2[1],degtran_2[2],Ver[sg_2]);
    printf("alpha: %d°%d\'%d\" %c\n",degtran_1[0],degtran_1[1],degtran_1[2],Hor[sg_1]);
    return 0;
}

int usage_print()
{
    printf("Usage: getGCS [option]\n\t[-h] or [--help] Display this information.\n");
    printf("\t[-i filename] Specify the input file to convert.\n");
    printf("\t[-o filename] Costomize output filename (optional).\n");
    return 0;
}

int main(int argc, char const *argv[])
{
    double TT_s,TT_d,X,Y,Z;
    int pf=0,namepf=0,nameip=0;
    FILE *outf,*inf;

    double r,theta=ref.theta,t=ref.t,alpha,delta;

    get_alpha_X0();

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
            }else{
                printf("unknown argument\n\n");
                usage_print();
                return -1;
            }
            pf++;
        }
        
    }else{usage_print();return -1;}
    if(!nameip){usage_print();return -1;}
    if(NULL==(inf=fopen(argv[nameip],"r"))){printf("!Unable to read file %s\n",argv[nameip]);return -1;}
    
    //fscanf(inf,"%d,%d,%d@%d:%d:%lf",&c_time.year,&c_time.month,&c_time.day,&c_time.hr,&c_time.min,&c_time.sec);
    fscanf(inf,"%lf%lf%lf%lf",&t,&X,&Y,&Z);
    TTsec2UTC(t);
    //t=UTC2TTsec(c_time);
    r=sqrt(X*X+Y*Y+Z*Z);
    ecs2gcs(&alpha,&delta,t,X,Y,Z);

    TT_s=t;
    printf("Started at: %d %s %d @ %02d:%02d:%06.3lf UTC+0\n\n",d_time.year,monv[d_time.month-1],d_time.day,d_time.hr,d_time.min,d_time.sec);

    message_print(TT_s,alpha,delta);

    if(namepf){outf=fopen(argv[namepf],"w");}else{outf=fopen("SC_result","w");}
    dt=30;
    while(fscanf(inf,"%lf%lf%lf%lf",&t,&X,&Y,&Z)>3)
    {
        r=sqrt(X*X+Y*Y+Z*Z);
        ecs2gcs(&alpha,&delta,t,X,Y,Z);
        
        fprintf(outf,"%.3lf\t%lf\t%.10lf\t%.10lf\n",t,r,alpha,delta);
    }
    TT_d=t;

    printf("\n====================================================================\n\n");
    TTsec2UTC(TT_d);
    printf("Destination: %d %s %d @ %02d:%02d:%06.3lf UTC+0\n\n",d_time.year,monv[d_time.month-1],d_time.day,d_time.hr,d_time.min,d_time.sec);
    message_print(TT_d,alpha,delta);

    
    fclose(outf);
    fclose(inf);

    printf("Done!\n\n");
    
    return 0;
}
