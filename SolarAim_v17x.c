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

double dt,Exp_time,alpha_X0;

struct UTC_time
{
    int year,month,day,hr,min;
    double sec;
};

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

double K_theta(double theta)
{
    return(2*pi/T*pow(1-e*e,-1.5)*(1-e*cos(theta))*(1-e*cos(theta)));
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

int usage_print()
{
    printf("Usage: SolarAim [option]\n\t[-h] or [--help] Display this information.\n");
    printf("\t[-l steps_limit] or [--limit steps_limit] redefine the computational steps to reach starting point.\n");
    printf("\t[-o filename] Place the intermedium results into <filename>.\n");
    printf("\t[-f] or [--force-trace] Record the trace from reference event to starting point.\n");
    return 0;
}

int main(int argc, char const *argv[])
{
    double k1,k2,k3,k4,TT_s,TT_d,X,Y;
    int pf=0,pb=0,namepf=0,multipler=1,steplimit=1000000,degtran_0[3],degtran_1[3],degtran_2[3],sg_1=1,sg_2=1,FORCE_TRACE=0;
    FILE *outf,*inf;

    double r,theta=ref.theta,t=ref.t,alpha,delta;

    struct UTC_time c_time={2019,7,4,15,8,0},d_time={2034,4,25,15,8,0};

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
            }else if(!strcmp(argv[pf],"-l")){
                pf++;steplimit=atoi(argv[pf]);
            }else if(!strcmp(argv[pf],"--limit")){
                pf++;steplimit=atoi(argv[pf]);
            }else if(!strcmp(argv[pf],"-o")){
                pf++;namepf=pf;
            }else if(!strcmp(argv[pf],"--force-trace")){
                FORCE_TRACE=1;
            }else if(!strcmp(argv[pf],"-f")){
                FORCE_TRACE=1;
            }else{
                usage_print();
                return 0;
            }
            pf++;
        }
        
    }
char timezone[3]="UTC";
if(NULL==(inf=fopen("settings.ini","r"))){printf("!Failed to open settings.ini\n");return -1;}
fscanf(inf,"%d,%d,%d@%d:%d:%lf %c%c%c%d",&c_time.year,&c_time.month,&c_time.day,&c_time.hr,&c_time.min,&c_time.sec,&timezone[0],&timezone[1],&timezone[2],&pf);
fscanf(inf,"%d,%d,%d@%d:%d:%lf %c%c%c%d",&d_time.year,&d_time.month,&d_time.day,&d_time.hr,&d_time.min,&d_time.sec,&timezone[0],&timezone[1],&timezone[2],&pb);
fclose(inf);


    printf("Start from : %d %s %d @ %02d:%02d:%06.3lf UTC%+d\n",c_time.year,monv[c_time.month-1],c_time.day,c_time.hr,c_time.min,c_time.sec,pf);
    TT_s=UTC2TTsec(c_time)-3600*pf;printf("Terrestrial Time: %.3lf\n\n",TT_s);pf=0;
    TT_d=UTC2TTsec(d_time)-3600*pb;
    if(TT_d<TT_s){printf("Time error: destination smaller than origin!\n");return -1;}

    multipler=(int) fabs(TT_s-t)/T;
    steplimit*=++multipler;printf("Multipler=%d\n",multipler);
    Exp_time=TT_s-t;

    if(FORCE_TRACE){printf("Force Trace mode activated!\n");outf=fopen("SA_trace","w");fprintf(outf,"This file record the results from the reference event to the input starting point\n\nSpecified meaning of each column: \n==============================\nt: Terrestial time in the units of second\nr: Distance to the sun in the units of meter\ntheta: Postion angle ahead the Aphelion point in the units of rad\nlon: Longitude of the subsolar point in the units of rad +(E)/-(W)\nlat: Latitude of the subsolar point in the units of rad +(N)/-(S)\n==============================\nt\tr\ttheta\tlon\tlat\n\n\n");}
    dt=Exp_time/steplimit;
    for(pf=0;pf<steplimit;pf++)
    {
        k1=dt*K_theta(theta);
        k2=dt*K_theta(theta+0.5*k1);
        k3=dt*K_theta(theta+0.5*k2);
        k4=dt*K_theta(theta+k3);
        theta+=(k1/6+k2/3+k3/3+k4/6);
        t+=dt;
        r=au_a*(1-e*e)/(1-e*cos(theta));
        
        X=-cos(theta-theta_X);
        Y=-sin(theta-theta_X);
        ecs2gcs(&alpha,&delta,t,X,Y,0);
        
        if(FORCE_TRACE){fprintf(outf,"%.3lf\t%lf\t%.10lf\t%.10lf\t%.10lf\n",t,r,theta,alpha,delta);}
    }

    if(FORCE_TRACE){fclose(outf);printf("Trace from reference event saved in file <SA_trace>.\n");}

    if(namepf){outf=fopen(argv[namepf],"w");}else{outf=fopen("SA_result","w");}
    fprintf(outf,"%.3lf\t%lf\t%.10lf\t%.10lf\t%.10lf\n",t,r,theta,alpha,delta);

    while(theta>2*pi){theta-=2*pi;}
    while(alpha+pi<0){alpha+=2*pi;}

    char Hor[2]="WE",Ver[2]="SN";
    sg_1=(alpha>0)?(1):(0);
    sg_2=(delta>0)?(1):(0);
    alpha=fabs(alpha);delta=fabs(delta);
    degtran_0[0]=(int) (theta*180/pi);degtran_0[1]=(int) 60*(theta*180/pi-degtran_0[0]);degtran_0[2]=(int) 60*(60*(theta*180/pi-degtran_0[0])-degtran_0[1]);
    degtran_1[0]=(int) (alpha*180/pi);degtran_1[1]=(int) 60*(alpha*180/pi-degtran_1[0]);degtran_1[2]=(int) 60*(60*(alpha*180/pi-degtran_1[0])-degtran_1[1]);
    degtran_2[0]=(int) (delta*180/pi);degtran_2[1]=(int) 60*(delta*180/pi-degtran_2[0]);degtran_2[2]=(int) 60*(60*(delta*180/pi-degtran_2[0])-degtran_2[1]);
    printf("Terrestrial Time footprint: %.3lf\nDistance to Sun: %.3lf km\nTheta from Aphelion: %d°%d\'%d\"\n\n",t,r/1000,degtran_0[0],degtran_0[1],degtran_0[2]);  
    printf("Delta: %d°%d\'%d\" %c\n",degtran_2[0],degtran_2[1],degtran_2[2],Ver[sg_2]);
    printf("Alpha: %d°%d\'%d\" %c\n",degtran_1[0],degtran_1[1],degtran_1[2],Hor[sg_1]);
    printf("Time step from reference to starting point: %.3lf(s)\n",dt);

    printf("\n====================================================================\n\n");
    printf("Destination: %d %s %d @ %02d:%02d:%06.3lf UTC%+d\n",d_time.year,monv[d_time.month-1],d_time.day,d_time.hr,d_time.min,d_time.sec,pb);
    printf("Terrestrial Time: %.3lf\n\n",TT_d);pb=0;

    printf("Computing...");

    dt=30;
    while(t<TT_d)
    {
        k1=dt*K_theta(theta);
        k2=dt*K_theta(theta+0.5*k1);
        k3=dt*K_theta(theta+0.5*k2);
        k4=dt*K_theta(theta+k3);
        theta+=(k1/6+k2/3+k3/3+k4/6);
        t+=dt;
        r=au_a*(1-e*e)/(1-e*cos(theta));
        
        X=-cos(theta-theta_X);
        Y=-sin(theta-theta_X);
        ecs2gcs(&alpha,&delta,t,X,Y,0);
        
        fprintf(outf,"%.3lf\t%lf\t%.10lf\t%.10lf\t%.10lf\n",t,r,theta,alpha,delta);
    }
    fclose(outf);

    printf("\b\b\b\b\b\b\b\b\b\b\b\b");
    while(theta>2*pi){theta-=2*pi;}
    while(alpha+pi<0){alpha+=2*pi;}

    sg_1=(alpha>0)?(1):(0);
    sg_2=(delta>0)?(1):(0);
    alpha=fabs(alpha);delta=fabs(delta);
    degtran_0[0]=(int) (theta*180/pi);degtran_0[1]=(int) 60*(theta*180/pi-degtran_0[0]);degtran_0[2]=(int) 60*(60*(theta*180/pi-degtran_0[0])-degtran_0[1]);
    degtran_1[0]=(int) (alpha*180/pi);degtran_1[1]=(int) 60*(alpha*180/pi-degtran_1[0]);degtran_1[2]=(int) 60*(60*(alpha*180/pi-degtran_1[0])-degtran_1[1]);
    degtran_2[0]=(int) (delta*180/pi);degtran_2[1]=(int) 60*(delta*180/pi-degtran_2[0]);degtran_2[2]=(int) 60*(60*(delta*180/pi-degtran_2[0])-degtran_2[1]);
    printf("Terrestrial Time footprint: %.3lf\nDistance to Sun: %.3lf km\nTheta from Aphelion: %d°%d\'%d\"\n\n",t,r/1000,degtran_0[0],degtran_0[1],degtran_0[2]);  
    printf("Delta: %d°%d\'%d\" %c\n",degtran_2[0],degtran_2[1],degtran_2[2],Ver[sg_2]);
    printf("Alpha: %d°%d\'%d\" %c\n",degtran_1[0],degtran_1[1],degtran_1[2],Hor[sg_1]);
    printf("Time step from starting point to destination: %.3lf(s)\n",dt);
    
    printf("Done!\n\n");
    
    return 0;
}
