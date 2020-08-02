#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define TTref_sec 43135.816

struct UTC_time
{
    int year,month,day,hr,min;
    double sec;
} d_time;

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

int main(int argc, char const *argv[])
{
    //struct UTC_time c_time={2000,1,1,11,58,55.816};
    TTsec2UTC(615553864.184);

    printf("%d %s %d @ %d:%d:%.3lf\n",d_time.year,monv[d_time.month-1],d_time.day,d_time.hr,d_time.min,d_time.sec);
    return 0;
}
