#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define TTref_sec 43135.816

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

int main(int argc, char const *argv[])
{
    //Setup c_time as [int year],[int month],[int day],[int hr],[int min],[double sec]
    struct UTC_time c_time={2019,1,1,0,0,0};

    printf("%d %s %d @ %d:%d:%.3lf\n",c_time.year,monv[c_time.month-1],c_time.day,c_time.hr,c_time.min,c_time.sec);
    printf("%.3lf\n",UTC2TTsec(c_time));
    return 0;
}
