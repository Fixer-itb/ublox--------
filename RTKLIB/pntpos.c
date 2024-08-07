/*------------------------------------------------------------------------------
 * pntpos.c : standard positioning
 *
 *          Copyright (C) 2007-2015 by T.TAKASU, All rights reserved.
 *
 * version : $Revision:$ $Date:$
 * history : 2010/07/28 1.0  moved from rtkcmn.c
 *                           changed api:
 *                               pntpos()
 *                           deleted api:
 *                               pntvel()
 *           2011/01/12 1.1  add option to include unhealthy satellite
 *                           reject duplicated observation data
 *                           changed api: ionocorr()
 *           2011/11/08 1.2  enable snr mask for single-mode (rtklib_2.4.1_p3)
 *           2012/12/25 1.3  add variable snr mask
 *           2014/05/26 1.4  support galileo and beidou
 *           2015/03/19 1.5  fix bug on ionosphere correction for GLO and BDS
 *-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants -----------------------------------------------------------------*/

#define SQR(x) ((x) * (x))

#define NX (4 + 3) /* # of estimated parameters */

#define MAXITR 10     /* max number of iteration for point pos */
#define ERR_ION 5.0   /* ionospheric delay std (m) */
#define ERR_TROP 3.0  /* tropspheric delay std (m) */
#define ERR_SAAS 0.3  /* saastamoinen model error std (m) */
#define ERR_BRDCI 0.5 /* broadcast iono model error factor */
#define ERR_CBIAS 0.3 /* code bias error std (m) */
#define REL_HUMI 0.7  /* relative humidity for saastamoinen model */

const double chisqr[100] = {/* chi-sqr(n) (alpha=0.001) */
                            10.8, 13.8, 16.3, 18.5, 20.5, 22.5, 24.3, 26.1, 27.9, 29.6,
                            31.3, 32.9, 34.5, 36.1, 37.7, 39.3, 40.8, 42.3, 43.8, 45.3,
                            46.8, 48.3, 49.7, 51.2, 52.6, 54.1, 55.5, 56.9, 58.3, 59.7,
                            61.1, 62.5, 63.9, 65.2, 66.6, 68.0, 69.3, 70.7, 72.1, 73.4,
                            74.7, 76.0, 77.3, 78.6, 80.0, 81.3, 82.6, 84.0, 85.4, 86.7,
                            88.0, 89.3, 90.6, 91.9, 93.3, 94.7, 96.0, 97.4, 98.7, 100,
                            101, 102, 103, 104, 105, 107, 108, 109, 110, 112,
                            113, 114, 115, 116, 118, 119, 120, 122, 123, 125,
                            126, 127, 128, 129, 131, 132, 133, 134, 135, 137,
                            138, 139, 140, 142, 143, 144, 145, 147, 148, 149};

/* pseudorange measurement error variance ------------------------------------*/
/**
 * @brief 计算导航系统伪距测量值的误差
 * 
 * @param opt 
 * @param el 
 * @param sys 
 * @return double 
 */
static double varerr(const prcopt_t *opt, double el, int sys)
{
    double fact, varr;
    //* 1、确定 sys系统的误差因子。
    fact = sys == SYS_GLO ? EFACT_GLO : (sys == SYS_SBS ? EFACT_SBS : EFACT_GPS);
    //* 2、计算由导航系统本身所带来的误差的方差。
    varr = SQR(opt->err[0]) * (SQR(opt->err[1]) + SQR(opt->err[2]) / sin(el));
    //* 3、如果 ionoopt==IONOOPT_IFLC时，IFLC模型的方差也会添加到最终计算得到的方差中。
    if (opt->ionoopt == IONOOPT_IFLC)
        varr *= SQR(3.0); /* iono-free */
    return SQR(fact) * varr;
}
/* get tgd parameter (m) -----------------------------------------------------*/
/**
 * @brief 获取TGD参数
 * 
 * @param sat 
 * @param nav 
 * @return double 
 */
static double gettgd(int sat, const nav_t *nav)
{
    int i;
    for (i = 0; i < nav->n; i++)
    {
        //* 1、从导航数据的星历中选择卫星号与 sat相同的那个星历，读取 tgd[0]参数后乘上光速。
        if (nav->eph[i].sat != sat)
            continue;
        return CLIGHT * nav->eph[i].tgd[0];
    }
    return 0.0;
}
/* psendo range with code bias correction -------------------------------------*/
/**
 * @brief 
 * 
 * obsd_t   *obs      I   observation data
 * nav_t    *nav      I   navigation data
 * double   *azel     I   对于当前定位值，每一颗观测卫星的 {方位角、高度角}
 * int      iter      I   迭代次数
 * prcopt_t *opt      I   processing options
 * double   *vare     O   伪距测量的码偏移误差
 * 返回类型：
 * double             O   最终能参与定位解算的伪距值
 */
static double prange(const obsd_t *obs, const nav_t *nav, const double *azel,
                     int iter, const prcopt_t *opt, double *var)
{
    const double *lam = nav->lam[obs->sat - 1];
    double PC, P1, P2, P1_P2, P1_C1, P2_C2, gamma;
    int i = 0, j = 1, sys;

    *var = 0.0;
    //* 1、定该卫星属于 RTKLIB设计时给定的几个导航系统之中。
    if (!(sys = satsys(obs->sat, NULL)))
        return 0.0;

    /* L1-L2 for GPS/GLO/QZS, L1-L5 for GAL/SBS */
    //* 2、如果 NFREQ>=3且该卫星属于 GAL/SBS系统，则 j=2。
    //*     而如果出现 NFREQ<2||lam[i]==0.0||lam[j]==0.0中的其中一个，直接返回 0.
    if (NFREQ >= 3 && (sys & (SYS_GAL | SYS_SBS)))
        j = 2;

    if (NFREQ < 2 || lam[i] == 0.0 || lam[j] == 0.0)
        return 0.0;

    /* test snr mask */
    if (iter > 0)
    {
        //* 3、调用 testsnr函数，测试第i、j(IFLC)个频率信号的载噪比是否符合要求。
        if (testsnr(0, i, azel[1], obs->SNR[i] * 0.25, &opt->snrmask))
        {
            //            trace(4,"snr mask: %s sat=%2d el=%.1f snr=%.1f\n",
            //                  time_str(obs->time,0),obs->sat,azel[1]*R2D,obs->SNR[i]*0.25);
            return 0.0;
        }
        //* 4、计算出γ值(f1^2/f2^2，见ICD-GPS-200C P90)，从 obs和 nav数据中读出测量伪距值和 码偏移值(?)。
        if (opt->ionoopt == IONOOPT_IFLC)
        {
            if (testsnr(0, j, azel[1], obs->SNR[j] * 0.25, &opt->snrmask))
                return 0.0;
        }
    }
    gamma = SQR(lam[j]) / SQR(lam[i]); /* f1^2/f2^2 */
    P1 = obs->P[i];
    P2 = obs->P[j];
    P1_P2 = nav->cbias[obs->sat - 1][0];
    P1_C1 = nav->cbias[obs->sat - 1][1];
    P2_C2 = nav->cbias[obs->sat - 1][2];

    /* if no P1-P2 DCB, use TGD instead */
    //* 5、从数据中读出的P1_P2==0，则使用 TGD代替，TGD值由 gettgd函数计算得到。
    if (P1_P2 == 0.0 && (sys & (SYS_GPS | SYS_GAL | SYS_QZS)))
    {
        P1_P2 = (1.0 - gamma) * gettgd(obs->sat, nav);
    }
    //* 6、如果 ionoopt==IONOOPT_IFLC，根据 obs->code的值来决定是否对 P1、P2进行修正，
    //*     之后再组合出 IFLC时的伪距值(ICD-GPS-200C P91)。
    //*     否则，则是针对单频接收即进行的数据处理。先对 P1进行修正，然后再计算出伪距值（PC）
    if (opt->ionoopt == IONOOPT_IFLC)
    { /* dual-frequency */

        if (P1 == 0.0 || P2 == 0.0)
            return 0.0;
        if (obs->code[i] == CODE_L1C)
            P1 += P1_C1; /* C1->P1 */
        if (obs->code[j] == CODE_L2C)
            P2 += P2_C2; /* C2->P2 */

        /* iono-free combination */
        PC = (gamma * P1 - P2) / (gamma - 1.0);
    }
    else
    { /* single-frequency */

        if (P1 == 0.0)
            return 0.0;
        if (obs->code[i] == CODE_L1C)
            P1 += P1_C1; /* C1->P1 */
        PC = P1 - P1_P2 / (1.0 - gamma);
    }
    //* 7、如果 sateph==EPHOPT_SBAS，则还要对 PC进行处理。之后给该函数计算出的伪距值的方差赋值。
    if (opt->sateph == EPHOPT_SBAS)
        PC -= P1_C1; /* sbas clock based C1 */

    *var = SQR(ERR_CBIAS);

    return PC;
}
/* ionospheric correction ------------------------------------------------------
 * compute ionospheric correction 计算给定电离层选项时的电离层延时(m)。
 * args   : gtime_t time     I   time
 *          nav_t  *nav      I   navigation data
 *          int    sat       I   satellite number
 *          double *pos      I   receiver position {lat,lon,h} (rad|m)
 *          double *azel     I   azimuth/elevation angle {az,el} (rad)
 *          int    ionoopt   I   ionospheric correction option (IONOOPT_???)
 *          double *ion      O   ionospheric delay (L1) (m)
 *          double *var      O   ionospheric delay (L1) variance (m^2)
 * return : status(1:ok,0:error)
 *-----------------------------------------------------------------------------*/
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var)
{
    //    trace(4,"ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
    //          time_str(time,3),ionoopt,sat,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
    //          azel[1]*R2D);

    /* broadcast model */
    //* 1\根据 opt的值，选用不同的电离层模型计算方法。
    //*     当 ionoopt==IONOOPT_BRDC时，调用 ionmodel，计算 Klobuchar模型时的电离层延时 (L1，m)；
    //*     当 ionoopt==IONOOPT_TEC时，调用 iontec，计算 TEC网格模型时的电离层延时 (L1，m)。


    //! 当 ionoopt==IONOOPT_IFLC时，此时通过此函数计算得到的延时和方差都为 0。
    //!     其实，对于 IFLC模型，其延时值在 prange函数中计算伪距时已经包括在里面了，
    //!     而方差是在 varerr函数中计算的，并且会作为导航系统误差的一部分给出。
    if (ionoopt == IONOOPT_BRDC)
    {
        *ion = ionmodel(time, nav->ion_gps, pos, azel);
        *var = SQR(*ion * ERR_BRDCI);
        return 1;
    }
    /* sbas ionosphere model */
    //    if (ionoopt==IONOOPT_SBAS) {
    //        return sbsioncorr(time,nav,pos,azel,ion,var);
    //    }
    /* ionex tec model */
    //    if (ionoopt==IONOOPT_TEC) {
    //        return iontec(time,nav,pos,azel,1,ion,var);
    //    }
    /* qzss broadcast model */
    if (ionoopt == IONOOPT_QZS && norm(nav->ion_qzs, 8) > 0.0)
    {
        *ion = ionmodel(time, nav->ion_qzs, pos, azel);
        *var = SQR(*ion * ERR_BRDCI);
        return 1;
    }
    /* lex ionosphere model */
    //    if (ionoopt==IONOOPT_LEX) {
    //        return lexioncorr(time,nav,pos,azel,ion,var);
    //    }
    *ion = 0.0;
    *var = ionoopt == IONOOPT_OFF ? SQR(ERR_ION) : 0.0;
    return 1;
}
/* tropospheric correction -----------------------------------------------------
 * compute tropospheric correction 计算对流层延时(m)。
 !貌似对流层延时与信号频率无关，所以这里计算得到的值并不是只针对于 L1信号！
 * args   : gtime_t time     I   time
 *          nav_t  *nav      I   navigation data
 *          double *pos      I   receiver position {lat,lon,h} (rad|m)
 *          double *azel     I   azimuth/elevation angle {az,el} (rad)
 *          int    tropopt   I   tropospheric correction option (TROPOPT_???)
 *          double *trp      O   tropospheric delay (m)
 *          double *var      O   tropospheric delay variance (m^2)
 * return : status(1:ok,0:error)
 *-----------------------------------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var)
{
    //    trace(4,"tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
    //          time_str(time,3),tropopt,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
    //          azel[1]*R2D);

    /* saastamoinen model */
    //* 1、当 tropopt==TROPOPT_SAAS或一些其它情况时，调用 tropmodel函数，计算 Saastamoinen模型下的对流层延时。
    if (tropopt == TROPOPT_SAAS || tropopt == TROPOPT_EST || tropopt == TROPOPT_ESTG)
    {
        *trp = tropmodel(time, pos, azel, REL_HUMI);
        *var = SQR(ERR_SAAS / (sin(azel[1]) + 0.1));
        return 1;
    }
    /* sbas troposphere model */
    //    if (tropopt==TROPOPT_SBAS) {
    //        *trp=sbstropcorr(time,pos,azel,var);
    //        return 1;
    //    }
    /* no correction */
    *trp = 0.0;
    *var = tropopt == TROPOPT_OFF ? SQR(ERR_TROP) : 0.0;
    return 1;
}
/* pseudorange residuals -----------------------------------------------------*/
/**
 * @brief 计算在当前接收机位置和钟差值的情况下，
 *          定位方程的右端部分 v(nv\*1)
 *          几何矩阵 H(NX*nv)
 *          此时所得的伪距残余的方差 var
 *          所有观测卫星的 azel{方位角、仰角}
 *          定位时有效性 vsat
 *          定位后伪距残差 resp
 *          参与定位的卫星个数 ns和方程个数 nv。
 *
 * 函数参数，17个
 * int      iter      I    迭代次数
 * obsd_t   *obs      I    observation data
 * int      n         I    number of observation data
 * double   *rs       I   satellite positions and velocities，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 * double   *dts      I   satellite clocks，长度为2*n， {bias,drift} (s|s/s)
 * double   *vare     I   sat position and clock error variances (m^2)
 * int      *svh      I   sat health flag (-1:correction not available)
 * nav_t    *nav      I   navigation data
 * double   *x        I   本次迭代开始之前的定位值
 * prcopt_t *opt      I   processing options
 * double   *v        O   定位方程的右端部分，伪距残余
 * double   *H        O   定位方程中的几何矩阵
 * double   *var      O   参与定位的伪距残余方差
 * double   *azel     O   对于当前定位值，每一颗观测卫星的 {方位角、高度角}
 * int      *vsat     O   每一颗观测卫星在当前定位时是否有效
 * double   *resp     O   每一颗观测卫星的伪距残余， (P-(r+c*dtr-c*dts+I+T))
 * int      *ns       O   参与定位的卫星的个数
 * 返回类型：
 * int                O   定位方程组的方程个数
 */
static int rescode(int iter, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   double *v, double *H, double *var, double *azel, int *vsat,
                   double *resp, int *ns)
{
    double r, dion, dtrp, vmeas, vion, vtrp, rr[3], pos[3], dtr, e[3], P, lam_L1;
    int i, j, nv = 0, sys, mask[4] = {0};

    // trace(3,"resprng : n=%d\n",n);
    //* 1、将之前得到的定位解信息赋值给rr和dtr数组，以进行关于当前解的伪距残余的相关计算。
    for (i = 0; i < 3; i++)
        rr[i] = x[i];
    dtr = x[3];

    //* 2、将得到的ecef位置信息转换为大地坐标系信息
    ecef2pos(rr, pos);

    for (i = *ns = 0; i < n && i < MAXOBS; i++)
    {
        //* 3、将 vsat、azel和 resp数组置 0，因为在前后两次定位结果中，每颗卫星的上述信息都会发生变化。
        vsat[i] = 0;
        azel[i * 2] = azel[1 + i * 2] = resp[i] = 0.0;
        //* 4、调用 satsys函数，验证卫星编号是否合理及其所属的导航系统。
        if (!(sys = satsys(obs[i].sat, NULL)))
            continue;

        /* reject duplicated observation data */
        //* 5、检测当前观测卫星是否和下一个相邻数据重复。是，则 i++后继续下一次循环；否，则进入下一步。
        if (i < n - 1 && i < MAXOBS - 1 && obs[i].sat == obs[i + 1].sat)
        {
            //            trace(2,"duplicated observation data %s sat=%2d\n",
            //                  time_str(obs[i].time,3),obs[i].sat);
            i++;
            continue;
        }
        /* geometric distance/azimuth/elevation angle */
        //* 6、计算卫星和当前接收机位置之间的几何距离和 receiver-to-satellite方向的单位向量。然后检验几何距离是否 >0。
        if ((r = geodist(rs + i * 6, rr, e)) <= 0.0 ||
            satazel(pos, e, azel + i * 2) < opt->elmin) //* 7、调用 satazel函数，计算在接收机位置处的站心坐标系中卫星的方位角和仰角，检验仰角是否≥截断值。
            continue;

        /* psudo range with code bias correction */
        //* 8、计算伪距值。
        if ((P = prange(obs + i, nav, azel + i * 2, iter, opt, &vmeas)) == 0.0)
            continue;

        /* excluded satellite? */
        //* 9、可以在处理选项中事先指定只选用哪些导航系统或卫星来进行定位，这是通过调用satexclude函数完成的。
        //        if (satexclude(obs[i].sat,svh[i],opt)) continue;

        /* ionospheric corrections */
        //* 10、计算电离层延时(m)已经转换为m的单位
        if (!ionocorr(obs[i].time, nav, obs[i].sat, pos, azel + i * 2,
                      iter > 0 ? opt->ionoopt : IONOOPT_BRDC, &dion, &vion))
            continue;

        /* GPS-L1 -> L1/B1 */
        //* 11、10中所得的电离层延时是建立在 L1信号上的，当使用其它频率信号时，
        //*     依据所用信号频组中第一个频率的波长与 L1波长的关系，对上一步得到的电离层延时进行修正。
        if ((lam_L1 = nav->lam[obs[i].sat - 1][0]) > 0.0)
        {
            dion *= SQR(lam_L1 / lam_carr[0]);
        }
        /* tropospheric corrections */
        //* 12、计算对流层延时(m)
        if (!tropcorr(obs[i].time, nav, pos, azel + i * 2,
                      iter > 0 ? opt->tropopt : TROPOPT_SAAS, &dtrp, &vtrp))
        {
            continue;
        }
        /* pseudorange residual */
        //* 13、计算伪距残差
        v[nv] = P - (r + dtr - CLIGHT * dts[i * 2] + dion + dtrp);

        /* design matrix */
        //* 14、组装几何矩阵，前 3行为 6中计算得到的视线单位向量的反向，第 4行为 1，其它行为 0。
        for (j = 0; j < NX; j++)
            H[j + nv * NX] = j < 3 ? -e[j] : (j == 3 ? 1.0 : 0.0);

        /* time system and receiver bias offset correction */
        if (sys == SYS_GLO)
        {
            v[nv] -= x[4];
            H[4 + nv * NX] = 1.0;
            mask[1] = 1;
        }
        else if (sys == SYS_GAL)
        {
            v[nv] -= x[5];
            H[5 + nv * NX] = 1.0;
            mask[2] = 1;
        }
        else if (sys == SYS_CMP)
        {
            v[nv] -= x[6];
            H[6 + nv * NX] = 1.0;
            mask[3] = 1;
        }
        else
            mask[0] = 1;

        //* 15、将参与定位的卫星的定位有效性标志设为1，给当前卫星的伪距残余赋值，参与定位的卫星个数 ns加 1.
        vsat[i] = 1;
        resp[i] = v[nv];
        (*ns)++;

        /* error variance */
        //* 16、调用 varerr函数，计算此时的导航系统误差（可能会包括 IFLC选项时的电离层延时），然后累加计算用户测距误差(URE)。
        var[nv++] = varerr(opt, azel[1 + i * 2], sys) + vare[i] + vmeas + vion + vtrp;

        //        trace(4,"sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n",obs[i].sat,
        //              azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],sqrt(var[nv-1]));
    }
    /* constraint to avoid rank-deficient */
    //* 17、为了防止亏秩，人为的添加了几组观测方程。
    for (i = 0; i < 4; i++)
    {
        if (mask[i])
            continue;
        v[nv] = 0.0;
        for (j = 0; j < NX; j++)
            H[j + nv * NX] = j == i + 3 ? 1.0 : 0.0;
        var[nv++] = 0.01;
    }
    return nv;
}
/* validate solution ---------------------------------------------------------*/
/**
 * @brief 确认当前解是否符合要求，即伪距残差小于某个χ^2值和GDOP小于某个门限值）
 * 
 * double   *azel     I   azimuth/elevation angle (rad)
 * int      *vsat     I   表征卫星在定位时是否有效
 * int      n         I   number of observation data
 * prcopt_t *opt      I   processing options
 * double   *v        I   定位后伪距残差 (P-(r+c*dtr-c*dts+I+T))
 * int      nv        I   定位方程的方程个数
 * int      nx        I   未知数的个数
 * char     *msg      O   error message for error exit
 * 返回类型：
 * int                O    (1:ok,0:error)
 */
static int valsol(const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, int nv, int nx,
                  char *msg)
{
    double azels[MAXOBS * 2], dop[4], vv;
    int i, ns;

    // trace(3,"valsol  : n=%d nv=%d\n",n,nv);

    /* chi-square validation of residuals */
    //* 1、计算定位后伪距残差平方加权和和vv
    vv = dot(v, v, nv);
    //* 2、卡方检测，定位误差是否过大
    if (nv > nx && vv > chisqr[nv - nx - 1])
    {
        sprintf(msg, "chi-square error nv=%d vv=%.1f cs=%.1f", nv, vv, chisqr[nv - nx - 1]);
        return 0;
    }
    /* large gdop check */
    //* 3、复制azel，这里只复制那些对于定位结果有贡献的卫星的 zael值，并且统计实现定位时所用卫星的数目。
    for (i = ns = 0; i < n; i++)
    {
        if (!vsat[i])
            continue;
        azels[ns * 2] = azel[i * 2];
        azels[1 + ns * 2] = azel[1 + i * 2];
        ns++;
    }
    //* 4、调用 dops函数，计算各种精度因子(DOP)，检验是否有 0<GDOP<max。
    //*     否，则说明该定位解的精度不符合要求，返回 0；是，则返回 1。
    dops(ns, azels, opt->elmin, dop);
    if (dop[0] <= 0.0 || dop[0] > opt->maxgdop)
    {
        sprintf(msg, "gdop error nv=%d gdop=%.1f", nv, dop[0]);
        return 0;
    }
    return 1;
}
/* estimate receiver position ------------------------------------------------*/
/**
 * @brief 通过伪距实现绝对定位，计算出接收机的位置和钟差，
 * 顺带返回实现定位后每颗卫星的{方位角、仰角}、定位时有效性、定位后伪距残差。
 *
 * obsd_t   *obs      I   observation data
 * int      n         I   number of observation data
 * double   *rs       I   satellite positions and velocities，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 * double   *dts      I   satellite clocks，长度为2*n， {bias,drift} (s|s/s)
 * double   *vare     I   sat position and clock error variances (m^2)
 * int      *svh      I   sat health flag (-1:correction not available)
 * nav_t    *nav      I   navigation data
 * prcopt_t *opt      I   processing options
 * sol_t    *sol      IO  solution
 * double   *azel     IO  azimuth/elevation angle (rad)
 * int      *vsat     IO  表征卫星在定位时是否有效
 * double   *resp     IO  定位后伪距残差 (P-(r+c*dtr-c*dts+I+T))
 * char     *msg      O   error message for error exit
 * 返回类型:
 * int                O     (1:ok,0:error)
 */
static int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, int *vsat,
                  double *resp, char *msg)
{
    double x[NX] = {0}, dx[NX], Q[NX * NX], sig;
    int i, j, k, info, stat, nv, ns;

    // trace(3,"estpos  : n=%d\n",n);
    double H[(MAXOBS + 4)], v[7 * (MAXOBS + 4)], var[MAXOBS + 4];
    //    v=mat(n+4,1); H=mat(NX,n+4); var=mat(n+4,1);

    //* 1、将 sol->rr的前 3项赋值给 x数组
    for (i = 0; i < 3; i++)
        x[i] = sol->rr[i];

    //* 2、开始迭代，最大迭代10次，用的是梯度下降方法求解最小二乘
    for (i = 0; i < MAXITR; i++)
    {

        /* pseudo range residuals */
        //* 首先调用 rescode函数，计算在当前接收机位置和钟差值的情况下，
        //*     定位方程的右端部分 v(nv\*1)、几何矩阵 H(NX*nv)、
        //*     此时所得的
        //*             伪距残余的方差 var
        //*             所有观测卫星的 azel{方位角、仰角}
        //*             定位时有效性 vsat
        //*             定位后伪距残差 resp
        //*             参与定位的卫星个数 ns和方程个数 nv。
        nv = rescode(i, obs, n, rs, dts, vare, svh, nav, x, opt, v, H, var, azel, vsat, resp,
                     &ns);
        //* 3、确定方程组中方程的个数要大于未知数的个数。
        if (nv < NX)
        {
            sprintf(msg, "lack of valid sats ns=%d", nv);
            break;
        }
        /* weight by variance */
        //* 4、以伪距残余的标准差的倒数作为权重，
        //*     对 H和 v分别左乘权重对角阵，得到加权之后的 H和 v。
        for (j = 0; j < nv; j++)
        {
            sig = sqrt(var[j]);
            v[j] /= sig;
            for (k = 0; k < NX; k++)
                H[k + j * NX] /= sig;
        }
        /* least square estimation */
        //* 5、调用 lsq函数，根据 Δx = (HH')^-1Hv和Q = (HH')^-1，
        //*     得到当前 x的修改量和定位误差协方差矩阵中的权系数阵。
        if ((info = lsq(H, v, NX, nv, dx, Q)))
        {
            sprintf(msg, "lsq error info=%d", info);
            break;
        }
        //* 6、将5中求得的x加入到当前x值中，得到更新之后的x值。
        for (j = 0; j < NX; j++)
            x[j] += dx[j];

        //* 7、如果 5中求得的修改量小于截断因子(目前是1e-4)，
        //*     则将 6中得到的 x值作为最终的定位结果，对 sol的相应参数赋值，
        //*     之后再调用 valsol函数确认当前解是否符合要求（伪距残余小于某个 值和 GDOP小于某个门限值）。
        //*     否则，进行下一次循环。
        if (norm(dx, NX) < 1E-4)
        {
            sol->type = 0;
            sol->time = timeadd(obs[0].time, -x[3] / CLIGHT);
            sol->dtr[0] = x[3] / CLIGHT; /* receiver clock bias (s) */
            sol->dtr[1] = x[4] / CLIGHT; /* glo-gps time offset (s) */
            sol->dtr[2] = x[5] / CLIGHT; /* gal-gps time offset (s) */
            sol->dtr[3] = x[6] / CLIGHT; /* bds-gps time offset (s) */
            for (j = 0; j < 6; j++)
                sol->rr[j] = j < 3 ? x[j] : 0.0;
            for (j = 0; j < 3; j++)
                sol->qr[j] = (float)Q[j + j * NX];
            sol->qr[3] = (float)Q[1];      /* cov xy */
            sol->qr[4] = (float)Q[2 + NX]; /* cov yz */
            sol->qr[5] = (float)Q[2];      /* cov zx */
            sol->ns = (unsigned char)ns;
            sol->age = sol->ratio = 0.0;

            /* validate solution */
            //! 如果某次迭代过程中步长小于门限值(1e-4)，但经 valsol函数检验后该解无效，
            //!     则会直接返回 0，并不会再进行下一次迭代计算。
            if ((stat = valsol(azel, vsat, n, opt, v, nv, NX, msg)))
            {
                sol->stat = opt->sateph == EPHOPT_SBAS ? SOLQ_SBAS : SOLQ_SINGLE;
            }
            //            free(v); free(H); free(var);

            return stat;
        }
    }
    //* 8、如果超过了规定的循环次数，则输出发散信息后，返回 0。
    if (i >= MAXITR)
        sprintf(msg, "iteration divergent i=%d", i);

    // free(v); free(H); free(var);

    return 0;
}
/* raim fde (failure detection and exclution) -------------------------------*/
/**
 * 函数参数，13个：
 * obsd_t   *obs      I   observation data
 * int      n         I   number of observation data
 * double   *rs       I   satellite positions and velocities，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 * double   *dts      I   satellite clocks，长度为2*n， {bias,drift} (s|s/s)
 * double   *vare     I   sat position and clock error variances (m^2)
 * int      *svh      I   sat health flag (-1:correction not available)
 * nav_t    *nav      I   navigation data
 * prcopt_t *opt      I   processing options
 * sol_t    *sol      IO  solution
 * double   *azel     IO  azimuth/elevation angle (rad)
 * int      *vsat     IO  表征卫星在定位时是否有效
 * double   *resp     IO  定位后伪距残差 (P-(r+c*dtr-c*dts+I+T))
 * char     *msg      O   error message for error exit
 * 返回类型:
 * int                O     (1:ok,0:error)
 */
static int raim_fde(const obsd_t *obs, int n, const double *rs,
                    const double *dts, const double *vare, const int *svh,
                    const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                    double *azel, int *vsat, double *resp, char *msg)
{
    obsd_t obs_e[sizeof(obsd_t) * MAXOBS];
    sol_t sol_e = {{0}};
    char tstr[32], name[16], msg_e[128];
    double rs_e[6 * MAXOBS], dts_e[2 * MAXOBS], vare_e[1 * MAXOBS], azel_e[2 * MAXOBS], resp_e[1 * MAXOBS], rms_e, rms = 100.0;
    int i, j, k, nvsat, stat = 0, svh_e[1 * MAXOBS], vsat_e[1 * MAXOBS], sat = 0;

    // trace(3,"raim_fde: %s n=%2d\n",time_str(obs[0].time,0),n);
    //    if (!(obs_e=(obsd_t *)malloc(sizeof(obsd_t)*n))) return 0;
    //    rs_e = mat(6,n); dts_e = mat(2,n); vare_e=mat(1,n); azel_e=zeros(2,n);
    //    svh_e=imat(1,n); vsat_e=imat(1,n); resp_e=mat(1,n);
    //* 1、关于观测卫星数目的循环，每次舍弃一颗卫星，计算使用余下卫星进行定位的定位值。
    for (i = 0; i < n; i++)
    {

        /* satellite exclution */
        for (j = k = 0; j < n; j++)
        {
            if (j == i)
                continue;
            obs_e[k] = obs[j];
            matcpy(rs_e + 6 * k, rs + 6 * j, 6, 1);
            matcpy(dts_e + 2 * k, dts + 2 * j, 2, 1);
            vare_e[k] = vare[j];
            svh_e[k++] = svh[j];
        }
        /* estimate receiver position without a satellite */
        //* 2、舍弃一颗卫星后，将剩下卫星的数据复制到一起，调用 estpos函数，计算使用余下卫星进行定位的定位值。
        if (!estpos(obs_e, n - 1, rs_e, dts_e, vare_e, svh_e, nav, opt, &sol_e, azel_e,
                    vsat_e, resp_e, msg_e))
        {
            // trace(3,"raim_fde: exsat=%2d (%s)\n",obs[i].sat,msg);
            continue;
        }
        //* 3、累加使用当前卫星实现定位后的伪距残差平方和与可用微信数目
        for (j = nvsat = 0, rms_e = 0.0; j < n - 1; j++)
        {
            if (!vsat_e[j])
                continue;
            rms_e += SQR(resp_e[j]);
            nvsat++;
        }
        //* 如果 nvsat<5，则说明当前卫星数目过少，无法进行 RAIM_FDE操作。
        if (nvsat < 5)
        {
            //            trace(3,"raim_fde: exsat=%2d lack of satellites nvsat=%2d\n",
            //                  obs[i].sat,nvsat);
            continue;
        }
        //* 4、计算伪距残差平方和的标准平均值
        rms_e = sqrt(rms_e / nvsat);

        //        trace(3,"raim_fde: exsat=%2d rms=%8.3f\n",obs[i].sat,rms_e);

        if (rms_e > rms)
            continue;

        /* save result */
        //*     如果小于 rms，则说明当前定位结果更合理，
        //*     将 stat置为 1，重新更新 sol、azel、vsat(当前被舍弃的卫星，此值置为0)、resp等值
        //*     并将当前的 rms_e更新到 `rms'中。
        for (j = k = 0; j < n; j++)
        {
            if (j == i)
                continue;
            matcpy(azel + 2 * j, azel_e + 2 * k, 2, 1);
            vsat[j] = vsat_e[k];
            resp[j] = resp_e[k++];
        }
        stat = 1;
        *sol = sol_e;
        sat = obs[i].sat;
        rms = rms_e;
        vsat[i] = 0;
        strcpy(msg, msg_e);
    }
    //* 6、如果 stat不为 0，则说明在弃用卫星的前提下有更好的解出现，输出信息，指出弃用了哪科卫星。
    if (stat)
    {
        time2str(obs[0].time, tstr, 2);
        satno2id(sat, name);
        //        trace(2,"%s: %s excluded by raim\n",tstr+11,name);
    }

    //* 5、继续弃用下一颗卫星，重复 2-4操作。总而言之，将同样是弃用一颗卫星条件下，
    //*     伪距残差标准平均值最小的组合所得的结果作为最终的结果输出。
    //    free(obs_e);
    //    free(rs_e ); free(dts_e ); free(vare_e); free(azel_e);
    //    free(svh_e); free(vsat_e); free(resp_e);
    return stat;
}
/* doppler residuals ---------------------------------------------------------*/
/**
 * @brief 计算定速方程组左边的几何矩阵和右端的速度残余，返回定速时所使用的卫星数目
 * 函数参数，11个：
 * obsd_t   *obs      I   observation data
 * int      n         I   number of observation data
 * double   *rs       I   satellite positions and velocities，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 * double   *dts      I   satellite clocks，长度为2*n， {bias,drift} (s|s/s)
 * nav_t    *nav      I   navigation data
 * double   *rr       I   receiver positions and velocities，长度为6，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 * double   *x        I   本次迭代开始之前的定速值
 * double   *azel     I   azimuth/elevation angle (rad)
 * int      *vsat     I   表征卫星在定速时是否有效
 * double   *v        O   定速方程的右端部分，速度残余
 * double   *H        O   定速方程中的几何矩阵
 * 返回类型:
 * int                O    定速时所使用的卫星数目
 */
static int resdop(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const nav_t *nav, const double *rr, const double *x,
                  const double *azel, const int *vsat, double *v, double *H)
{
    double lam, rate, pos[3], E[9], a[3], e[3], vs[3], cosel;
    int i, j, nv = 0;

    //    trace(3,"resdop  : n=%d\n",n);
    //* 1、调用 ecef2pos函数，将接收机位置由 ECEF转换为大地坐标系
    ecef2pos(rr, pos);
    //* 2、调用 xyz2enu函数，计算此时的坐标转换矩阵。
    xyz2enu(pos, E);
    for (i = 0; i < n && i < MAXOBS; i++)
    {

        lam = nav->lam[obs[i].sat - 1][0];
        //* 3、当前卫星在定速时不可用，直接进行下一次循环。
        if (obs[i].D[0] == 0.0 || lam == 0.0 || !vsat[i] || norm(rs + 3 + i * 6, 3) <= 0.0)
        {
            continue;
        }
        //* 4、计算当前接收机位置下 ENU中的视向量，然后转换得到 ECEF中视向量的值。
        /* line-of-sight vector in ecef */
        cosel = cos(azel[1 + i * 2]);
        a[0] = sin(azel[i * 2]) * cosel;
        a[1] = cos(azel[i * 2]) * cosel;
        a[2] = sin(azel[1 + i * 2]);
        matmul("TN", 3, 1, 3, 1.0, E, a, 0.0, e);

        /* satellite velocity relative to receiver in ecef */
        //* 5、计算 ECEF中卫星相对于接收机的速度，然后再计算出考虑了地球自转的用户和卫星之间的几何距离变化率，校正公式见 RTKLIB manual P159 (E.6.29)
        for (j = 0; j < 3; j++)
            vs[j] = rs[j + 3 + i * 6] - x[j];

        /* range rate with earth rotation correction */
        //* 6、根据公式 计算出方程右端项的多普勒残余，然后再构建左端项的几何矩阵。最后再将观测方程数增 1.
        rate = dot(vs, e, 3) + OMGE / CLIGHT * (rs[4 + i * 6] * rr[0] + rs[1 + i * 6] * x[0] - rs[3 + i * 6] * rr[1] - rs[i * 6] * x[1]);

        /* doppler residual */
        v[nv] = -lam * obs[i].D[0] - (rate + x[3] - CLIGHT * dts[1 + i * 2]);

        /* design matrix */
        for (j = 0; j < 4; j++)
            H[j + nv * 4] = j < 3 ? -e[j] : 1.0;

        nv++;
    }
    return nv;
}
/* estimate receiver velocity ------------------------------------------------*/
/**
 * obsd_t   *obs      I   observation data
 * int      n         I   number of observation data
 * double   *rs       I   satellite positions and velocities，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 * double   *dts      I   satellite clocks，长度为2*n， {bias,drift} (s|s/s)
 * nav_t    *nav      I   navigation data
 * prcopt_t *opt      I   processing options
 * sol_t    *sol      IO  solution
 * double   *azel     IO  azimuth/elevation angle (rad)
 * int      *vsat     IO  表征卫星在定位时是否有效
 * 返回类型:
 * int                O     (1:ok,0:error)
 * 不像定位时，初值为上一历元的位置，定速直接给的0
 */
static void estvel(const obsd_t *obs, int n, const double *rs, const double *dts,
                   const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                   const double *azel, const int *vsat)
{
    double x[4] = {0}, dx[4], Q[16], v[MAXOBS], H[4 * MAXOBS];
    int i, j, nv;

    //    trace(3,"estvel  : n=%d\n",n);

    //    v=mat(n,1); H=mat(4,n);
    //* 1、在最大迭代次数限制内，调用resdop，计算定速方程组左边的几何矩阵和右端的速度残余，返回定速时所使用的卫星数目。
    for (i = 0; i < MAXITR; i++)
    {

        /* doppler residuals */
        if ((nv = resdop(obs, n, rs, dts, nav, sol->rr, x, azel, vsat, v, H)) < 4)
        {
            break;
        }
        /* least square estimation */
        //* 2、调用 lsq函数，解出 {速度、频漂}的步长，累加到 x中。
        if (lsq(H, v, 4, nv, dx, Q))
            break;

        for (j = 0; j < 4; j++)
            x[j] += dx[j];
        //* 3、检查当前计算出的步长的绝对值是否小于 1E-6。是，则说明当前解已经很接近真实值了，
        //*     将接收机三个方向上的速度存入到 sol_t.rr中；否，则进行下一次循环。
        if (norm(dx, 4) < 1E-6)
        {
            for (i = 0; i < 3; i++)
                sol->rr[i + 3] = x[i]; //! 注意，最后并没有存储计算出的接收器的时钟漂移
            break;
        }
    }
    //    free(v); free(H);
}
/* single-point positioning ----------------------------------------------------
 * compute receiver position, velocity, clock bias by single-point positioning
 * with pseudorange and doppler observables
 * 依靠多普勒频移测量值和伪距来进行单点定位，给出接收机的位置、速度和钟差
 * args   : obsd_t *obs      I   observation data
 *          int    n         I   number of observation data
 *          nav_t  *nav      I   navigation data
 *          prcopt_t *opt    I   processing options
 *          sol_t  *sol      IO  solution
 *          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
 *          ssat_t *ssat     IO  satellite status              (NULL: no output)
 *          char   *msg      O   error message for error exit
 * return : status(1:ok,0:error)
 * notes  : assuming sbas-gps, galileo-gps, qzss-gps, compass-gps time offset and
 *          receiver bias are negligible (only involving glonass-gps time offset
 *          and receiver bias)
 *-----------------------------------------------------------------------------*/
extern int pntpos(const obsd_t *obs, int n, nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, ssat_t *ssat,
                  char *msg)
{
    prcopt_t opt_ = *opt;
    double rs[6 * MAXOBS], dts[2 * MAXOBS], var[MAXOBS], azel_[2 * MAXOBS], resp[MAXOBS];
    //    double *rs,*dts,*var,*azel_,*resp;
    int i, stat, vsat[MAXOBS] = {0}, svh[MAXOBS];

    sol->stat = SOLQ_NONE;
    //* 1、检查卫星个数是否>0
    if (n <= 0)
    {
        strcpy(msg, "no observation data");
        return 0;
    }

    // trace(3,"pntpos  : tobs=%s n=%d\n",time_str(obs[0].time,3),n);

    sol->time = obs[0].time;
    msg[0] = '\0';

    //    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel_=zeros(2,n); resp=mat(1,n);
    //* 2、当模式不是单点模式时
    if (opt_.mode != PMODE_SINGLE)
    { /* for precise positioning */
#if 0
        opt_.sateph =EPHOPT_BRDC;
#endif
        opt_.ionoopt = IONOOPT_BRDC; //*电离层矫正模型
        opt_.tropopt = TROPOPT_SAAS; //*对流层矫正模型
    }
    /* satellite positons, velocities and clocks */
    //* 3、按照所观测到的卫星顺序计算出没课卫星的位置、速度、（钟差，频漂）
    satposs(sol->time, obs, n, nav, opt_.sateph, rs, dts, var, svh);

    /* estimate receiver position with pseudorange */
    //* 4、通过伪距实现绝对定位，计算出接收机的位置和钟差，顺带返回实现定位后每颗卫星的(\
    //*     方位角，仰角)、定位时有效性、定位后的伪距残差
    stat = estpos(obs, n, rs, dts, var, svh, nav, &opt_, sol, azel_, vsat, resp, msg);

    //* 5、对上一步得到的定位结果进行接收机自主正直性检测（RAIM）。通过再次使用 vsat数组，
    //*     这里只会在对定位结果有贡献的卫星数据进行检测。
    /* raim fde */
    if (!stat && n >= 6 && opt->posopt[4])
    {
        // if (!stat&&n>=6) {
        stat = raim_fde(obs, n, rs, dts, var, svh, nav, &opt_, sol, azel_, vsat, resp, msg);
    }
    /* estimate receiver velocity with doppler */
    //* 6、 调用 estvel函数，依靠多普勒频移测量值计算接收机的速度。
    //*     这里只使用通过了上一步RAIM_FDE操作的卫星数据，所以对于计算出的速度就没有再次进行 RAIM了。
    //* 这里只计算了接收机的钟差，而没有计算接收机的频漂，
    //*     原因在于 estvel函数中虽然计算得到了接收机频漂，但并没有将其输出到 sol_t:dtr中。
    if (stat)
        estvel(obs, n, rs, dts, nav, &opt_, sol, azel_, vsat);

    if (azel)
    {
        for (i = 0; i < n * 2; i++)
            azel[i] = azel_[i];
    }
    if (ssat)
    {
        //* 7、首先将 ssat_t结构体数组的
        //*     vs(定位时有效性)
        //*     azel（方位角、仰角）
        //*     resp(伪距残余)
        //*     resc(载波相位残余)
        //*     snr(信号强度)           都置为 0，
        for (i = 0; i < MAXSAT; i++)
        {
            ssat[i].vs = 0;
            ssat[i].azel[0] = ssat[i].azel[1] = 0.0;
            ssat[i].resp[0] = ssat[i].resc[0] = 0.0;
            ssat[i].snr[0] = 0;
        }
        //* 将实现定位后的 azel、snr赋予 ssat_t结构体数组，
        //*     而 vs、resp则只赋值给那些对定位有贡献的卫星，没有参与定位的卫星，这两个属性值为 0。
        for (i = 0; i < n; i++)
        {
            ssat[obs[i].sat - 1].azel[0] = azel_[i * 2];
            ssat[obs[i].sat - 1].azel[1] = azel_[1 + i * 2];
            ssat[obs[i].sat - 1].snr[0] = obs[i].SNR[0];
            if (!vsat[i])
                continue;
            ssat[obs[i].sat - 1].vs = 1;
            ssat[obs[i].sat - 1].resp[0] = resp[i];
        }
    }
    // 修改后只使用了局部变量没有使用malloc，因此不用free
    //    free(rs); free(dts); free(var); free(azel_); free(resp);
    return stat;
}
