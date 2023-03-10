#ifndef __MIE_H__
#define __MIE_H__
/*3:*//*10:*/
#line 120 "mie.w"
struct complex Lentz_Dn(struct complex z,long n)/*:10*/
#line 48 "mie.w"
;/*19:*/
#line 241 "mie.w"
void Dn_down(struct complex z,long nstop,struct complex*D)/*:19*/
#line 49 "mie.w"
;/*16:*/
#line 210 "mie.w"
void Dn_up(struct complex z,long nstop,struct complex*D)/*:16*/
#line 50 "mie.w"
;/*22:*/
#line 274 "mie.w"
void small_Mie(double x,struct complex m,double*mu,
long nangles,struct complex*s1,
struct complex*s2,double*qext,double*qsca,
double*qback,double*g)/*:22*/
#line 51 "mie.w"
;/*31:*/
#line 455 "mie.w"
void Mie(double x,struct complex m,double*mu,long nangles,struct complex*s1,
struct complex*s2,double*qext,double*qsca,double*qback,double*g)/*:31*/
#line 52 "mie.w"
;/*47:*/
#line 751 "mie.w"
void ez_Mie(double x,double n,double*qsca,double*g)/*:47*/
#line 53 "mie.w"
;/*:3*/
#endif