#ifndef __COMPLEX_H__
#define __COMPLEX_H__
/*3:*/
#line 66 "complex.w"
struct complex{double re,im;};/*7:*/
#line 129 "complex.w"
struct complex cset(double a,double b)/*:7*/
#line 68 "complex.w"
;/*9:*/
#line 144 "complex.w"
struct complex cpolarset(double r,double theta)/*:9*/
#line 69 "complex.w"
;/*11:*/
#line 157 "complex.w"
double cabbs(struct complex z)/*:11*/
#line 71 "complex.w"
;/*15:*/
#line 190 "complex.w"
double carg(struct complex z)/*:15*/
#line 72 "complex.w"
;/*21:*/
#line 235 "complex.w"
struct complex csqr(struct complex z)/*:21*/
#line 73 "complex.w"
;/*13:*/
#line 181 "complex.w"
struct complex conj(struct complex z)/*:13*/
#line 74 "complex.w"
;/*17:*/
#line 201 "complex.w"
double cnorm(struct complex z)/*:17*/
#line 75 "complex.w"
;/*19:*/
#line 210 "complex.w"
struct complex csqrt(struct complex z)/*:19*/
#line 76 "complex.w"
;/*23:*/
#line 246 "complex.w"
struct complex cinv(struct complex w)/*:23*/
#line 77 "complex.w"
;/*26:*/
#line 271 "complex.w"
struct complex cadd(struct complex z,struct complex w)/*:26*/
#line 79 "complex.w"
;/*28:*/
#line 286 "complex.w"
struct complex csub(struct complex z,struct complex w)/*:28*/
#line 80 "complex.w"
;/*30:*/
#line 300 "complex.w"
struct complex cmul(struct complex z,struct complex w)/*:30*/
#line 81 "complex.w"
;/*32:*/
#line 315 "complex.w"
struct complex cdiv(struct complex z,struct complex w)/*:32*/
#line 82 "complex.w"
;/*34:*/
#line 343 "complex.w"
double crdiv(struct complex z,struct complex w)/*:34*/
#line 84 "complex.w"
;/*36:*/
#line 367 "complex.w"
double crmul(struct complex z,struct complex w)/*:36*/
#line 85 "complex.w"
;/*41:*/
#line 394 "complex.w"
struct complex csadd(double x,struct complex z)/*:41*/
#line 86 "complex.w"
;/*43:*/
#line 409 "complex.w"
struct complex csdiv(double x,struct complex w)/*:43*/
#line 87 "complex.w"
;/*39:*/
#line 380 "complex.w"
struct complex csmul(double x,struct complex z)/*:39*/
#line 88 "complex.w"
;/*46:*/
#line 438 "complex.w"
struct complex csin(struct complex z)/*:46*/
#line 90 "complex.w"
;/*48:*/
#line 449 "complex.w"
struct complex ccos(struct complex z)/*:48*/
#line 91 "complex.w"
;/*50:*/
#line 480 "complex.w"
struct complex ctan(struct complex z)/*:50*/
#line 92 "complex.w"
;/*52:*/
#line 502 "complex.w"
struct complex casin(struct complex z)/*:52*/
#line 93 "complex.w"
;/*54:*/
#line 515 "complex.w"
struct complex cacos(struct complex z)/*:54*/
#line 94 "complex.w"
;/*56:*/
#line 528 "complex.w"
struct complex catan(struct complex z)/*:56*/
#line 95 "complex.w"
;/*61:*/
#line 550 "complex.w"
struct complex csinh(struct complex z)/*:61*/
#line 97 "complex.w"
;/*59:*/
#line 541 "complex.w"
struct complex ccosh(struct complex z)/*:59*/
#line 98 "complex.w"
;/*63:*/
#line 559 "complex.w"
struct complex ctanh(struct complex z)/*:63*/
#line 99 "complex.w"
;/*65:*/
#line 572 "complex.w"
struct complex catanh(struct complex z)/*:65*/
#line 100 "complex.w"
;/*67:*/
#line 581 "complex.w"
struct complex casinh(struct complex z)/*:67*/
#line 101 "complex.w"
;/*70:*/
#line 592 "complex.w"
struct complex cexp(struct complex z)/*:70*/
#line 103 "complex.w"
;/*72:*/
#line 602 "complex.w"
struct complex clog(struct complex z)/*:72*/
#line 104 "complex.w"
;/*74:*/
#line 611 "complex.w"
struct complex clog10(struct complex z)/*:74*/
#line 105 "complex.w"
;/*77:*/
#line 624 "complex.w"
struct complex*new_carray(long size)/*:77*/
#line 107 "complex.w"
;/*79:*/
#line 640 "complex.w"
void free_carray(struct complex*a)/*:79*/
#line 108 "complex.w"
;/*81:*/
#line 652 "complex.w"
struct complex*copy_carray(struct complex*a,long size)/*:81*/
#line 109 "complex.w"
;/*83:*/
#line 669 "complex.w"
void set_carray(struct complex*a,long size,struct complex z)/*:83*/
#line 110 "complex.w"
;/*:3*/
#endif