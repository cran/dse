
/*

NOTE  This C code is out of date. The fortran version of the code has been
       modified and this file has not been updated!!!!!!!!
       
  -- translated by f2c (version 20010224).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__48 = 48;
static integer c__1 = 1;
static integer c__100 = 100;
static integer c__16 = 16;
static integer c__40 = 40;
static integer c__37 = 37;
static doublereal c_b236 = 4.;
static integer c__11 = 11;
static integer c__22 = 22;

/*           Copyright 1993, 1994, 1995, 1996  Bank of Canada. */
/*           Copyright 1996, 1997  Paul Gilbert. */
/*           Copyright 1998, 2001  Bank of Canada. */
/*  Any person using this software has the right to use, */
/*        reproduce and distribute it. */
/*  The Bank does not warrant the accuracy of the information contained in the */
/*  software. User assumes all liability for any loss or damages, either direct */
/*  or indirect, arising from the use of the software. */
/* ------------------------------------------------------------------- */
/*  A C code version of this code is also distribute. It has been generated */
/*   using the following (extracted from the f2c readme). */
/* NOTE: You may exercise f2c by sending netlib@netlib.bell-labs.com */
/*  a message whose first line is "execute f2c" and whose remaining */
/*  lines are the Fortran 77 source that you wish to have converted. */
/*  Return mail brings you the resulting C, with f2c's error */
/*  messages between #ifdef uNdEfInEd and #endif at the end. */
/* Compile with: f77 -c -o dsefor.Sun4.o dsefor.f */
/*           or */
/*               f77 -c -o dsefor.Sun5.o dsefor.f */
/*           or */
/*             /opt/SUNWspro/f77 -c -o dsefor.S3.3.Sun5.o dsefor.f */
/*           or */
/*               f77.SC3.01.Sun4       -c -o dsefor.S3.3.Sun4.o dsefor.f */
/*     (Splus 3.3 requires compiler SC3.0.1 not SC3.0 as for Splus 3.1 and 3.2) */
/*           or  (with Solaris f77  (SunPro) )   for R */
/*               f77  -G -pic -o dse.so dsefor.f */
/*   or */
/*               f77  -G -pic -ansi -o dse.so dsefor.f */
/*   or */
/*        f77 -G -pic -ansi -L/home/asd3/opt/SUNWspro.old/lib -lF77 -lM77 -lm */
/*               -lc -lucb -o dse.so dsefor.f */
/*   or */
/*               f77  -G -pic -fast -o dse.so dsefor.f   bombs */
/*   or */
/*               f77  -c -pic -o dsefor.o dsefor.f */
/*               ld -G -o dse.so dsefor.o */
/*                  -pic:  Generate pic code with short offset */
/*                  -c:    Produce '.o' file. Do not run ld. */
/*                  -P:    Generate optimized code */
/*                  -fast:  Bundled set of options for best performance */
/*                  -G:     Create shared object */
/*                  -O1:    Generate optimized code */
/*                  -O2:    Generate optimized code */
/*                  -O3:    Generate optimized code */
/*                  -O4:    Generate optimized code */
/*                  -O:     Generate optimized code */
/*                  Suffix 'so':    Shared object */
/*           or  (the following seems to be the preferred method for Splus) */
/*               splus COMPILE dsefor.f */
/*                 and then mv dsefor.o to dsefor.Sunx.o */

/*   The pararameter IS=xxx controls the maximum size of the state in KF models. */
/*   A larger state makes S take more memory. */
/* Compile with: f77 -c -o dsefor.Sun4.large.o dsefor.f */
/*           or */
/*               f77 -c -o dsefor.Sun5.large.o dsefor.f */


/*        1         2         3         4         5         6         712       8 */
/* Subroutine */ int error_(str, l, is, n, str_len)
char *str;
integer *l, *is, *n;
ftnlen str_len;
{
/*  STR is a string of length L and IS is an integer vector of length N */
/*      CALL INTPR(STR,L,IS,N) */
/*      WRITE(STR,IS) */
    /* Parameter adjustments */
    --str;
    --is;

    /* Function Body */
    return 0;
} /* error_ */

/* Subroutine */ int dbpr_(str, l, is, n, str_len)
char *str;
integer *l, *is, *n;
ftnlen str_len;
{
/*  STR is a string of length L and IS is an integer vector of length N */
/*      CALL INTPR(STR,L,IS,N) */
/*      WRITE(STR,IS) */
    /* Parameter adjustments */
    --str;
    --is;

    /* Function Body */
    return 0;
} /* dbpr_ */

/* Subroutine */ int dbprdb_(str, l, r__, n, str_len)
char *str;
integer *l;
doublereal *r__;
integer *n;
ftnlen str_len;
{
/*  STR is a string of length L and R is an double real vector of length N */
/*      CALL DBLEPR(STR,L, R, N) */
/*      WRITE(STR,R) */
    /* Parameter adjustments */
    --str;
    --r__;

    /* Function Body */
    return 0;
} /* dbprdb_ */

/* Subroutine */ int simss_(y, z__, m, n, p, nsmpl, u, w, e, f, g, h__, fk, q,
	 r__, gain)
doublereal *y, *z__;
integer *m, *n, *p, *nsmpl;
doublereal *u, *w, *e, *f, *g, *h__, *fk, *q, *r__;
integer *gain;
{
    /* System generated locals */
    integer z_dim1, z_offset, y_dim1, y_offset, u_dim1, u_offset, w_dim1, 
	    w_offset, e_dim1, e_offset, f_dim1, f_offset, g_dim1, g_offset, 
	    h_dim1, h_offset, r_dim1, r_offset, q_dim1, q_offset, fk_dim1, 
	    fk_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, it;


/*  Simulate a state space model model. */
/*  See s code simulate.ss for details */

/* Note: the first period is done in the calling S routine so that */
/*   intial conditions can be handled there. */

    /* Parameter adjustments */
    q_dim1 = *n;
    q_offset = 1 + q_dim1 * 1;
    q -= q_offset;
    g_dim1 = *n;
    g_offset = 1 + g_dim1 * 1;
    g -= g_offset;
    f_dim1 = *n;
    f_offset = 1 + f_dim1 * 1;
    f -= f_offset;
    r_dim1 = *p;
    r_offset = 1 + r_dim1 * 1;
    r__ -= r_offset;
    fk_dim1 = *n;
    fk_offset = 1 + fk_dim1 * 1;
    fk -= fk_offset;
    h_dim1 = *p;
    h_offset = 1 + h_dim1 * 1;
    h__ -= h_offset;
    e_dim1 = *nsmpl;
    e_offset = 1 + e_dim1 * 1;
    e -= e_offset;
    w_dim1 = *nsmpl;
    w_offset = 1 + w_dim1 * 1;
    w -= w_offset;
    u_dim1 = *nsmpl;
    u_offset = 1 + u_dim1 * 1;
    u -= u_offset;
    z_dim1 = *nsmpl;
    z_offset = 1 + z_dim1 * 1;
    z__ -= z_offset;
    y_dim1 = *nsmpl;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;

    /* Function Body */
    i__1 = *nsmpl;
    for (it = 2; it <= i__1; ++it) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[it + i__ * z_dim1] = (float)0.;
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
/* L2001: */
		z__[it + i__ * z_dim1] += f[i__ + k * f_dim1] * z__[it - 1 + 
			k * z_dim1];
	    }
	}
	if (*m != 0) {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__2 = *m;
		for (k = 1; k <= i__2; ++k) {
/* L2003: */
		    z__[it + i__ * z_dim1] += g[i__ + k * g_dim1] * u[it + k *
			     u_dim1];
		}
	    }
	}
	if (*gain == 1) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__3 = *p;
		for (k = 1; k <= i__3; ++k) {
/* L2004: */
		    z__[it + i__ * z_dim1] += fk[i__ + k * fk_dim1] * w[it - 
			    1 + k * w_dim1];
		}
	    }
	} else {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__2 = *n;
		for (k = 1; k <= i__2; ++k) {
/* L2005: */
		    z__[it + i__ * z_dim1] += q[i__ + k * q_dim1] * e[it - 1 
			    + k * e_dim1];
		}
	    }
	}
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[it + i__ * y_dim1] = (float)0.;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
/* L2010: */
		y[it + i__ * y_dim1] += h__[i__ + j * h_dim1] * z__[it + j * 
			z_dim1];
	    }
	}
	if (*gain == 1) {
	    i__3 = *p;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L2014: */
		y[it + i__ * y_dim1] += w[it + i__ * w_dim1];
	    }
	} else {
	    i__3 = *p;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__2 = *p;
		for (k = 1; k <= i__2; ++k) {
/* L2015: */
		    y[it + i__ * y_dim1] += r__[i__ + k * r_dim1] * e[it + k *
			     e_dim1];
		}
	    }
	}
/* L1000: */
    }
    return 0;
} /* simss_ */


/* Subroutine */ int smooth_(z__, trkerr, u, y, n, m, p, nsmpl, f, g, h__, rr,
	 a, d__, l, pt1, zt)
doublereal *z__, *trkerr, *u, *y;
integer *n, *m, *p, *nsmpl;
doublereal *f, *g, *h__, *rr, *a, *d__, *l, *pt1, *zt;
{
    /* System generated locals */
    integer z_dim1, z_offset, trkerr_dim1, trkerr_dim2, trkerr_offset, u_dim1,
	     u_offset, y_dim1, y_offset, f_dim1, f_offset, g_dim1, g_offset, 
	    h_dim1, h_offset, rr_dim1, rr_offset, a_dim1, a_offset, d_dim1, 
	    d_offset, l_dim1, l_offset, pt1_dim1, pt1_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal detom;
    extern /* Subroutine */ int error_();
    static integer is, it;
    extern /* Subroutine */ int invers_();


/*   Calculate the smoothed state for the model: */

/*        z(t) = Fz(t-1) + Gu(t) + Qe(t) */
/*        y(t) = Hz(t)  + Rw(t) */

/*   see KF and KF.s for details. */
/*   Z should be supplied as the filtered estimate of the state and is */
/*     returned as the smoothed estimate, and similarily for the */
/*     tracking error TRKERR. */

/*   The following are for scratch space */

/*      CALL DBPR('M     ',6, M,1) */
/*      CALL DBPR('N     ',6, N,1) */
/*      CALL DBPR('P     ',6, P,1) */
/*      CALL DBPR('NSMPL  ',6, NSMPL,1) */
    /* Parameter adjustments */
    --zt;
    pt1_dim1 = *n;
    pt1_offset = 1 + pt1_dim1 * 1;
    pt1 -= pt1_offset;
    l_dim1 = *n;
    l_offset = 1 + l_dim1 * 1;
    l -= l_offset;
    d_dim1 = *n;
    d_offset = 1 + d_dim1 * 1;
    d__ -= d_offset;
    a_dim1 = *n;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    f_dim1 = *n;
    f_offset = 1 + f_dim1 * 1;
    f -= f_offset;
    g_dim1 = *n;
    g_offset = 1 + g_dim1 * 1;
    g -= g_offset;
    rr_dim1 = *p;
    rr_offset = 1 + rr_dim1 * 1;
    rr -= rr_offset;
    h_dim1 = *p;
    h_offset = 1 + h_dim1 * 1;
    h__ -= h_offset;
    y_dim1 = *nsmpl;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    u_dim1 = *nsmpl;
    u_offset = 1 + u_dim1 * 1;
    u -= u_offset;
    trkerr_dim1 = *nsmpl;
    trkerr_dim2 = *n;
    trkerr_offset = 1 + trkerr_dim1 * (1 + trkerr_dim2 * 1);
    trkerr -= trkerr_offset;
    z_dim1 = *nsmpl;
    z_offset = 1 + z_dim1 * 1;
    z__ -= z_offset;

    /* Function Body */
    if (is < *n) {
	error_("ERROR: state dimension cannot exceed ", &c__48, &is, &c__1, (
		ftnlen)37);
	return 0;
    }
    if (is < *p) {
	error_("ERROR: output dimensions cannot exceed ", &c__48, &is, &c__1, 
		(ftnlen)39);
	return 0;
    }


/*        RR= RR' */

/*  Next period state and tracking error get clobbered in */
/*    backwards recursion, so: */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* L902: */
	    pt1[i__ + j * pt1_dim1] = trkerr[*nsmpl + (i__ + j * trkerr_dim2) 
		    * trkerr_dim1];
	}
    }
    for (it = *nsmpl - 1; it >= 1; --it) {
/*                                          D=H*P(t|t-1)*H' + RR */
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = *p;
	    for (j = 1; j <= i__1; ++j) {
		a[i__ + j * a_dim1] = 0.;
		i__3 = *n;
		for (k = 1; k <= i__3; ++k) {
/* L2: */
		    a[i__ + j * a_dim1] += trkerr[it + (i__ + k * trkerr_dim2)
			     * trkerr_dim1] * h__[j + k * h_dim1];
		}
	    }
	}
	i__3 = *p;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__1 = *p;
	    for (j = 1; j <= i__1; ++j) {
		d__[i__ + j * d_dim1] = rr[i__ + j * rr_dim1];
		i__2 = *n;
		for (k = 1; k <= i__2; ++k) {
/* L3: */
		    d__[i__ + j * d_dim1] += h__[i__ + k * h_dim1] * a[k + j *
			     a_dim1];
		}
	    }
	}
/*     INVERS inverts in place and RETURNS THE DETERMINANT. */
	invers_(&d__[d_offset], p, &is, &detom);
/*                        Kalman gain  (A=)  K=P(t|t-1)*H'*inv(D) */
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = *p;
	    for (j = 1; j <= i__1; ++j) {
		l[i__ + j * l_dim1] = 0.;
		i__3 = *p;
		for (k = 1; k <= i__3; ++k) {
/* L4: */
		    l[i__ + j * l_dim1] += a[i__ + k * a_dim1] * d__[k + j * 
			    d_dim1];
		}
	    }
	}
/*                       L now contains the Kalman gain K */
/*      IF (IT.EQ.NSMPL-1) THEN */
/*         CALL DBPRDB('K in L ',7,L,IS*IS) */
/*      ENDIF */
/*  E(z(t)|y(t),u(t+1) ZT = z(t|t) = Z(t|t-1) + K*(y-H*Z(t|t-1) - G*u(t) */
	i__3 = *p;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    a[i__ + a_dim1] = y[it + i__ * y_dim1];
	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {
/* L107: */
		a[i__ + a_dim1] -= h__[i__ + k * h_dim1] * z__[it + k * 
			z_dim1];
	    }
	}
	if (*m != 0) {
	    i__1 = *p;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__3 = *m;
		for (k = 1; k <= i__3; ++k) {
/* L108: */
		    a[i__ + a_dim1] -= g[i__ + k * g_dim1] * u[it + k * 
			    u_dim1];
		}
	    }
	}
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    zt[i__] = z__[it + i__ * z_dim1];
	    i__1 = *p;
	    for (k = 1; k <= i__1; ++k) {
/* L109: */
		zt[i__] += l[i__ + k * l_dim1] * a[k + a_dim1];
	    }
	}
/*                                  P(t|t) = P(t|t-1) - K*H*P(t|t-1) */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		a[i__ + j * a_dim1] = 0.;
		i__2 = *p;
		for (k = 1; k <= i__2; ++k) {
/* L7: */
		    a[i__ + j * a_dim1] += l[i__ + k * l_dim1] * h__[k + j * 
			    h_dim1];
		}
	    }
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		l[i__ + j * l_dim1] = trkerr[it + (i__ + j * trkerr_dim2) * 
			trkerr_dim1];
		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
/* L8: */
		    l[i__ + j * l_dim1] -= a[i__ + k * a_dim1] * trkerr[it + (
			    k + j * trkerr_dim2) * trkerr_dim1];
		}
	    }
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
/* L9: */
		a[i__ + j * a_dim1] = (l[i__ + j * l_dim1] + l[j + i__ * 
			l_dim1]) / 2.;
	    }
	}
/*                                        A now contains P(t|t) */
/*      IF (IT.EQ.NSMPL-1) THEN */
/*         CALL DBPRDB('P(t|t) in A ',12,A,IS*IS) */
/*      ENDIF */
/*                                   J = P(t|t)*F'*inv(P(t+1|t)) */
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/* L51: */
		l[i__ + j * l_dim1] = pt1[i__ + j * pt1_dim1];
	    }
	}
	invers_(&l[l_offset], n, &is, &detom);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		d__[i__ + j * d_dim1] = 0.;
		i__2 = *n;
		for (k = 1; k <= i__2; ++k) {
/* L52: */
		    d__[i__ + j * d_dim1] += f[k + i__ * f_dim1] * l[k + j * 
			    l_dim1];
		}
	    }
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		l[i__ + j * l_dim1] = 0.;
		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
/* L53: */
		    l[i__ + j * l_dim1] += a[i__ + k * a_dim1] * d__[k + j * 
			    d_dim1];
		}
	    }
	}
/*                         L now contains J and A now contains P(t|t). */
/*            smoothed state sm[t] = ZT + J*(sm[t+1] - F*ZT - G*u(t+1)) */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d__[i__ + d_dim1] = z__[it + 1 + i__ * z_dim1];
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
/* L16: */
		d__[i__ + d_dim1] -= f[i__ + k * f_dim1] * zt[k];
	    }
	}
	if (*m != 0) {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__1 = *m;
		for (k = 1; k <= i__1; ++k) {
/* L17: */
		    d__[i__ + d_dim1] -= g[i__ + k * g_dim1] * u[it + 1 + k * 
			    u_dim1];
		}
	    }
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    z__[it + i__ * z_dim1] = zt[i__];
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
/* L18: */
		z__[it + i__ * z_dim1] += l[i__ + k * l_dim1] * d__[k + 
			d_dim1];
	    }
	}
/*  smoothed tracking error strk[t]= P[t|t] + J*(strk[t+1]-trk[t+1])*J' */
/*     L contains J and A contains P(t|t) and PT1 contains trk[t+1]. */
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/* L26: */
		pt1[i__ + j * pt1_dim1] = trkerr[it + 1 + (i__ + j * 
			trkerr_dim2) * trkerr_dim1] - pt1[i__ + j * pt1_dim1];
	    }
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		d__[i__ + j * d_dim1] = 0.;
		i__2 = *n;
		for (k = 1; k <= i__2; ++k) {
/* L27: */
		    d__[i__ + j * d_dim1] += pt1[i__ + k * pt1_dim1] * l[j + 
			    k * l_dim1];
		}
	    }
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
/* L29: */
		pt1[i__ + j * pt1_dim1] = trkerr[it + (i__ + j * trkerr_dim2) 
			* trkerr_dim1];
	    }
	}
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		trkerr[it + (i__ + j * trkerr_dim2) * trkerr_dim1] = a[i__ + 
			j * a_dim1];
		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
/* L28: */
		    trkerr[it + (i__ + j * trkerr_dim2) * trkerr_dim1] += l[
			    i__ + k * l_dim1] * d__[k + j * d_dim1];
		}
	    }
	}
/* L1000: */
    }
    return 0;
} /* smooth_ */


/* Subroutine */ int kfp_(ey, hperr, prderr, errwt, m, n, p, nsmpl, npred, 
	nacc, u, y, f, g, h__, fk, q, r__, gain, z0, p0, ith, parm, ap, ip, 
	jp, ict, const__, an, in, jn)
doublereal *ey;
integer *hperr;
doublereal *prderr, *errwt;
integer *m, *n, *p, *nsmpl, *npred, *nacc;
doublereal *u, *y, *f, *g, *h__, *fk, *q, *r__;
integer *gain;
doublereal *z0, *p0;
integer *ith;
doublereal *parm;
integer *ap, *ip, *jp, *ict;
doublereal *const__;
integer *an, *in, *jn;
{
    /* System generated locals */
    integer ey_dim1, ey_offset, prderr_dim1, prderr_offset, y_dim1, y_offset, 
	    u_dim1, u_offset, f_dim1, f_offset, g_dim1, g_offset, h_dim1, 
	    h_offset, fk_dim1, fk_offset, r_dim1, r_offset, q_dim1, q_offset, 
	    p0_dim1, p0_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal state[1]	/* was [1][1] */;
    extern /* Subroutine */ int kf_();
    static integer lstate, ltkerr;
    static doublereal trkerr[1]	/* was [1][1][1] */;


/* Put parameters into arrays (as in S function setArrays) and call KF */

/*  The state and tracking error are not calculated. */
/*       Use KF if these are needed. */

/*  It is assummed that M,N,P, the dimensions of the parameter */
/*    arrays - are given. Trying to calculate these causes problems. */

/*      DOUBLE PRECISION STATE(NPRED,N),TRKERR(NPRED,N,N) */

/* ..bug in S: passing characters is unreliable */
/*   use integer for AP and AN... */
/*  state and trkerr are not used but must be passed to KF */


    /* Parameter adjustments */
    --errwt;
    p0_dim1 = *n;
    p0_offset = 1 + p0_dim1 * 1;
    p0 -= p0_offset;
    --z0;
    q_dim1 = *n;
    q_offset = 1 + q_dim1 * 1;
    q -= q_offset;
    g_dim1 = *n;
    g_offset = 1 + g_dim1 * 1;
    g -= g_offset;
    f_dim1 = *n;
    f_offset = 1 + f_dim1 * 1;
    f -= f_offset;
    r_dim1 = *p;
    r_offset = 1 + r_dim1 * 1;
    r__ -= r_offset;
    fk_dim1 = *n;
    fk_offset = 1 + fk_dim1 * 1;
    fk -= fk_offset;
    h_dim1 = *p;
    h_offset = 1 + h_dim1 * 1;
    h__ -= h_offset;
    prderr_dim1 = *nsmpl;
    prderr_offset = 1 + prderr_dim1 * 1;
    prderr -= prderr_offset;
    ey_dim1 = *npred;
    ey_offset = 1 + ey_dim1 * 1;
    ey -= ey_offset;
    y_dim1 = *nacc;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    u_dim1 = *nacc;
    u_offset = 1 + u_dim1 * 1;
    u -= u_offset;
    --jp;
    --ip;
    --ap;
    --parm;
    --jn;
    --in;
    --an;
    --const__;

    /* Function Body */
    lstate = 0;
    ltkerr = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* L1: */
	    f[i__ + j * f_dim1] = 0.;
	}
    }
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
/* L2: */
	    g[i__ + j * g_dim1] = 0.;
	}
    }
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* L3: */
	    h__[i__ + j * h_dim1] = 0.;
	}
    }
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = *p;
	for (j = 1; j <= i__1; ++j) {
/* L4: */
	    fk[i__ + j * fk_dim1] = 0.;
	}
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* L5: */
	    q[i__ + j * q_dim1] = 0.;
	}
    }
    i__2 = *p;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = *p;
	for (j = 1; j <= i__1; ++j) {
/* L6: */
	    r__[i__ + j * r_dim1] = 0.;
	}
    }
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L7: */
	z0[i__] = 0.;
    }
    if (*ith > 0) {
	i__1 = *ith;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (ap[i__] == 1) {
		f[ip[i__] + jp[i__] * f_dim1] = parm[i__];
	    } else if (ap[i__] == 2) {
		g[ip[i__] + jp[i__] * g_dim1] = parm[i__];
	    } else if (ap[i__] == 3) {
		h__[ip[i__] + jp[i__] * h_dim1] = parm[i__];
	    } else if (ap[i__] == 4) {
		fk[ip[i__] + jp[i__] * fk_dim1] = parm[i__];
	    } else if (ap[i__] == 5) {
		q[ip[i__] + jp[i__] * q_dim1] = parm[i__];
	    } else if (ap[i__] == 6) {
		r__[ip[i__] + jp[i__] * r_dim1] = parm[i__];
	    } else if (ap[i__] == 7) {
		z0[ip[i__]] = parm[i__];
	    } else if (ap[i__] == 8) {
		p0[ip[i__] + jp[i__] * p0_dim1] = parm[i__];
	    }
/* L101: */
	}
    }
    if (*ict > 0) {
	i__1 = *ict;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (an[i__] == 1) {
		f[in[i__] + jn[i__] * f_dim1] = const__[i__];
	    } else if (an[i__] == 2) {
		g[in[i__] + jn[i__] * g_dim1] = const__[i__];
	    } else if (an[i__] == 3) {
		h__[in[i__] + jn[i__] * h_dim1] = const__[i__];
	    } else if (an[i__] == 4) {
		fk[in[i__] + jn[i__] * fk_dim1] = const__[i__];
	    } else if (an[i__] == 5) {
		q[in[i__] + jn[i__] * q_dim1] = const__[i__];
	    } else if (an[i__] == 6) {
		r__[in[i__] + jn[i__] * r_dim1] = const__[i__];
	    } else if (an[i__] == 7) {
		z0[in[i__]] = const__[i__];
	    } else if (an[i__] == 8) {
		p0[in[i__] + jn[i__] * p0_dim1] = const__[i__];
	    }
/* L102: */
	}
    }
    kf_(&ey[ey_offset], hperr, &prderr[prderr_offset], &errwt[1], &lstate, 
	    state, &ltkerr, trkerr, m, n, p, nsmpl, npred, nacc, &u[u_offset],
	     &y[y_offset], &f[f_offset], &g[g_offset], &h__[h_offset], &fk[
	    fk_offset], &q[q_offset], &r__[r_offset], gain, &z0[1], &p0[
	    p0_offset]);
    return 0;
} /* kfp_ */


/* Subroutine */ int kfprj_(proj, dscard, horiz, nho, ey, m, n, p, nacc, u, y,
	 f, g, h__, fk, q, r__, gain, z0, p0)
doublereal *proj;
integer *dscard, *horiz, *nho;
doublereal *ey;
integer *m, *n, *p, *nacc;
doublereal *u, *y, *f, *g, *h__, *fk, *q, *r__;
integer *gain;
doublereal *z0, *p0;
{
    /* System generated locals */
    integer proj_dim1, proj_dim2, proj_offset, ey_dim1, ey_offset, y_dim1, 
	    y_offset, u_dim1, u_offset, f_dim1, f_offset, g_dim1, g_offset, 
	    h_dim1, h_offset, fk_dim1, fk_offset, r_dim1, r_offset, q_dim1, 
	    q_offset, p0_dim1, p0_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, hperr;
    static doublereal state[1]	/* was [1][1] */, errwt[1];
    extern /* Subroutine */ int kf_();
    static integer ho, it, lstate;
    static doublereal prderr[1]	/* was [1][1] */;
    static integer ltkerr, mhoriz;
    static doublereal trkerr[1]	/* was [1][1][1] */;

/*  multiple calls to KF for prediction at given horizons. */
/*     See S program project. */

/*  The state and tracking error are not calculate. */
/*  state and trkerr are not used but must be passed to KF */

    /* Parameter adjustments */
    --horiz;
    p0_dim1 = *n;
    p0_offset = 1 + p0_dim1 * 1;
    p0 -= p0_offset;
    --z0;
    q_dim1 = *n;
    q_offset = 1 + q_dim1 * 1;
    q -= q_offset;
    g_dim1 = *n;
    g_offset = 1 + g_dim1 * 1;
    g -= g_offset;
    f_dim1 = *n;
    f_offset = 1 + f_dim1 * 1;
    f -= f_offset;
    r_dim1 = *p;
    r_offset = 1 + r_dim1 * 1;
    r__ -= r_offset;
    fk_dim1 = *n;
    fk_offset = 1 + fk_dim1 * 1;
    fk -= fk_offset;
    h_dim1 = *p;
    h_offset = 1 + h_dim1 * 1;
    h__ -= h_offset;
    y_dim1 = *nacc;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    u_dim1 = *nacc;
    u_offset = 1 + u_dim1 * 1;
    u -= u_offset;
    ey_dim1 = *nacc;
    ey_offset = 1 + ey_dim1 * 1;
    ey -= ey_offset;
    proj_dim1 = *nho;
    proj_dim2 = *nacc;
    proj_offset = 1 + proj_dim1 * (1 + proj_dim2 * 1);
    proj -= proj_offset;

    /* Function Body */
    lstate = 0;
    ltkerr = 0;
    hperr = 0;
    mhoriz = horiz[1];
    i__1 = *nho;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L1: */
/* Computing MIN */
	i__2 = mhoriz, i__3 = horiz[i__];
	mhoriz = min(i__2,i__3);
    }
    i__2 = *nacc - mhoriz;
    for (it = *dscard; it <= i__2; ++it) {
/*       this assumes HORIZ is sorted in ascending order */
	if (it > *nacc - horiz[*nho]) {
	    --(*nho);
/*            CALL DBPR('NHO   ',6, NHO,1) */
	}
	kf_(&ey[ey_offset], &hperr, prderr, errwt, &lstate, state, &ltkerr, 
		trkerr, m, n, p, &it, nacc, nacc, &u[u_offset], &y[y_offset], 
		&f[f_offset], &g[g_offset], &h__[h_offset], &fk[fk_offset], &
		q[q_offset], &r__[r_offset], gain, &z0[1], &p0[p0_offset]);
	i__3 = *nho;
	for (ho = 1; ho <= i__3; ++ho) {
	    i__1 = *p;
	    for (j = 1; j <= i__1; ++j) {
		proj[ho + (it + horiz[ho] + j * proj_dim2) * proj_dim1] = ey[
			it + horiz[ho] + j * ey_dim1];
/* L4: */
	    }
	}
/* L10: */
    }
    return 0;
} /* kfprj_ */

/* Subroutine */ int kfepr_(cov, dscard, horiz, nh, nt, ey, m, n, p, npred, 
	nacc, u, y, f, g, h__, fk, q, r__, gain, z0, p0)
doublereal *cov;
integer *dscard, *horiz, *nh, *nt;
doublereal *ey;
integer *m, *n, *p, *npred, *nacc;
doublereal *u, *y, *f, *g, *h__, *fk, *q, *r__;
integer *gain;
doublereal *z0, *p0;
{
    /* System generated locals */
    integer cov_dim1, cov_dim2, cov_offset, ey_dim1, ey_offset, y_dim1, 
	    y_offset, u_dim1, u_offset, f_dim1, f_offset, g_dim1, g_offset, 
	    h_dim1, h_offset, fk_dim1, fk_offset, r_dim1, r_offset, q_dim1, 
	    q_offset, p0_dim1, p0_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, hperr;
    static doublereal state[1]	/* was [1][1] */, errwt[1];
    static integer hi;
    static doublereal mf;
    extern /* Subroutine */ int kf_();
    static integer it, lstate;
    static doublereal prderr[1]	/* was [1][1] */;
    static integer ltkerr;
    static doublereal trkerr[1]	/* was [1][1][1] */;

/*  multiple calls to KF for prediction evaluation. */
/*     See S program predictions.cov.TSmodel */

/*  The state and tracking error are not calculate. */
/*  state and trkerr are not used but must be passed to KF */


/*        CALL DBPR('NPRED ',6, NPRED,1) */
/*        CALL DBPR('HORIZ ',6, HORIZ(1),1) */
/*        CALL DBPR('DSCARD',7, DSCARD,1) */
    /* Parameter adjustments */
    --nt;
    --horiz;
    p0_dim1 = *n;
    p0_offset = 1 + p0_dim1 * 1;
    p0 -= p0_offset;
    --z0;
    q_dim1 = *n;
    q_offset = 1 + q_dim1 * 1;
    q -= q_offset;
    g_dim1 = *n;
    g_offset = 1 + g_dim1 * 1;
    g -= g_offset;
    f_dim1 = *n;
    f_offset = 1 + f_dim1 * 1;
    f -= f_offset;
    r_dim1 = *p;
    r_offset = 1 + r_dim1 * 1;
    r__ -= r_offset;
    fk_dim1 = *n;
    fk_offset = 1 + fk_dim1 * 1;
    fk -= fk_offset;
    h_dim1 = *p;
    h_offset = 1 + h_dim1 * 1;
    h__ -= h_offset;
    cov_dim1 = *nh;
    cov_dim2 = *p;
    cov_offset = 1 + cov_dim1 * (1 + cov_dim2 * 1);
    cov -= cov_offset;
    ey_dim1 = *npred;
    ey_offset = 1 + ey_dim1 * 1;
    ey -= ey_offset;
    y_dim1 = *nacc;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    u_dim1 = *nacc;
    u_offset = 1 + u_dim1 * 1;
    u -= u_offset;

    /* Function Body */
    lstate = 0;
    ltkerr = 0;
    hperr = 0;
    i__1 = *nh;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L1: */
	nt[i__] = 0;
    }
    i__1 = *nh;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *p;
	    for (j = 1; j <= i__3; ++j) {
/* L2: */
		cov[k + (i__ + j * cov_dim2) * cov_dim1] = 0.;
	    }
	}
    }
    i__3 = *npred + 1 - horiz[1];
    for (it = *dscard; it <= i__3; ++it) {
	kf_(&ey[ey_offset], &hperr, prderr, errwt, &lstate, state, &ltkerr, 
		trkerr, m, n, p, &it, npred, nacc, &u[u_offset], &y[y_offset],
		 &f[f_offset], &g[g_offset], &h__[h_offset], &fk[fk_offset], &
		q[q_offset], &r__[r_offset], gain, &z0[1], &p0[p0_offset]);
/*       this assumes HORIZ is sorted in ascending order */
	if (it - 1 + horiz[*nh] > *npred) {
	    --(*nh);
	}
	i__2 = *nh;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    hi = it - 1 + horiz[i__];
	    i__1 = *p;
	    for (j = 1; j <= i__1; ++j) {
/* L4: */
		ey[i__ + j * ey_dim1] = ey[hi + j * ey_dim1] - y[hi + j * 
			y_dim1];
	    }
	}
	i__1 = *nh;
	for (k = 1; k <= i__1; ++k) {
	    ++nt[k];
	    mf = (doublereal) (nt[k] - 1) / (doublereal) nt[k];
	    i__2 = *p;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__4 = *p;
		for (j = 1; j <= i__4; ++j) {
		    cov[k + (i__ + j * cov_dim2) * cov_dim1] = cov[k + (i__ + 
			    j * cov_dim2) * cov_dim1] * mf + ey[k + i__ * 
			    ey_dim1] * ey[k + j * ey_dim1] / nt[k];
/* L5: */
		}
	    }
	}
/* L10: */
    }
    return 0;
} /* kfepr_ */

/* Subroutine */ int kf_(ey, hperr, prderr, errwt, lstate, state, ltkerr, 
	trkerr, m, n, p, nsmpl, npred, nacc, u, y, f, g, h__, fk, q, r__, 
	gain, z0, p0)
doublereal *ey;
integer *hperr;
doublereal *prderr, *errwt;
integer *lstate;
doublereal *state;
integer *ltkerr;
doublereal *trkerr;
integer *m, *n, *p, *nsmpl, *npred, *nacc;
doublereal *u, *y, *f, *g, *h__, *fk, *q, *r__;
integer *gain;
doublereal *z0, *p0;
{
    /* System generated locals */
    integer ey_dim1, ey_offset, prderr_dim1, prderr_offset, state_dim1, 
	    state_offset, trkerr_dim1, trkerr_dim2, trkerr_offset, p0_dim1, 
	    p0_offset, y_dim1, y_offset, u_dim1, u_offset, f_dim1, f_offset, 
	    g_dim1, g_offset, h_dim1, h_offset, r_dim1, r_offset, q_dim1, 
	    q_offset, fk_dim1, fk_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static doublereal a[10000]	/* was [100][100] */;
    static integer i__, j, k;
    static doublereal z__[100], detom;
    static integer lperr;
    extern /* Subroutine */ int error_();
    static doublereal aa[10000]	/* was [100][100] */;
    static integer it;
    static doublereal hw, pp[10000]	/* was [100][100] */, qq[10000]	/* 
	    was [100][100] */, rr[10000]	/* was [100][100] */, ww[100],
	     zz[100];
    extern /* Subroutine */ int invers_();


/*   Calculate the likelihood value for the model: */

/*        z(t) = Fz(t-1) + Gu(t) + Qe(t-1) */
/*        y(t) = Hz(t)  + Rw(t) */

/* or the innovations model: */
/*        z(t) = Fz(t-1) + Gu(t) + FKw(t-1) */
/*        y(t) = Hz(t)  + w(t) */


/*  FK is the Kalman gain */
/*  If GAIN is true then FK is taken as given (innovations model) */

/*  M is the dimension of the input u. */
/*  N is the dimension of the state z and the system noise e. */
/*  P is the dimension of the output y and the ouput noise w. */
/*  NSMPL is the length of the data series to use for residual */
/*      and likelihood calculations. */
/*  NPRED is the period to predict ahead (past NSMPL) (for z and y) */
/*  NACC is the actual first (time) dimension of Y and U. */

/*   STATE is the one step ahead ESTIMATE OF STATE. */
/*   It is returned only if LSTATE is TRUE. */
/*   Z0  is the initial state (often set to zero). */
/*   P0  is the initial state tracking error (often set to I and */
/*   totally ignored in innovations models). */
/*   PP is  the one step ahead est. cov matrix of the state estimation error. */
/*   TRKERR is the history of PP at each period. */
/*   It is returned only for non-innovations models (GAIN=FALSE) and */
/*     then only if LTKERR is TRUE. */
/*   EY is the output prediction. EY is used to store WW during computation! */
/*   The prediction error at each period is WW (innovations) = Y - EY. */
/*     If HPERR is equal or greater than one then weighted prediction */
/*     errors are calculated up to the horizon indicated */
/*     by HPERR. The weights taken from ERRWT are applied to the squared */
/*     error at each period ahead. */
/*   If HPERR is zero and LSTATE and LTKERR are false then ERRWT, */
/*    PRDERR, STATE, and TRKERR are not referenced, */
/*    so KF can be called with dummy arguments as */
/*    in KFP and GEND. */

/*   IS is the maximum state dimension  and the maximum output */
/*     dimension  (used for working arrays) */
/*    NSTART is not properly implemented and must be set to 1. */




/*      CALL DBPR('M     ',6, M,1) */
/*      CALL DBPR('N     ',6, N,1) */
/*      CALL DBPR('P     ',6, P,1) */
/*      CALL DBPR('NSMPL  ',6, NSMPL,1) */
/*      CALL DBPR('NPRED  ',6, NPRED,1) */
/*      CALL DBPR('NACC  ',6, NACC,1) */
/*      CALL DBPRDB('R  ',3,R,9) */
/*      IF (GAIN.EQ.1) THEN */
/*         CALL DBPR('GAIN is T',9, 1,0) */
/*      ELSE */
/*         CALL DBPR('GAIN is F',9, 1,0) */
/*      ENDIF */
/*      CALL DBPRDB('Q  ',3,Q,N*N) */
/*      CALL DBPR('HPERR ',6, HPERR,1) */
    /* Parameter adjustments */
    --errwt;
    p0_dim1 = *n;
    p0_offset = 1 + p0_dim1 * 1;
    p0 -= p0_offset;
    --z0;
    q_dim1 = *n;
    q_offset = 1 + q_dim1 * 1;
    q -= q_offset;
    g_dim1 = *n;
    g_offset = 1 + g_dim1 * 1;
    g -= g_offset;
    f_dim1 = *n;
    f_offset = 1 + f_dim1 * 1;
    f -= f_offset;
    r_dim1 = *p;
    r_offset = 1 + r_dim1 * 1;
    r__ -= r_offset;
    fk_dim1 = *n;
    fk_offset = 1 + fk_dim1 * 1;
    fk -= fk_offset;
    h_dim1 = *p;
    h_offset = 1 + h_dim1 * 1;
    h__ -= h_offset;
    prderr_dim1 = *nsmpl;
    prderr_offset = 1 + prderr_dim1 * 1;
    prderr -= prderr_offset;
    trkerr_dim1 = *npred;
    trkerr_dim2 = *n;
    trkerr_offset = 1 + trkerr_dim1 * (1 + trkerr_dim2 * 1);
    trkerr -= trkerr_offset;
    state_dim1 = *npred;
    state_offset = 1 + state_dim1 * 1;
    state -= state_offset;
    ey_dim1 = *npred;
    ey_offset = 1 + ey_dim1 * 1;
    ey -= ey_offset;
    y_dim1 = *nacc;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    u_dim1 = *nacc;
    u_offset = 1 + u_dim1 * 1;
    u -= u_offset;

    /* Function Body */
    if (100 < *n) {
	error_("ERROR: state dimension cannot exceed ", &c__48, &c__100, &
		c__1, (ftnlen)37);
	return 0;
    }
    if (100 < *p) {
	error_("ERROR: output dimensions cannot exceed ", &c__48, &c__100, &
		c__1, (ftnlen)39);
	return 0;
    }
    lperr = 0;
    if (*hperr > 0) {
	i__1 = *hperr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (errwt[i__] > (float)0.) {
		lperr = 1;
	    }
/* L500: */
	}
    }

/*     initial innovation. */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L209: */
	ww[i__ - 1] = (float)0.;
    }
/*     initial state. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L210: */
	zz[i__ - 1] = z0[i__];
    }

    if (*gain != 1) {
/*        initial tracking error. */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
/* L220: */
		pp[i__ + j * 100 - 101] = p0[i__ + j * p0_dim1];
	    }
	}

/*        RR= RR'    QQ= QQ' */

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		qq[i__ + j * 100 - 101] = 0.;
		i__3 = *n;
		for (k = 1; k <= i__3; ++k) {
/* L230: */
		    qq[i__ + j * 100 - 101] += q[i__ + k * q_dim1] * q[j + k *
			     q_dim1];
		}
	    }
	}
	i__3 = *p;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__1 = *p;
	    for (j = 1; j <= i__1; ++j) {
		rr[i__ + j * 100 - 101] = 0.;
		i__2 = *p;
		for (k = 1; k <= i__2; ++k) {
/* L236: */
		    rr[i__ + j * 100 - 101] += r__[i__ + k * r_dim1] * r__[j 
			    + k * r_dim1];
		}
	    }
	}
    }

/*    Start Time loop */

    i__2 = *nsmpl;
    for (it = 1; it <= i__2; ++it) {
	if (*gain != 1) {


/*   Kalman gain  FK = F*P(t|t-1)*H'* inv( H*P(t|t-1)*H' + RR') */

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__3 = *p;
		for (j = 1; j <= i__3; ++j) {
		    aa[i__ + j * 100 - 101] = 0.;
		    i__4 = *n;
		    for (k = 1; k <= i__4; ++k) {
/* L5: */
			aa[i__ + j * 100 - 101] += pp[i__ + k * 100 - 101] * 
				h__[j + k * h_dim1];
		    }
		}
	    }
	    i__4 = *n;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__3 = *p;
		for (j = 1; j <= i__3; ++j) {
		    fk[i__ + j * fk_dim1] = 0.;
		    i__1 = *n;
		    for (k = 1; k <= i__1; ++k) {
/* L6: */
			fk[i__ + j * fk_dim1] += f[i__ + k * f_dim1] * aa[k + 
				j * 100 - 101];
		    }
		}
	    }
	    i__1 = *p;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__3 = *p;
		for (j = 1; j <= i__3; ++j) {
		    a[i__ + j * 100 - 101] = rr[i__ + j * 100 - 101];
		    i__4 = *n;
		    for (k = 1; k <= i__4; ++k) {
/* L7: */
			a[i__ + j * 100 - 101] += h__[i__ + k * h_dim1] * aa[
				k + j * 100 - 101];
		    }
		}
	    }
/*      CALL DBPRDB('DO 7 A ',7, A,(IS*IS)) */
/*     force symetry. */
	    i__4 = *p;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__3 = *p;
		for (j = 1; j <= i__3; ++j) {
/* L9: */
		    aa[i__ + j * 100 - 101] = (a[i__ + j * 100 - 101] + a[j + 
			    i__ * 100 - 101]) / 2.;
		}
	    }

/*     INVERS inverts in place and RETURNS THE DETERMINANT. */
	    invers_(aa, p, &c__100, &detom);
/*     DETOM SHOULD BE POSITIVE */

	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__4 = *p;
		for (j = 1; j <= i__4; ++j) {
		    a[i__ + j * 100 - 101] = 0.;
		    i__1 = *p;
		    for (k = 1; k <= i__1; ++k) {
/* L14: */
			a[i__ + j * 100 - 101] += fk[i__ + k * fk_dim1] * aa[
				k + j * 100 - 101];
		    }
		}
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__4 = *p;
		for (j = 1; j <= i__4; ++j) {
/* L15: */
		    fk[i__ + j * fk_dim1] = a[i__ + j * 100 - 101];
		}
	    }
/*      IF (IT.EQ.1) THEN */
/*         CALL DBPRDB('K ',2,FK,N*P) */
/*      ENDIF */

/*     P(t|t-1)= F*P(t-|t-2)*F' -  K*H*P(t-1|t-2)*F' + Q*Q' */

	    i__4 = *n;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    a[i__ + j * 100 - 101] = 0.;
		    i__3 = *n;
		    for (k = 1; k <= i__3; ++k) {
/* L2: */
			a[i__ + j * 100 - 101] += pp[i__ + k * 100 - 101] * f[
				j + k * f_dim1];
		    }
		}
	    }
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    aa[i__ + j * 100 - 101] = qq[i__ + j * 100 - 101];
		    i__4 = *n;
		    for (k = 1; k <= i__4; ++k) {
/* L3: */
			aa[i__ + j * 100 - 101] += f[i__ + k * f_dim1] * a[k 
				+ j * 100 - 101];
		    }
		}
	    }
	    i__4 = *p;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    pp[i__ + j * 100 - 101] = 0.;
		    i__3 = *n;
		    for (k = 1; k <= i__3; ++k) {
/* L18: */
			pp[i__ + j * 100 - 101] += h__[i__ + k * h_dim1] * a[
				k + j * 100 - 101];
		    }
		}
	    }
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__4 = *p;
		    for (k = 1; k <= i__4; ++k) {
/* L19: */
			aa[i__ + j * 100 - 101] -= fk[i__ + k * fk_dim1] * pp[
				k + j * 100 - 101];
		    }
		}
	    }
/*     force symetry to avoid numerical round off problems */
	    i__4 = *n;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
/* L20: */
		    pp[i__ + j * 100 - 101] = (aa[i__ + j * 100 - 101] + aa[j 
			    + i__ * 100 - 101]) / 2.;
		}
	    }
/*      IF (IT.EQ.1) THEN */
/*         CALL DBPRDB('PP  ',4,PP,IS*IS) */
/*      ENDIF */
	    if (*ltkerr == 1) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__4 = *n;
		    for (j = 1; j <= i__4; ++j) {
/* L21: */
			trkerr[it + (i__ + j * trkerr_dim2) * trkerr_dim1] = 
				pp[i__ + j * 100 - 101];
		    }
		}
	    }
	}
/*   end of Kalman gain and tracking error update ( if NOT GAIN ) */

/*   one step ahead state estimate */
/*     z(t|t-1)= Fz(t-1|t-2) + FK*WW(t-1) + Gu(t) */

	i__4 = *n;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    z__[i__ - 1] = (float)0.;
	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {
/* L1: */
		z__[i__ - 1] += f[i__ + k * f_dim1] * zz[k - 1];
	    }
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__4 = *p;
	    for (k = 1; k <= i__4; ++k) {
/* L22: */
		z__[i__ - 1] += fk[i__ + k * fk_dim1] * ww[k - 1];
	    }
	}
	i__4 = *n;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    i__1 = *m;
	    for (k = 1; k <= i__1; ++k) {
/* L23: */
		z__[i__ - 1] += g[i__ + k * g_dim1] * u[it + k * u_dim1];
	    }
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L24: */
	    zz[i__ - 1] = z__[i__ - 1];
	}
	if (*lstate == 1) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L25: */
		state[it + i__ * state_dim1] = z__[i__ - 1];
	    }
	}
/*      CALL DBPRDB('Z  ',3,Z,N) */
/*  one step ahead prediction  EY(t)=H*z(t|t-1)) */

/*  innovations  WW(t)=y(t)-H*z(t|t-1)) */

/*   EY stores history of predition error WW to reconstruct predictions */
	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ww[i__ - 1] = y[it + i__ * y_dim1];
	    i__4 = *n;
	    for (j = 1; j <= i__4; ++j) {
/* L10: */
		ww[i__ - 1] -= h__[i__ + j * h_dim1] * z__[j - 1];
	    }
	}
	i__4 = *p;
	for (i__ = 1; i__ <= i__4; ++i__) {
/* L11: */
	    ey[it + i__ * ey_dim1] = ww[i__ - 1];
	}
/*      CALL DBPRDB('WW ',3,WW,P) */
/*   Return weighted prediction error */
	if (lperr == 1) {
	    i__4 = *p;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/* L401: */
/* Computing 2nd power */
		d__1 = ww[i__ - 1];
		prderr[it + i__ * prderr_dim1] = errwt[1] * (d__1 * d__1);
	    }
	    if (*hperr > 1) {
		i__4 = *hperr;
		for (k = 2; k <= i__4; ++k) {
		    if (it + k - 1 <= *nsmpl) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
/* L402: */
			    aa[i__ * 100 - 100] = z__[i__ - 1];
			}
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
/* L409: */
			    z__[i__ - 1] = (float)0.;
			}
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    i__3 = *n;
			    for (j = 1; j <= i__3; ++j) {
/* L403: */
				z__[i__ - 1] += f[i__ + j * f_dim1] * aa[j * 
					100 - 100];
			    }
			}
			if (k == 2) {
			    i__3 = *n;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				i__1 = *p;
				for (j = 1; j <= i__1; ++j) {
/* L406: */
				    z__[i__ - 1] += fk[i__ + j * fk_dim1] * 
					    ww[j - 1];
				}
			    }
			}
			if (*m != 0) {
			    i__1 = *n;
			    for (i__ = 1; i__ <= i__1; ++i__) {
				i__3 = *m;
				for (j = 1; j <= i__3; ++j) {
/* L404: */
				    z__[i__ - 1] += g[i__ + j * g_dim1] * u[
					    it + k - 1 + j * u_dim1];
				}
			    }
			}
			i__3 = *p;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    hw = (float)0.;
			    i__1 = *n;
			    for (j = 1; j <= i__1; ++j) {
/* L407: */
				hw += h__[i__ + j * h_dim1] * z__[j - 1];
			    }
/* L408: */
/* Computing 2nd power */
			    d__1 = y[it + k - 1 + i__ * y_dim1] - hw;
			    prderr[it + i__ * prderr_dim1] += errwt[k] * (
				    d__1 * d__1);
			}
/*                 IF (IT.GT.97) THEN */
/*                    CALL DBPR('K     ',6, K,1) */
/*                    CALL DBPRDB('PRDERR ',8,PRDERR,NSMPL*P) */
/*                 ENDIF */
		    }
/* L400: */
		}
	    }
	}
/* L1000: */
    }
/*  end of time loop */
/*     reconstruct predictions */
    i__2 = *nsmpl;
    for (it = 1; it <= i__2; ++it) {
	i__4 = *p;
	for (i__ = 1; i__ <= i__4; ++i__) {
/* L41: */
	    ey[it + i__ * ey_dim1] = y[it + i__ * y_dim1] - ey[it + i__ * 
		    ey_dim1];
	}
    }

/*    Start multi-step prediction loop */

    if (*npred > *nsmpl) {
	i__4 = *npred;
	for (it = *nsmpl + 1; it <= i__4; ++it) {

	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		z__[i__ - 1] = (float)0.;
		i__3 = *n;
		for (k = 1; k <= i__3; ++k) {
/* L2001: */
		    z__[i__ - 1] += f[i__ + k * f_dim1] * zz[k - 1];
		}
	    }
	    if (*m != 0) {
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__2 = *m;
		    for (k = 1; k <= i__2; ++k) {
/* L2003: */
			z__[i__ - 1] += g[i__ + k * g_dim1] * u[it + k * 
				u_dim1];
		    }
		}
	    }
/*   use prediction error from last sample point (i.e. for first prediction) */
	    if (it == *nsmpl + 1) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__3 = *p;
		    for (k = 1; k <= i__3; ++k) {
/* L2004: */
			z__[i__ - 1] += fk[i__ + k * fk_dim1] * ww[k - 1];
		    }
		}
	    }
	    if (*lstate == 1) {
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
/* L2005: */
		    state[it + i__ * state_dim1] = z__[i__ - 1];
		}
	    }
	    i__3 = *p;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ey[it + i__ * ey_dim1] = (float)0.;
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
/* L2010: */
		    ey[it + i__ * ey_dim1] += h__[i__ + j * h_dim1] * z__[j - 
			    1];
		}
	    }
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L2024: */
		zz[i__ - 1] = z__[i__ - 1];
	    }
/* L2000: */
	}
    }
/*  end of multi-step prediction loop */
    return 0;
} /* kf_ */



/*        1         2         3         4         5         6         7         8 */
/* Subroutine */ int simrma_(y, y0, m, p, ia, ib, ic, nsmpl, u, u0, w, w0, a, 
	b, c__, trend)
doublereal *y, *y0;
integer *m, *p, *ia, *ib, *ic, *nsmpl;
doublereal *u, *u0, *w, *w0, *a, *b, *c__, *trend;
{
    /* System generated locals */
    integer y_dim1, y_offset, u_dim1, u_offset, y0_dim1, y0_offset, u0_dim1, 
	    u0_offset, w_dim1, w_offset, w0_dim1, w0_offset, a_dim1, a_dim2, 
	    a_offset, b_dim1, b_dim2, b_offset, c_dim1, c_dim2, c_offset, 
	    i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, l, it;

/*      Simulate an ARMA model. See documentation in ARMA and in the S version. */
/*       CALL DBPRDB('inSIMARMA ',7, 1,1) */
/*       CALL DBPR('M     ',6, M,1) */
/*       CALL DBPRDB('A     ',6, A,(IA*P*P)) */
/*      IF (IT.LE.5) CALL DBPRDB('step1 ',6, Y(IT,3),1) */
/*      CALL DBPRDB('C ',2, C,(IC*P*M) ) */
/*      CALL DBPRDB('U0 ',3, U0,(IC*M) ) */
/*      CALL DBPRDB('U  ',3, U,(NSMPL*M) ) */
    /* Parameter adjustments */
    --trend;
    a_dim1 = *ia;
    a_dim2 = *p;
    a_offset = 1 + a_dim1 * (1 + a_dim2 * 1);
    a -= a_offset;
    y0_dim1 = *ia;
    y0_offset = 1 + y0_dim1 * 1;
    y0 -= y0_offset;
    b_dim1 = *ib;
    b_dim2 = *p;
    b_offset = 1 + b_dim1 * (1 + b_dim2 * 1);
    b -= b_offset;
    w0_dim1 = *ib;
    w0_offset = 1 + w0_dim1 * 1;
    w0 -= w0_offset;
    c_dim1 = *ic;
    c_dim2 = *p;
    c_offset = 1 + c_dim1 * (1 + c_dim2 * 1);
    c__ -= c_offset;
    u0_dim1 = *ic;
    u0_offset = 1 + u0_dim1 * 1;
    u0 -= u0_offset;
    w_dim1 = *nsmpl;
    w_offset = 1 + w_dim1 * 1;
    w -= w_offset;
    u_dim1 = *nsmpl;
    u_offset = 1 + u_dim1 * 1;
    u -= u_offset;
    y_dim1 = *nsmpl;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *nsmpl;
	for (it = 1; it <= i__2; ++it) {
/* L2001: */
	    y[it + i__ * y_dim1] = 0.;
	}
    }
    i__2 = *nsmpl;
    for (it = 1; it <= i__2; ++it) {
	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L1: */
	    y[it + i__ * y_dim1] = trend[i__];
	}
/*      IF (IT.LE.5) CALL DBPRDB('step1 ',6, Y(IT,3),1) */
	i__1 = *ia;
	for (l = 2; l <= i__1; ++l) {
	    if (it + 1 <= l) {
		i__3 = *p;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__4 = *p;
		    for (j = 1; j <= i__4; ++j) {
/* L2: */
			y[it + i__ * y_dim1] -= a[l + (i__ + j * a_dim2) * 
				a_dim1] * y0[l - it + j * y0_dim1];
		    }
		}
	    } else {
		i__4 = *p;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    i__3 = *p;
		    for (j = 1; j <= i__3; ++j) {
/* L3: */
			y[it + i__ * y_dim1] -= a[l + (i__ + j * a_dim2) * 
				a_dim1] * y[it + 1 - l + j * y_dim1];
		    }
		}
	    }
/* L5: */
	}
/*      IF (IT.LE.5) CALL DBPRDB('step2 ',6, Y(IT,3),1) */
	i__1 = *ib;
	for (l = 1; l <= i__1; ++l) {
	    if (it + 1 <= l) {
		i__3 = *p;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__4 = *p;
		    for (j = 1; j <= i__4; ++j) {
/* L12: */
			y[it + i__ * y_dim1] += b[l + (i__ + j * b_dim2) * 
				b_dim1] * w0[l - it + j * w0_dim1];
		    }
		}
	    } else {
		i__4 = *p;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    i__3 = *p;
		    for (j = 1; j <= i__3; ++j) {
/* L13: */
			y[it + i__ * y_dim1] += b[l + (i__ + j * b_dim2) * 
				b_dim1] * w[it + 1 - l + j * w_dim1];
		    }
		}
	    }
/* L15: */
	}
/*      IF (IT.LE.5) CALL DBPRDB('step3 ',6, Y(IT,3),1) */
	if (*m > 0) {
	    i__1 = *ic;
	    for (l = 1; l <= i__1; ++l) {
		if (it + 1 <= l) {
		    i__3 = *p;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			i__4 = *m;
			for (j = 1; j <= i__4; ++j) {
/* L22: */
			    y[it + i__ * y_dim1] += c__[l + (i__ + j * c_dim2)
				     * c_dim1] * u0[l - it + j * u0_dim1];
			}
		    }
/*         IF (IT.LE.5) CALL DBPRDB('stepa ',6, Y(IT,3),1) */
		} else {
		    i__4 = *p;
		    for (i__ = 1; i__ <= i__4; ++i__) {
			i__3 = *m;
			for (j = 1; j <= i__3; ++j) {
/* L23: */
			    y[it + i__ * y_dim1] += c__[l + (i__ + j * c_dim2)
				     * c_dim1] * u[it + 1 - l + j * u_dim1];
			}
		    }
/*           IF (IT.LE.5) CALL DBPRDB('stepb ',6, Y(IT,3),1) */
		}
/* L25: */
	    }
	}
/*      IF (IT.LE.5) CALL DBPRDB('step4 ',6, Y(IT,3),1) */
/* L1000: */
    }
    return 0;
} /* simrma_ */

/* Subroutine */ int armap_(ey, hperr, prderr, errwt, m, p, ia, ib, ic, nsmpl,
	 npred, nacc, u, y, a, b, c__, trend, ith, parm, ap, lp, ip, jp, ict, 
	const__, an, ln, in, jn, is, aa, bb, ww)
doublereal *ey;
integer *hperr;
doublereal *prderr, *errwt;
integer *m, *p, *ia, *ib, *ic, *nsmpl, *npred, *nacc;
doublereal *u, *y, *a, *b, *c__, *trend;
integer *ith;
doublereal *parm;
integer *ap, *lp, *ip, *jp, *ict;
doublereal *const__;
integer *an, *ln, *in, *jn, *is;
doublereal *aa, *bb, *ww;
{
    /* System generated locals */
    integer ey_dim1, ey_offset, prderr_dim1, prderr_offset, y_dim1, y_offset, 
	    u_dim1, u_offset, a_dim1, a_dim2, a_offset, b_dim1, b_dim2, 
	    b_offset, c_dim1, c_dim2, c_offset, aa_dim1, aa_offset, bb_dim1, 
	    bb_offset, i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int arma_();
    static integer i__, j, l;


/* Put parameters into arrays (as in S function setArrays) and call ARMA */

/*  It is assummed that M,P,IA,IB, and IC - the dimensions of the parameter */
/*    arrays - are given. Trying to calculate these causes problems. */


/* ..bug in S: passing characters is unreliable */
/*   use integer for AP and AN... */


    /* Parameter adjustments */
    --errwt;
    --trend;
    a_dim1 = *ia;
    a_dim2 = *p;
    a_offset = 1 + a_dim1 * (1 + a_dim2 * 1);
    a -= a_offset;
    b_dim1 = *ib;
    b_dim2 = *p;
    b_offset = 1 + b_dim1 * (1 + b_dim2 * 1);
    b -= b_offset;
    c_dim1 = *ic;
    c_dim2 = *p;
    c_offset = 1 + c_dim1 * (1 + c_dim2 * 1);
    c__ -= c_offset;
    prderr_dim1 = *nsmpl;
    prderr_offset = 1 + prderr_dim1 * 1;
    prderr -= prderr_offset;
    ey_dim1 = *npred;
    ey_offset = 1 + ey_dim1 * 1;
    ey -= ey_offset;
    y_dim1 = *nacc;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    u_dim1 = *nacc;
    u_offset = 1 + u_dim1 * 1;
    u -= u_offset;
    --jp;
    --ip;
    --lp;
    --ap;
    --parm;
    --jn;
    --in;
    --ln;
    --an;
    --const__;
    --ww;
    bb_dim1 = *is;
    bb_offset = 1 + bb_dim1 * 1;
    bb -= bb_offset;
    aa_dim1 = *is;
    aa_offset = 1 + aa_dim1 * 1;
    aa -= aa_offset;

    /* Function Body */
    i__1 = *ia;
    for (l = 1; l <= i__1; ++l) {
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *p;
	    for (j = 1; j <= i__3; ++j) {
/* L1: */
		a[l + (i__ + j * a_dim2) * a_dim1] = (float)0.;
	    }
	}
    }
    i__3 = *ib;
    for (l = 1; l <= i__3; ++l) {
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = *p;
	    for (j = 1; j <= i__1; ++j) {
/* L2: */
		b[l + (i__ + j * b_dim2) * b_dim1] = (float)0.;
	    }
	}
    }
    i__1 = *ic;
    for (l = 1; l <= i__1; ++l) {
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *m;
	    for (j = 1; j <= i__3; ++j) {
/* L3: */
		c__[l + (i__ + j * c_dim2) * c_dim1] = (float)0.;
	    }
	}
    }
    i__3 = *p;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L4: */
	trend[i__] = (float)0.;
    }
/*       IA=0 */
/*       IB=0 */
/*       IC=0 */
    if (*ith > 0) {
	i__3 = *ith;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    if (ap[i__] == 1) {
		a[lp[i__] + (ip[i__] + jp[i__] * a_dim2) * a_dim1] = parm[i__]
			;
/*               IA=MAX(IA,LP(I)) */
	    } else if (ap[i__] == 2) {
		b[lp[i__] + (ip[i__] + jp[i__] * b_dim2) * b_dim1] = parm[i__]
			;
/*               IB=MAX(IB,LP(I)) */
	    } else if (ap[i__] == 3) {
		c__[lp[i__] + (ip[i__] + jp[i__] * c_dim2) * c_dim1] = parm[
			i__];
/*               IC=MAX(IC,LP(I)) */
	    } else if (ap[i__] == 4) {
		trend[ip[i__]] = parm[i__];
	    }
/* L101: */
	}
    }
    if (*ict > 0) {
	i__3 = *ict;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    if (an[i__] == 1) {
		a[ln[i__] + (in[i__] + jn[i__] * a_dim2) * a_dim1] = const__[
			i__];
/*               IA=MAX(IA,LN(I)) */
	    } else if (an[i__] == 2) {
		b[ln[i__] + (in[i__] + jn[i__] * b_dim2) * b_dim1] = const__[
			i__];
/*               IB=MAX(IB,LN(I)) */
	    } else if (an[i__] == 3) {
		c__[ln[i__] + (in[i__] + jn[i__] * c_dim2) * c_dim1] = 
			const__[i__];
/*               IC=MAX(IC,LN(I)) */
	    } else if (an[i__] == 4) {
		trend[in[i__]] = const__[i__];
	    }
/* L102: */
	}
    }
    arma_(&ey[ey_offset], hperr, &prderr[prderr_offset], &errwt[1], m, p, ia, 
	    ib, ic, nsmpl, npred, nacc, &u[u_offset], &y[y_offset], &a[
	    a_offset], &b[b_offset], &c__[c_offset], &trend[1], is, &aa[
	    aa_offset], &bb[bb_offset], &ww[1]);
    return 0;
} /* armap_ */

/* Subroutine */ int rmaprj_(proj, dscard, horiz, nho, ey, m, p, ia, ib, ic, 
	nacc, u, y, a, b, c__, trend, is, aa, bb, ww)
doublereal *proj;
integer *dscard, *horiz, *nho;
doublereal *ey;
integer *m, *p, *ia, *ib, *ic, *nacc;
doublereal *u, *y, *a, *b, *c__, *trend;
integer *is;
doublereal *aa, *bb, *ww;
{
    /* System generated locals */
    integer proj_dim1, proj_dim2, proj_offset, ey_dim1, ey_offset, y_dim1, 
	    y_offset, u_dim1, u_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, c_dim1, c_dim2, c_offset, aa_dim1, aa_offset, 
	    bb_dim1, bb_offset, i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int arma_();
    static integer i__, j, hperr;
    static doublereal errwt[1];
    static integer ho, it;
    static doublereal prderr[1]	/* was [1][1] */;
    static integer mhoriz;

/*  multiple calls to ARMA for for prediction at given horizons. */
/*     See S program horizonForecasts.TSmodel */

/*  Note: If DSCARD is too small then forecasting starts based on little (or */
/*          no) data and the results will be spurious. */
    /* Parameter adjustments */
    --horiz;
    --trend;
    a_dim1 = *ia;
    a_dim2 = *p;
    a_offset = 1 + a_dim1 * (1 + a_dim2 * 1);
    a -= a_offset;
    b_dim1 = *ib;
    b_dim2 = *p;
    b_offset = 1 + b_dim1 * (1 + b_dim2 * 1);
    b -= b_offset;
    c_dim1 = *ic;
    c_dim2 = *p;
    c_offset = 1 + c_dim1 * (1 + c_dim2 * 1);
    c__ -= c_offset;
    y_dim1 = *nacc;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    u_dim1 = *nacc;
    u_offset = 1 + u_dim1 * 1;
    u -= u_offset;
    ey_dim1 = *nacc;
    ey_offset = 1 + ey_dim1 * 1;
    ey -= ey_offset;
    proj_dim1 = *nho;
    proj_dim2 = *nacc;
    proj_offset = 1 + proj_dim1 * (1 + proj_dim2 * 1);
    proj -= proj_offset;
    --ww;
    bb_dim1 = *is;
    bb_offset = 1 + bb_dim1 * 1;
    bb -= bb_offset;
    aa_dim1 = *is;
    aa_offset = 1 + aa_dim1 * 1;
    aa -= aa_offset;

    /* Function Body */
    hperr = 0;
    mhoriz = horiz[1];
    i__1 = *nho;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L1: */
/* Computing MIN */
	i__2 = mhoriz, i__3 = horiz[i__];
	mhoriz = min(i__2,i__3);
    }
    i__2 = *nacc - mhoriz;
    for (it = *dscard; it <= i__2; ++it) {
/*       this assumes HORIZ is sorted in ascending order */
	if (it > *nacc - horiz[*nho]) {
	    --(*nho);
	}
	arma_(&ey[ey_offset], &hperr, prderr, errwt, m, p, ia, ib, ic, &it, 
		nacc, nacc, &u[u_offset], &y[y_offset], &a[a_offset], &b[
		b_offset], &c__[c_offset], &trend[1], is, &aa[aa_offset], &bb[
		bb_offset], &ww[1]);
	i__3 = *nho;
	for (ho = 1; ho <= i__3; ++ho) {
	    i__1 = *p;
	    for (j = 1; j <= i__1; ++j) {
/* L4: */
		proj[ho + (it + horiz[ho] + j * proj_dim2) * proj_dim1] = ey[
			it + horiz[ho] + j * ey_dim1];
	    }
	}
/* L10: */
    }
    return 0;
} /* rmaprj_ */

/* Subroutine */ int rmaepr_(cov, dscard, horiz, nh, nt, ey, m, p, ia, ib, ic,
	 npred, nacc, u, y, a, b, c__, trend, is, aa, bb, ww)
doublereal *cov;
integer *dscard, *horiz, *nh, *nt;
doublereal *ey;
integer *m, *p, *ia, *ib, *ic, *npred, *nacc;
doublereal *u, *y, *a, *b, *c__, *trend;
integer *is;
doublereal *aa, *bb, *ww;
{
    /* System generated locals */
    integer cov_dim1, cov_dim2, cov_offset, ey_dim1, ey_offset, y_dim1, 
	    y_offset, u_dim1, u_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, c_dim1, c_dim2, c_offset, aa_dim1, aa_offset, 
	    bb_dim1, bb_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    extern /* Subroutine */ int arma_();
    static integer i__, j, k, hperr;
    static doublereal errwt[1];
    static integer hi;
    static doublereal mf;
    static integer it;
    static doublereal prderr[1]	/* was [1][1] */;

/*  multiple calls to ARMA for prediction analysis. */
/*     See S program forecastCov */

/*  Note: If DSCARD is too small then forecasting starts based on little (or */
/*          no) data and the results will be spurious. */

/*        CALL DBPR('NPRED ',6, NPRED,1) */
/*        CALL DBPR('HORIZ ',6, HORIZ(1),1) */
/*        CALL DBPR('DSCARD',7, DSCARD,1) */
    /* Parameter adjustments */
    --nt;
    --horiz;
    --trend;
    cov_dim1 = *nh;
    cov_dim2 = *p;
    cov_offset = 1 + cov_dim1 * (1 + cov_dim2 * 1);
    cov -= cov_offset;
    a_dim1 = *ia;
    a_dim2 = *p;
    a_offset = 1 + a_dim1 * (1 + a_dim2 * 1);
    a -= a_offset;
    b_dim1 = *ib;
    b_dim2 = *p;
    b_offset = 1 + b_dim1 * (1 + b_dim2 * 1);
    b -= b_offset;
    c_dim1 = *ic;
    c_dim2 = *p;
    c_offset = 1 + c_dim1 * (1 + c_dim2 * 1);
    c__ -= c_offset;
    ey_dim1 = *npred;
    ey_offset = 1 + ey_dim1 * 1;
    ey -= ey_offset;
    y_dim1 = *nacc;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    u_dim1 = *nacc;
    u_offset = 1 + u_dim1 * 1;
    u -= u_offset;
    --ww;
    bb_dim1 = *is;
    bb_offset = 1 + bb_dim1 * 1;
    bb -= bb_offset;
    aa_dim1 = *is;
    aa_offset = 1 + aa_dim1 * 1;
    aa -= aa_offset;

    /* Function Body */
    hperr = 0;
    i__1 = *nh;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L1: */
	nt[i__] = 0;
    }
    i__1 = *nh;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *p;
	    for (j = 1; j <= i__3; ++j) {
/* L2: */
		cov[k + (i__ + j * cov_dim2) * cov_dim1] = 0.;
	    }
	}
    }
    i__3 = *npred + 1 - horiz[1];
    for (it = *dscard; it <= i__3; ++it) {
	arma_(&ey[ey_offset], &hperr, prderr, errwt, m, p, ia, ib, ic, &it, 
		npred, nacc, &u[u_offset], &y[y_offset], &a[a_offset], &b[
		b_offset], &c__[c_offset], &trend[1], is, &aa[aa_offset], &bb[
		bb_offset], &ww[1]);
/*       Eliminate longer horizons as date runs out. */
/*       This assumes HORIZ is sorted in ascending order. */
	if (it - 1 + horiz[*nh] > *npred) {
	    --(*nh);
	}
	i__2 = *nh;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    hi = it - 1 + horiz[i__];
	    i__1 = *p;
	    for (j = 1; j <= i__1; ++j) {
/* L4: */
		ey[i__ + j * ey_dim1] = ey[hi + j * ey_dim1] - y[hi + j * 
			y_dim1];
	    }
	}
	i__1 = *nh;
	for (k = 1; k <= i__1; ++k) {
	    ++nt[k];
	    mf = (doublereal) (nt[k] - 1) / (doublereal) nt[k];
	    i__2 = *p;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__4 = *p;
		for (j = 1; j <= i__4; ++j) {
		    cov[k + (i__ + j * cov_dim2) * cov_dim1] = cov[k + (i__ + 
			    j * cov_dim2) * cov_dim1] * mf + ey[k + i__ * 
			    ey_dim1] * ey[k + j * ey_dim1] / nt[k];
/* L5: */
		}
	    }
	}
/* L10: */
    }
    return 0;
} /* rmaepr_ */

/* Subroutine */ int arma_(ey, hperr, prderr, errwt, m, p, ia, ib, ic, nsmpl, 
	npred, nacc, u, y, a, b, c__, trend, is, aa, bb, ww)
doublereal *ey;
integer *hperr;
doublereal *prderr, *errwt;
integer *m, *p, *ia, *ib, *ic, *nsmpl, *npred, *nacc;
doublereal *u, *y, *a, *b, *c__, *trend;
integer *is;
doublereal *aa, *bb, *ww;
{
    /* System generated locals */
    integer ey_dim1, ey_offset, prderr_dim1, prderr_offset, y_dim1, y_offset, 
	    u_dim1, u_offset, a_dim1, a_dim2, a_offset, b_dim1, b_dim2, 
	    b_offset, c_dim1, c_dim2, c_offset, aa_dim1, aa_offset, bb_dim1, 
	    bb_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal detom;
    static integer lperr, it;
    extern /* Subroutine */ int invers_();

/* sampleT is the length of data which should be used for estimation. */
/* Calculate the one-step ahead predictions, and likelihood value for the model: */

/*       A(L)y(t) =  B(L)w(t) + C(L)u(t)  + TREND */

/* A(L) (axpxp) is the auto-regressive polynomial array. */
/* B(L) (bxpxp) is the moving-average polynomial array. */
/* C(L) (cxpxm) is the  input polynomial array. */
/* TREND is a constant vector added at each period. */
/* y is the p dimensional output data. */
/* u is the m dimensional control (input) data. */
/*  M is the dimension of the input u. */
/*  P is the dimension of the output y and the ouput noise w. */
/*  NSMPL is the length of the data series to use for residual */
/*      and likelihood calculations. */
/*  NPRED is the period to predict (past NSMPL) */
/*  NACC is the actual first (time) dimension of Y and U. */

/*   EY is the output prediction. Initially EY is used to store WW. */
/*   The prediction error WW (innovations) = Y - EY. */
/*   Weighted prediction errors are returned in PRDERR. */
/*     If HPERR is equal or greater than one then weighted prediction */
/*     errors are calculated up to the horizon indicated */
/*     by HPERR. The weights taken from ERRWT are applied to the squared */
/*     error at each period ahead. */
/*   If HPERR is zero or all elements of ERRWT are zero then */
/*   LPERR is set false then PRDERR is not referenced, so ARMA can be called */
/*    with dummy arguments as in GEND. */

/*    NSTART is not properly implemented and must be set to 1. */

/*   IS should be max(P,M) */



/*       CALL DBPRDB('inARMA ',7, 1,1) */
/*      CALL DBPR('M     ',6, M,1) */
/*      CALL DBPR('P     ',6, P,1) */
/*      CALL DBPR('IA    ',6, IA,1) */
/*      CALL DBPR('IB    ',6, IB,1) */
/*      CALL DBPR('IC    ',6, IC,1) */
/*      CALL DBPR('NSMPL  ',6, NSMPL,1) */
/*      CALL DBPR('NPRED  ',6, NPRED,1) */
/*      CALL DBPR('NACC  ',6, NACC,1) */
/*       CALL DBPRDB('A     ',6, A,(IA*P*P)) */
/*       CALL DBPRDB('B     ',6, B,(IB*P*P)) */
/*       CALL DBPRDB('C     ',6, C,(IC*P*M)) */
    /* Parameter adjustments */
    --errwt;
    --trend;
    a_dim1 = *ia;
    a_dim2 = *p;
    a_offset = 1 + a_dim1 * (1 + a_dim2 * 1);
    a -= a_offset;
    b_dim1 = *ib;
    b_dim2 = *p;
    b_offset = 1 + b_dim1 * (1 + b_dim2 * 1);
    b -= b_offset;
    c_dim1 = *ic;
    c_dim2 = *p;
    c_offset = 1 + c_dim1 * (1 + c_dim2 * 1);
    c__ -= c_offset;
    prderr_dim1 = *nsmpl;
    prderr_offset = 1 + prderr_dim1 * 1;
    prderr -= prderr_offset;
    ey_dim1 = *npred;
    ey_offset = 1 + ey_dim1 * 1;
    ey -= ey_offset;
    y_dim1 = *nacc;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    u_dim1 = *nacc;
    u_offset = 1 + u_dim1 * 1;
    u -= u_offset;
    --ww;
    bb_dim1 = *is;
    bb_offset = 1 + bb_dim1 * 1;
    bb -= bb_offset;
    aa_dim1 = *is;
    aa_offset = 1 + aa_dim1 * 1;
    aa -= aa_offset;

    /* Function Body */
    lperr = 0;
    if (*hperr > 0) {
	i__1 = *hperr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (errwt[i__] > (float)0.) {
		lperr = 1;
	    }
/* L500: */
	}
    }

/*     Ensure B(0) = I by inverting B(0) and multiplying through */
/*        B(1,,) is not modified yet as B(0) is needed later, but */
/*        it is not referenced through the time loop, so effectively */
/*        asummed = I. */

    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *p;
	for (j = 1; j <= i__2; ++j) {
/* L300: */
	    bb[i__ + j * bb_dim1] = b[(i__ + j * b_dim2) * b_dim1 + 1];
	}
    }
    invers_(&bb[bb_offset], p, is, &detom);
    i__2 = *ia;
    for (l = 1; l <= i__2; ++l) {
	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = *p;
	    for (j = 1; j <= i__3; ++j) {
		aa[i__ + j * aa_dim1] = a[l + (i__ + j * a_dim2) * a_dim1];
/* L301: */
		a[l + (i__ + j * a_dim2) * a_dim1] = 0.;
	    }
	}
	i__3 = *p;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__1 = *p;
	    for (j = 1; j <= i__1; ++j) {
		i__4 = *p;
		for (k = 1; k <= i__4; ++k) {
/* L302: */
		    a[l + (i__ + j * a_dim2) * a_dim1] += bb[i__ + k * 
			    bb_dim1] * aa[k + j * aa_dim1];
		}
	    }
	}
    }
    i__4 = *ib;
    for (l = 2; l <= i__4; ++l) {
	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = *p;
	    for (j = 1; j <= i__3; ++j) {
		aa[i__ + j * aa_dim1] = b[l + (i__ + j * b_dim2) * b_dim1];
/* L303: */
		b[l + (i__ + j * b_dim2) * b_dim1] = 0.;
	    }
	}
	i__3 = *p;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__1 = *p;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *p;
		for (k = 1; k <= i__2; ++k) {
/* L304: */
		    b[l + (i__ + j * b_dim2) * b_dim1] += bb[i__ + k * 
			    bb_dim1] * aa[k + j * aa_dim1];
		}
	    }
	}
    }
    i__2 = *ic;
    for (l = 1; l <= i__2; ++l) {
	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = *m;
	    for (j = 1; j <= i__3; ++j) {
		aa[i__ + j * aa_dim1] = c__[l + (i__ + j * c_dim2) * c_dim1];
/* L305: */
		c__[l + (i__ + j * c_dim2) * c_dim1] = 0.;
	    }
	}
	i__3 = *p;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		i__4 = *p;
		for (k = 1; k <= i__4; ++k) {
/* L306: */
		    c__[l + (i__ + j * c_dim2) * c_dim1] += bb[i__ + k * 
			    bb_dim1] * aa[k + j * aa_dim1];
		}
	    }
	}
    }
    i__4 = *p;
    for (i__ = 1; i__ <= i__4; ++i__) {
	aa[i__ + aa_dim1] = trend[i__];
/* L307: */
	trend[i__] = 0.;
    }
    i__4 = *p;
    for (i__ = 1; i__ <= i__4; ++i__) {
	i__1 = *p;
	for (k = 1; k <= i__1; ++k) {
/* L308: */
	    trend[i__] += bb[i__ + k * bb_dim1] * aa[i__ + aa_dim1];
	}
    }
    i__1 = *npred;
    for (it = 1; it <= i__1; ++it) {
	i__4 = *p;
	for (j = 1; j <= i__4; ++j) {
/* L309: */
	    ey[it + j * ey_dim1] = 0.;
	}
    }

/*    Start Time loop */

    i__4 = *nsmpl;
    for (it = 1; it <= i__4; ++it) {

	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L1: */
	    ww[i__] = -trend[i__];
	}

	i__1 = *ia;
	for (l = 1; l <= i__1; ++l) {
	    if (l <= it) {
		i__3 = *p;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__2 = *p;
		    for (j = 1; j <= i__2; ++j) {
/* L2: */
			ww[i__] += a[l + (i__ + j * a_dim2) * a_dim1] * y[it 
				+ 1 - l + j * y_dim1];
		    }
		}
	    }
/* L22: */
	}

	if (*ib >= 2) {
	    i__1 = *ib;
	    for (l = 2; l <= i__1; ++l) {
		if (l <= it) {
		    i__2 = *p;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = *p;
			for (j = 1; j <= i__3; ++j) {
/* L3: */
			    ww[i__] -= b[l + (i__ + j * b_dim2) * b_dim1] * 
				    ey[it + 1 - l + j * ey_dim1];
			}
		    }
		}
/* L23: */
	    }
	}

	i__1 = *ic;
	for (l = 1; l <= i__1; ++l) {
	    if (l <= it) {
		i__3 = *p;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__2 = *m;
		    for (j = 1; j <= i__2; ++j) {
/*         CALL DBPRDB('C     ',6, C(L,I,J),1) */
/*         CALL DBPRDB('U     ',6,U(IT+1-L,J) ,1) */
/* L4: */
			ww[i__] -= c__[l + (i__ + j * c_dim2) * c_dim1] * u[
				it + 1 - l + j * u_dim1];
		    }
		}
	    }
/* L24: */
	}
/*      IF (IT.LE.3) THEN */
/*         CALL DBPRDB('ww  ',4, WW,P) */
/*      ENDIF */

/*   EY stores history of predition error WW to reconstruct predictions */
	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	    ey[it + i__ * ey_dim1] = ww[i__];
	}
/*   Return weighted prediction error */
	if (lperr == 1) {
	    i__1 = *p;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L410: */
/* Computing 2nd power */
		d__1 = ww[i__];
		prderr[it + i__ * prderr_dim1] = errwt[1] * (d__1 * d__1);
	    }
	    if (*hperr > 1) {
		i__1 = *hperr;
		for (k = 2; k <= i__1; ++k) {
		    if (it + k - 1 <= *nsmpl) {
			i__2 = *p;
			for (i__ = 1; i__ <= i__2; ++i__) {
/* L401: */
			    ww[i__] = -trend[i__];
			}

			i__2 = *ia;
			for (l = 1; l <= i__2; ++l) {
			    if (l < it + k) {
				i__3 = *p;
				for (i__ = 1; i__ <= i__3; ++i__) {
				    i__5 = *p;
				    for (j = 1; j <= i__5; ++j) {
/* L402: */
					ww[i__] += a[l + (i__ + j * a_dim2) * 
						a_dim1] * y[it + k - l + j * 
						y_dim1];
				    }
				}
			    }
/* L4022: */
			}

			if (*ib >= 2) {
			    i__2 = *ib;
			    for (l = 2; l <= i__2; ++l) {
				if (l < it + k) {
				    i__5 = *p;
				    for (i__ = 1; i__ <= i__5; ++i__) {
					i__3 = *p;
					for (j = 1; j <= i__3; ++j) {
/* L403: */
					    ww[i__] -= b[l + (i__ + j * 
						    b_dim2) * b_dim1] * ey[it 
						    + k - l + j * ey_dim1];
					}
				    }
				}
/* L4023: */
			    }
			}

			i__2 = *ic;
			for (l = 1; l <= i__2; ++l) {
			    if (l < it + k) {
				i__3 = *p;
				for (i__ = 1; i__ <= i__3; ++i__) {
				    i__5 = *m;
				    for (j = 1; j <= i__5; ++j) {
/* L404: */
					ww[i__] -= c__[l + (i__ + j * c_dim2) 
						* c_dim1] * u[it + k - l + j *
						 u_dim1];
				    }
				}
			    }
/* L4024: */
			}
/*              correction for WW by B0. */
/*      NB this has not been well tested and models with B0 != I may not work !! */
			i__2 = *p;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    i__5 = *p;
			    for (j = 1; j <= i__5; ++j) {
/* L407: */
				bb[i__ + bb_dim1] += b[(i__ + j * b_dim2) * 
					b_dim1 + 1] * ww[j];
			    }
			}
			i__5 = *p;
			for (i__ = 1; i__ <= i__5; ++i__) {
/* L408: */
/* Computing 2nd power */
			    d__1 = bb[i__ + bb_dim1];
			    prderr[it + i__ * prderr_dim1] += errwt[k] * (
				    d__1 * d__1);
			}
/*                 IF (IT.GT.97) THEN */
/*                    CALL DBPR('K     ',6, K,1) */
/*                    CALL DBPRDB('PRDERR ',8,PRDERR,NSMPL*P) */
/*                 ENDIF */
		    }
/* L400: */
		}
	    }
	}
/* L1000: */
    }
/*  end of time loop */

    i__4 = *nsmpl;
    for (it = 1; it <= i__4; ++it) {
	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	    ww[i__] = (float)0.;
	}
	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__5 = *p;
	    for (j = 1; j <= i__5; ++j) {
/* L41: */
		ww[i__] += b[(i__ + j * b_dim2) * b_dim1 + 1] * ey[it + j * 
			ey_dim1];
	    }
	}
	i__5 = *p;
	for (i__ = 1; i__ <= i__5; ++i__) {
/* L45: */
	    ey[it + i__ * ey_dim1] = y[it + i__ * y_dim1] - ww[i__];
	}
    }

/*    Start multi-step prediction loop */

    if (*nsmpl < *npred) {
/*  B(1,,) now needs to be filled in as I (storage was previously use for B0). */
	i__5 = *p;
	for (i__ = 1; i__ <= i__5; ++i__) {
	    i__4 = *p;
	    for (j = 1; j <= i__4; ++j) {
/* L2298: */
		b[(i__ + j * b_dim2) * b_dim1 + 1] = (float)0.;
	    }
	}
	i__4 = *p;
	for (i__ = 1; i__ <= i__4; ++i__) {
/* L2299: */
	    b[(i__ + i__ * b_dim2) * b_dim1 + 1] = (float)1.;
	}

/*     Ensure A(0) = I by inverting A(0) and multiplying through */

	i__4 = *p;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    i__5 = *p;
	    for (j = 1; j <= i__5; ++j) {
/* L2300: */
		bb[i__ + j * bb_dim1] = a[(i__ + j * a_dim2) * a_dim1 + 1];
	    }
	}
	invers_(&bb[bb_offset], p, is, &detom);
	i__5 = *ia;
	for (l = 1; l <= i__5; ++l) {
	    i__4 = *p;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__1 = *p;
		for (j = 1; j <= i__1; ++j) {
		    aa[i__ + j * aa_dim1] = a[l + (i__ + j * a_dim2) * a_dim1]
			    ;
/* L2301: */
		    a[l + (i__ + j * a_dim2) * a_dim1] = 0.;
		}
	    }
	    i__1 = *p;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__4 = *p;
		for (j = 1; j <= i__4; ++j) {
		    i__2 = *p;
		    for (k = 1; k <= i__2; ++k) {
/* L2302: */
			a[l + (i__ + j * a_dim2) * a_dim1] += bb[i__ + k * 
				bb_dim1] * aa[k + j * aa_dim1];
		    }
		}
	    }
	}
	i__2 = *ib;
	for (l = 1; l <= i__2; ++l) {
	    i__4 = *p;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__1 = *p;
		for (j = 1; j <= i__1; ++j) {
		    aa[i__ + j * aa_dim1] = b[l + (i__ + j * b_dim2) * b_dim1]
			    ;
/* L2303: */
		    b[l + (i__ + j * b_dim2) * b_dim1] = 0.;
		}
	    }
	    i__1 = *p;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__4 = *p;
		for (j = 1; j <= i__4; ++j) {
		    i__5 = *p;
		    for (k = 1; k <= i__5; ++k) {
/* L2304: */
			b[l + (i__ + j * b_dim2) * b_dim1] += bb[i__ + k * 
				bb_dim1] * aa[k + j * aa_dim1];
		    }
		}
	    }
	}
	i__5 = *ic;
	for (l = 1; l <= i__5; ++l) {
	    i__4 = *p;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__1 = *m;
		for (j = 1; j <= i__1; ++j) {
		    aa[i__ + j * aa_dim1] = c__[l + (i__ + j * c_dim2) * 
			    c_dim1];
/* L2305: */
		    c__[l + (i__ + j * c_dim2) * c_dim1] = 0.;
		}
	    }
	    i__1 = *p;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__4 = *m;
		for (j = 1; j <= i__4; ++j) {
		    i__2 = *p;
		    for (k = 1; k <= i__2; ++k) {
/* L2306: */
			c__[l + (i__ + j * c_dim2) * c_dim1] += bb[i__ + k * 
				bb_dim1] * aa[k + j * aa_dim1];
		    }
		}
	    }
	}
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    aa[i__ + aa_dim1] = trend[i__];
/* L2307: */
	    trend[i__] = 0.;
	}
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__4 = *p;
	    for (k = 1; k <= i__4; ++k) {
/* L2308: */
		trend[i__] += bb[i__ + k * bb_dim1] * aa[i__ + aa_dim1];
	    }
	}
	i__4 = *npred;
	for (it = *nsmpl + 1; it <= i__4; ++it) {

	    i__2 = *p;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L2001: */
		ey[it + i__ * ey_dim1] = trend[i__];
	    }

	    i__2 = *ia;
	    for (l = 2; l <= i__2; ++l) {
		i__1 = *p;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__5 = *p;
		    for (j = 1; j <= i__5; ++j) {
			if (it + 1 - l <= *nsmpl) {
			    ey[it + i__ * ey_dim1] -= a[l + (i__ + j * a_dim2)
				     * a_dim1] * y[it + 1 - l + j * y_dim1];
			} else {
			    ey[it + i__ * ey_dim1] -= a[l + (i__ + j * a_dim2)
				     * a_dim1] * ey[it + 1 - l + j * ey_dim1];
			}
/* L2002: */
		    }
		}
	    }
	    if (*ib >= 2) {
		i__5 = *ib;
		for (l = 2; l <= i__5; ++l) {
		    if (it + 1 - l <= *nsmpl) {
			i__1 = *p;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    i__2 = *p;
			    for (j = 1; j <= i__2; ++j) {
				ey[it + i__ * ey_dim1] += b[l + (i__ + j * 
					b_dim2) * b_dim1] * (y[it + 1 - l + j 
					* y_dim1] - ey[it + 1 - l + j * 
					ey_dim1]);
/* L2003: */
			    }
			}
		    }
/* L2004: */
		}
	    }

	    i__5 = *ic;
	    for (l = 1; l <= i__5; ++l) {
		i__2 = *p;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__1 = *m;
		    for (j = 1; j <= i__1; ++j) {
			ey[it + i__ * ey_dim1] += c__[l + (i__ + j * c_dim2) *
				 c_dim1] * u[it + 1 - l + j * u_dim1];
/* L2005: */
		    }
		}
	    }
/* L2000: */
	}
    }
/*  end of multi-step prediction loop */
    return 0;
} /* arma_ */

/* Subroutine */ int datepr_(cov, dscard, horiz, nh, nt, p, npred, err)
doublereal *cov;
integer *dscard, *horiz, *nh, *nt, *p, *npred;
doublereal *err;
{
    /* System generated locals */
    integer cov_dim1, cov_dim2, cov_offset, err_dim1, err_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, j, k, hi;
    static doublereal mf;
    static integer it;

/*     See S program predictions.cov.TSdata */


/*        CALL DBPR('NPRED ',6, NPRED,1) */
/*        CALL DBPR('HORIZ ',6, HORIZ(1),1) */
/*        CALL DBPR('DSCARD',7, DSCARD,1) */
    /* Parameter adjustments */
    --nt;
    --horiz;
    cov_dim1 = *nh;
    cov_dim2 = *p;
    cov_offset = 1 + cov_dim1 * (1 + cov_dim2 * 1);
    cov -= cov_offset;
    err_dim1 = *npred;
    err_offset = 1 + err_dim1 * 1;
    err -= err_offset;

    /* Function Body */
    i__1 = *nh;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L1: */
	nt[i__] = 0;
    }
    i__1 = *nh;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *p;
	    for (j = 1; j <= i__3; ++j) {
/* L2: */
		cov[k + (i__ + j * cov_dim2) * cov_dim1] = 0.;
	    }
	}
    }
    i__3 = *npred + 1 - horiz[1];
    for (it = *dscard; it <= i__3; ++it) {
/*       this assumes HORIZ is sorted in ascending order */
	if (it - 1 + horiz[*nh] > *npred) {
	    --(*nh);
	}
	i__2 = *nh;
	for (k = 1; k <= i__2; ++k) {
	    ++nt[k];
	    mf = (doublereal) (nt[k] - 1) / (doublereal) nt[k];
	    hi = it - 1 + horiz[k];
	    i__1 = *p;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__4 = *p;
		for (j = 1; j <= i__4; ++j) {
		    cov[k + (i__ + j * cov_dim2) * cov_dim1] = cov[k + (i__ + 
			    j * cov_dim2) * cov_dim1] * mf + err[hi + i__ * 
			    err_dim1] * err[hi + j * err_dim1] / nt[k];
/* L5: */
		}
	    }
	}
/* L10: */
    }
    return 0;
} /* datepr_ */

/* Subroutine */ int invers_(a, n, is, det)
doublereal *a;
integer *n, *is;
doublereal *det;
{
    /* Initialized data */

    static doublereal eps = 1e-20;
    static doublereal epss = 1e-100;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal work[1000];
    static integer i__, j, k, m;
    static doublereal w, y;
    static integer iwork[1000], jn;

/*    IS is the declared dimension of A and */
/*      N is the dimension of the matrix to invert. */
/*     COMMON /MAT/ A */
/*     COMMON /SIZE/ N */
/*     BOC_INVMAT is a real matrix inversion routine.It invert the */
/*     array A(N,N) in place; that is, the array A is destroyed by */
/*     the  routine  and  the inverse takes its place. */
    /* Parameter adjustments */
    a_dim1 = *is;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    *det = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	iwork[j - 1] = j;
/* L100: */
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = i__;
	y = a[i__ + i__ * a_dim1];
	m = i__ + 1;
	if (i__ == *n) {
	    goto L120;
	}
	i__2 = *n;
	for (j = m; j <= i__2; ++j) {
	    w = a[i__ + j * a_dim1];
	    if (abs(w) <= abs(y)) {
		goto L110;
	    }
	    k = j;
	    y = w;
L110:
	    ;
	}
L120:
	*det *= y;
	if (abs(y) < eps) {
	    goto L900;
	}
	if (abs(*det) < epss) {
	    goto L900;
	}
	y = (float)1. / y;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    work[j - 1] = a[j + k * a_dim1];
	    a[j + k * a_dim1] = a[j + i__ * a_dim1];
	    a[j + i__ * a_dim1] = -work[j - 1] * y;
	    jn = j + *n;
	    a[i__ + j * a_dim1] *= y;
	    work[jn - 1] = a[i__ + j * a_dim1];
/* L130: */
	}
	a[i__ + i__ * a_dim1] = y;
	j = iwork[i__ - 1];
	iwork[i__ - 1] = iwork[k - 1];
	iwork[k - 1] = j;
	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		if (k == i__ || j == i__) {
		    goto L140;
		}
		jn = j + *n;
		a[k + j * a_dim1] -= work[jn - 1] * work[k - 1];
L140:
		;
	    }
/* L150: */
	}
/* L160: */
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
L170:
	k = iwork[i__ - 1];
	if (k == i__) {
	    goto L190;
	}
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    w = a[i__ + j * a_dim1];
	    a[i__ + j * a_dim1] = a[k + j * a_dim1];
	    a[k + j * a_dim1] = w;
/* L180: */
	}
	m = iwork[i__ - 1];
	iwork[i__ - 1] = iwork[k - 1];
	iwork[k - 1] = m;
	*det = -(*det);
	goto L170;
L190:
	;
    }
/*      CALL DBPRDB('DET ',4, DET,1) */
    return 0;
L900:
    *det = (float)0.;
    return 0;
} /* invers_ */

/* routines for curvature calculation */

/* Subroutine */ int gend_(d__, fc, ith, x0, delta0, n, nd, f0, rd, haprox, 
	hdiag, daprox, x, delta, f1, f2, m, p, nsmpl, nacc, u, y, ap, ip, jp, 
	ict, const__, an, in, jn, lp, ln, ia, ib, ic, a, b, c__, ns, z0, p0, 
	f, g, h__, fk, q, r__, gain)
doublereal *d__;
integer *fc, *ith;
doublereal *x0, *delta0;
integer *n, *nd;
doublereal *f0;
integer *rd;
doublereal *haprox, *hdiag, *daprox, *x, *delta, *f1, *f2;
integer *m, *p, *nsmpl, *nacc;
doublereal *u, *y;
integer *ap, *ip, *jp, *ict;
doublereal *const__;
integer *an, *in, *jn, *lp, *ln, *ia, *ib, *ic;
doublereal *a, *b, *c__;
integer *ns;
doublereal *z0, *p0, *f, *g, *h__, *fk, *q, *r__;
integer *gain;
{
    /* System generated locals */
    integer d_dim1, d_offset, daprox_dim1, daprox_offset, hdiag_dim1, 
	    hdiag_offset, haprox_dim1, haprox_offset, y_dim1, y_offset, 
	    u_dim1, u_offset, a_dim1, a_dim2, a_offset, b_dim1, b_dim2, 
	    b_offset, c_dim1, c_dim2, c_offset, p0_dim1, p0_offset, f_dim1, 
	    f_offset, g_dim1, g_offset, h_dim1, h_offset, fk_dim1, fk_offset, 
	    q_dim1, q_offset, r_dim1, r_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_di();

    /* Local variables */
    extern /* Subroutine */ int dbpr_();
    static integer i__, j, k;
    static doublereal v;
    extern /* Subroutine */ int armap_();
    static integer hperr;
    extern /* Subroutine */ int error_();
    static doublereal errwt[1];
    static integer mc, ii;
    static doublereal md;
    static integer up;
    static doublereal prderr[1]	/* was [1][1] */;
    extern /* Subroutine */ int kfp_();


/*  Z0 is TREND for ARMA models. */

/*  FC  indicator of function (0=KF, 1=ARMA). It would be nice if this */
/*      could be the  name of the function as in C. */
/*      The function must have a single vector arguement X. */
/*  X0   the parameter vector. */
/*  X    is the working copy (altered by DELTA). */
/*  ITH   is the length of the parameter vector. */
/*  F0  is the value (in sample/residual space) of the function. */
/*      (only the space is needed, the function is calculated). */
/*  N   is the dimension of the sample space (length of F0). */
/*  DELTA0  gives the fraction of X to use for the initial */
/*           numerical approximation. */
/*  ND  is the number of columns of matrix D.( first */
/*               der. & lower triangle of Hessian) */
/*  EPS     is used for zero elements of X. */
/*  RD       the number of Richardson improvement iterations. */
/*  V=2       reduction factor for Richardson iterations. */
/*       V could be a parameter but the way the reduction formula is */
/*        coded assumes it is =2 */



/*  parameters passed directly to ARMAp and/or KFp: */


/*      F, Q, and R are used for scratch space in the call to ARMAP instead of: */
/*      DOUBLE PRECISION AA(IS,IS), BB(IS,IS), WW(IS) */
/*        this could cause some problems ... some checks are made. */
/*       Z0(NS) is used for TREND(P) in ARMA models */
/*     PARM(ITH) for ARMAP/KFP is X(ITH) in GEND, EY is function value */

/*      DOUBLE PRECISION Z(NSMPL,NS),TRKERR(NSMPL,NS,NS) */

/* ..bug in S: passing characters is unreliable */
/*   use integer for AP and AN... */
    /* Parameter adjustments */
    --lp;
    --jp;
    --ip;
    --ap;
    --delta;
    --x;
    --delta0;
    --x0;
    --f2;
    --f1;
    hdiag_dim1 = *n;
    hdiag_offset = 1 + hdiag_dim1 * 1;
    hdiag -= hdiag_offset;
    --f0;
    d_dim1 = *n;
    d_offset = 1 + d_dim1 * 1;
    d__ -= d_offset;
    daprox_dim1 = *n;
    daprox_offset = 1 + daprox_dim1 * 1;
    daprox -= daprox_offset;
    haprox_dim1 = *n;
    haprox_offset = 1 + haprox_dim1 * 1;
    haprox -= haprox_offset;
    r_dim1 = *p;
    r_offset = 1 + r_dim1 * 1;
    r__ -= r_offset;
    y_dim1 = *nacc;
    y_offset = 1 + y_dim1 * 1;
    y -= y_offset;
    u_dim1 = *nacc;
    u_offset = 1 + u_dim1 * 1;
    u -= u_offset;
    --ln;
    --jn;
    --in;
    --an;
    --const__;
    a_dim1 = *ia;
    a_dim2 = *p;
    a_offset = 1 + a_dim1 * (1 + a_dim2 * 1);
    a -= a_offset;
    b_dim1 = *ib;
    b_dim2 = *p;
    b_offset = 1 + b_dim1 * (1 + b_dim2 * 1);
    b -= b_offset;
    c_dim1 = *ic;
    c_dim2 = *p;
    c_offset = 1 + c_dim1 * (1 + c_dim2 * 1);
    c__ -= c_offset;
    q_dim1 = *ns;
    q_offset = 1 + q_dim1 * 1;
    q -= q_offset;
    fk_dim1 = *ns;
    fk_offset = 1 + fk_dim1 * 1;
    fk -= fk_offset;
    h_dim1 = *p;
    h_offset = 1 + h_dim1 * 1;
    h__ -= h_offset;
    g_dim1 = *ns;
    g_offset = 1 + g_dim1 * 1;
    g -= g_offset;
    f_dim1 = *ns;
    f_offset = 1 + f_dim1 * 1;
    f -= f_offset;
    p0_dim1 = *ns;
    p0_offset = 1 + p0_dim1 * 1;
    p0 -= p0_offset;
    --z0;

    /* Function Body */
    dbpr_("starting gend N=", &c__16, n, &c__1, (ftnlen)16);
    dbpr_("            ITH=", &c__16, ith, &c__1, (ftnlen)16);
    dbpr_("             ND=", &c__16, nd, &c__1, (ftnlen)16);
    if (*ns < max(*m,*p)) {
	error_("warning: scratch (NS too small) in GEND.", &c__40, ns, &c__1, 
		(ftnlen)40);
    }
    if (*ns < *p * *p) {
	error_("warning: scratch (P too big) in GEND.", &c__37, ns, &c__1, (
		ftnlen)37);
    }
    hperr = 0;
    v = (float)2.;
    i__1 = *ith;
    for (ii = 1; ii <= i__1; ++ii) {
/* L1: */
	x[ii] = x0[ii];
    }
    if (*fc == 0) {
	kfp_(&f0[1], &hperr, prderr, errwt, m, ns, p, nsmpl, nsmpl, nacc, &u[
		u_offset], &y[y_offset], &f[f_offset], &g[g_offset], &h__[
		h_offset], &fk[fk_offset], &q[q_offset], &r__[r_offset], gain,
		 &z0[1], &p0[p0_offset], ith, &x[1], &ap[1], &ip[1], &jp[1], 
		ict, &const__[1], &an[1], &in[1], &jn[1]);
    } else if (*fc == 1) {
	armap_(&f0[1], &hperr, prderr, errwt, m, p, ia, ib, ic, nsmpl, nsmpl, 
		nacc, &u[u_offset], &y[y_offset], &a[a_offset], &b[b_offset], 
		&c__[c_offset], &z0[1], ith, &x[1], &ap[1], &lp[1], &ip[1], &
		jp[1], ict, &const__[1], &an[1], &ln[1], &in[1], &jn[1], ns, &
		f[f_offset], &q[q_offset], &r__[r_offset]);
    }
/*                   each parameter  - first deriv. & hessian diagonal */
    i__1 = *ith;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ith;
	for (ii = 1; ii <= i__2; ++ii) {
/* L10: */
	    delta[ii] = delta0[ii];
	}
/*                                 successively reduce DELTA */

/*  This could be done without both X and X0 by adding and then subtracting */
/*   DELTA, but accumulated round off error seems to affect the result. */
	i__2 = *rd;
	for (k = 1; k <= i__2; ++k) {
	    x[i__] = x0[i__] + delta[i__];
	    if (*fc == 0) {
		kfp_(&f1[1], &hperr, prderr, errwt, m, ns, p, nsmpl, nsmpl, 
			nacc, &u[u_offset], &y[y_offset], &f[f_offset], &g[
			g_offset], &h__[h_offset], &fk[fk_offset], &q[
			q_offset], &r__[r_offset], gain, &z0[1], &p0[
			p0_offset], ith, &x[1], &ap[1], &ip[1], &jp[1], ict, &
			const__[1], &an[1], &in[1], &jn[1]);
	    } else if (*fc == 1) {
		armap_(&f1[1], &hperr, prderr, errwt, m, p, ia, ib, ic, nsmpl,
			 nsmpl, nacc, &u[u_offset], &y[y_offset], &a[a_offset]
			, &b[b_offset], &c__[c_offset], &z0[1], ith, &x[1], &
			ap[1], &lp[1], &ip[1], &jp[1], ict, &const__[1], &an[
			1], &ln[1], &in[1], &jn[1], ns, &f[f_offset], &q[
			q_offset], &r__[r_offset]);
	    }
	    x[i__] = x0[i__] - delta[i__];
	    if (*fc == 0) {
		kfp_(&f2[1], &hperr, prderr, errwt, m, ns, p, nsmpl, nsmpl, 
			nacc, &u[u_offset], &y[y_offset], &f[f_offset], &g[
			g_offset], &h__[h_offset], &fk[fk_offset], &q[
			q_offset], &r__[r_offset], gain, &z0[1], &p0[
			p0_offset], ith, &x[1], &ap[1], &ip[1], &jp[1], ict, &
			const__[1], &an[1], &in[1], &jn[1]);
	    } else if (*fc == 1) {
		armap_(&f2[1], &hperr, prderr, errwt, m, p, ia, ib, ic, nsmpl,
			 nsmpl, nacc, &u[u_offset], &y[y_offset], &a[a_offset]
			, &b[b_offset], &c__[c_offset], &z0[1], ith, &x[1], &
			ap[1], &lp[1], &ip[1], &jp[1], ict, &const__[1], &an[
			1], &ln[1], &in[1], &jn[1], ns, &f[f_offset], &q[
			q_offset], &r__[r_offset]);
	    }
	    x[i__] = x0[i__];
	    i__3 = *n;
	    for (ii = 1; ii <= i__3; ++ii) {
/* L15: */
		daprox[ii + k * daprox_dim1] = (f1[ii] - f2[ii]) / (delta[i__]
			 * (float)2.);
	    }
	    i__3 = *n;
	    for (ii = 1; ii <= i__3; ++ii) {
/* L16: */
/* Computing 2nd power */
		d__1 = delta[i__];
		haprox[ii + k * haprox_dim1] = (f1[ii] - f0[ii] * (float)2. + 
			f2[ii]) / (d__1 * d__1);
	    }
	    delta[i__] /= v;
/* L20: */
	}
	i__2 = *rd - 1;
	for (mc = 1; mc <= i__2; ++mc) {
	    md = pow_di(&c_b236, &mc);
	    i__3 = *rd - mc;
	    for (k = 1; k <= i__3; ++k) {
		i__4 = *n;
		for (ii = 1; ii <= i__4; ++ii) {
/* L25: */
		    daprox[ii + k * daprox_dim1] = (daprox[ii + (k + 1) * 
			    daprox_dim1] * md - daprox[ii + k * daprox_dim1]) 
			    / (md - 1);
		}
		i__4 = *n;
		for (ii = 1; ii <= i__4; ++ii) {
/* L26: */
		    haprox[ii + k * haprox_dim1] = (haprox[ii + (k + 1) * 
			    haprox_dim1] * md - haprox[ii + k * haprox_dim1]) 
			    / (md - 1);
		}
/* L30: */
	    }
	}
	i__3 = *n;
	for (ii = 1; ii <= i__3; ++ii) {
/* L31: */
	    d__[ii + i__ * d_dim1] = daprox[ii + daprox_dim1];
	}
	i__3 = *n;
	for (ii = 1; ii <= i__3; ++ii) {
/* L32: */
	    hdiag[ii + i__ * hdiag_dim1] = haprox[ii + haprox_dim1];
	}
/* L100: */
    }

/*                  2nd derivative  - do lower half of hessian only */
    up = *ith;
    dbpr_("2nd deriv. UP=\n", &c__16, &up, &c__1, (ftnlen)15);
    i__1 = *ith;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__3 = i__;
	for (j = 1; j <= i__3; ++j) {
	    ++up;
	    dbpr_("      UP=\n", &c__11, &up, &c__1, (ftnlen)10);
	    if (i__ == j) {
		i__2 = *n;
		for (ii = 1; ii <= i__2; ++ii) {
/* L120: */
		    d__[ii + up * d_dim1] = hdiag[ii + i__ * hdiag_dim1];
		}
	    } else {
		i__2 = *ith;
		for (ii = 1; ii <= i__2; ++ii) {
/* L121: */
		    delta[ii] = delta0[ii];
		}
/*                                successively reduce DELTA */
		i__2 = *rd;
		for (k = 1; k <= i__2; ++k) {
		    x[i__] = x0[i__] + delta[i__];
		    x[j] = x0[j] + delta[j];
		    if (*fc == 0) {
			kfp_(&f1[1], &hperr, prderr, errwt, m, ns, p, nsmpl, 
				nsmpl, nacc, &u[u_offset], &y[y_offset], &f[
				f_offset], &g[g_offset], &h__[h_offset], &fk[
				fk_offset], &q[q_offset], &r__[r_offset], 
				gain, &z0[1], &p0[p0_offset], ith, &x[1], &ap[
				1], &ip[1], &jp[1], ict, &const__[1], &an[1], 
				&in[1], &jn[1]);
		    } else if (*fc == 1) {
			dbpr_("calling armap.F1..K=\n", &c__22, &k, &c__1, (
				ftnlen)21);
			armap_(&f1[1], &hperr, prderr, errwt, m, p, ia, ib, 
				ic, nsmpl, nsmpl, nacc, &u[u_offset], &y[
				y_offset], &a[a_offset], &b[b_offset], &c__[
				c_offset], &z0[1], ith, &x[1], &ap[1], &lp[1],
				 &ip[1], &jp[1], ict, &const__[1], &an[1], &
				ln[1], &in[1], &jn[1], ns, &f[f_offset], &q[
				q_offset], &r__[r_offset]);
		    }
		    x[i__] = x0[i__] - delta[i__];
		    x[j] = x0[j] - delta[j];
		    if (*fc == 0) {
			kfp_(&f2[1], &hperr, prderr, errwt, m, ns, p, nsmpl, 
				nsmpl, nacc, &u[u_offset], &y[y_offset], &f[
				f_offset], &g[g_offset], &h__[h_offset], &fk[
				fk_offset], &q[q_offset], &r__[r_offset], 
				gain, &z0[1], &p0[p0_offset], ith, &x[1], &ap[
				1], &ip[1], &jp[1], ict, &const__[1], &an[1], 
				&in[1], &jn[1]);
		    } else if (*fc == 1) {
			armap_(&f2[1], &hperr, prderr, errwt, m, p, ia, ib, 
				ic, nsmpl, nsmpl, nacc, &u[u_offset], &y[
				y_offset], &a[a_offset], &b[b_offset], &c__[
				c_offset], &z0[1], ith, &x[1], &ap[1], &lp[1],
				 &ip[1], &jp[1], ict, &const__[1], &an[1], &
				ln[1], &in[1], &jn[1], ns, &f[f_offset], &q[
				q_offset], &r__[r_offset]);
		    }
		    x[i__] = x0[i__];
		    x[j] = x0[j];
		    i__4 = *n;
		    for (ii = 1; ii <= i__4; ++ii) {
/* L130: */
/* Computing 2nd power */
			d__1 = delta[i__];
/* Computing 2nd power */
			d__2 = delta[j];
			daprox[ii + k * daprox_dim1] = (f1[ii] - f0[ii] * (
				float)2. + f2[ii] - hdiag[ii + i__ * 
				hdiag_dim1] * (d__1 * d__1) - hdiag[ii + j * 
				hdiag_dim1] * (d__2 * d__2)) / (delta[i__] * (
				float)2. * delta[j]);
		    }
		    i__4 = *ith;
		    for (ii = 1; ii <= i__4; ++ii) {
/* L140: */
			delta[ii] /= v;
		    }
/* L150: */
		}
		i__2 = *rd - 1;
		for (mc = 1; mc <= i__2; ++mc) {
		    md = pow_di(&c_b236, &mc);
		    i__4 = *rd - mc;
		    for (k = 1; k <= i__4; ++k) {
			i__5 = *n;
			for (ii = 1; ii <= i__5; ++ii) {
/* L170: */
			    daprox[ii + k * daprox_dim1] = (daprox[ii + (k + 
				    1) * daprox_dim1] * md - daprox[ii + k * 
				    daprox_dim1]) / (md - (float)1.);
			}
		    }
		    i__5 = *n;
		    for (ii = 1; ii <= i__5; ++ii) {
/* L180: */
			d__[ii + up * d_dim1] = daprox[ii + daprox_dim1];
		    }
/* L190: */
		}
	    }
/* L200: */
	}
    }
/*      DBPRDB('gend returning D[1,1]=',24, D(1,1),1) */
    return 0;
} /* gend_ */

#ifdef uNdEfInEd
comments from the converter:  (stderr from f2c)
   error:
   dbpr:
   dbprdb:
   simss:
   smooth:
   kfp:
   kfprj:
   kfepr:
   kf:
   simrma:
   armap:
   rmaprj:
   rmaepr:
   arma:
   datepr:
   invers:
   gend:
#endif
