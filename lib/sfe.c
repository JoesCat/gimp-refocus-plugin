/* sequential formatted external common routines*/
#include "f2c.h"
#include "fio.h"

extern char *f__fmtbuf;

integer e_rsfe(Void)
{	int n;
	n=en_fio();
	f__fmtbuf=NULL;
	return(n);
}
#ifdef KR_headers
c_sfe(a) cilist *a; /* check */
#else
int c_sfe(cilist *a) /* check */
#endif
{	unit *p;
	f__curunit = p = &f__units[a->ciunit];
	if(a->ciunit >= MXUNIT || a->ciunit<0)
		err(a->cierr,101,"startio");
	if(p->ufd==NULL && fk_open(SEQ,FMT,a->ciunit)) err(a->cierr,114,"sfe")
	if(!p->ufmt) err(a->cierr,102,"sfe")
	return(0);
}
integer e_wsfe(Void)
{
	int n = en_fio();
	f__fmtbuf = NULL;
	return n;
}
