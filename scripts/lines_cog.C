//

#include <cstdio>
#include <math.h>

struct LineData_t {
  double x0;
  double y0;
  double nx;
  double ny;
};

void lines_cog() {


  LineData_t data[] = {
    {0,0,1,0},
    {2,0,0,1},
    {0.,5.,1./sqrt(2.),-1./sqrt(2.)},
    {0,0,99,99}
  };
  
    // line parameterization: x0, y0, nx, ny

  int nlines = 0;

  for (int i=0; data[i].nx < 10; i++) nlines++;

  printf("nlines = %i\n",nlines);

  double sx(0), sy(0), snx2(0),snxny(0), sny2(0), snxnr(0), snynr(0);
  
  for (int i=0; i<nlines; i++) {
    LineData_t* l = &data[i];
    
    double nr = l->nx*l->x0+l->ny*l->y0;

    sx    += l->x0;
    sy    += l->y0;
    snx2  += l->nx*l->nx;
    snxny += l->nx*l->ny;
    sny2  += l->ny*l->ny;
    snxnr += l->nx*nr;
    snynr += l->ny*nr;
  }

  double x_mean, y_mean, nxny_mean, nx2_mean, ny2_mean, nxnr_mean, nynr_mean;

  x_mean    = sx/nlines;
  y_mean    = sy/nlines;
  nxny_mean = snxny/nlines;
  nx2_mean  = snx2/nlines;
  ny2_mean  = sny2/nlines;
  nxnr_mean = snxnr/nlines;
  nynr_mean = snynr/nlines;


  double d = (1-nx2_mean)*(1-ny2_mean)-nxny_mean*nxny_mean;

  double x0 = ((x_mean-nxnr_mean)*(1-ny2_mean)+(y_mean-nynr_mean)*nxny_mean)/d;
  double y0 = ((y_mean-nynr_mean)*(1-nx2_mean)+(x_mean-nxnr_mean)*nxny_mean)/d;

  printf("x_mean, y_mean = %10.4f %10.4f\n",x_mean, y_mean);
  printf("nx2_mean, nxy_mean, ny2_mean = %10.4f %10.4f %10.4f\n",nx2_mean, nxny_mean, ny2_mean);
  printf("d              = %10.4f\n",d);
  printf("x0, y0         = %10.4f %10.4f\n",x0,y0);
}
