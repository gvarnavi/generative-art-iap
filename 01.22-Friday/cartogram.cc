#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#ifdef _OPENMP
#include "omp.h"
inline double wtime() {return omp_get_wtime();}
#else
#include <ctime>
inline double wtime() {return static_cast<double>(clock())*(1./CLOCKS_PER_SEC);}
#endif

// Open a file.
FILE* open_file(const char *filename,const char* mode) {
    FILE *fp=fopen(filename,mode);
    if(fp==NULL) {
        fprintf(stderr,"Can't open file '%s'\n",filename);
        exit(1);
    }
    return fp;
}

// Read contents from a file.
void read_file(void *ptr,size_t size,size_t count,FILE *fp) {
    if(fread(ptr,size,count,fp)!=count) {
        fprintf(stderr,"Can't read from file\n");
        exit(1);
    }
}

// Integrate feature and reference map fields forward by one timestep, using ADI method.
void step(int m,int n,double *u,double *x,double *y,double dt,double h) {
	int mn=m*n,i,j,i2,j2,k;
	double nu=dt/h/h,fac=nu/2.,fac2=fac/2.,a,w,
		  *b=y+mn,*c=b+mn,*d=c+mn,*d2=d+mn,
		  vx,vy;

	// Upwinded reference map update.
	// First treat x derivative implicitly:
	*b=1.;*c=0.;*d=*x;*d2=*y;
	for(k=1;k<mn;k++) {
		i=k%m;j=k/m;
		if((i>0)&&(i<m-1)) {
			vx=-(u[n*(i+1)+j]-u[n*(i-1)+j])*fac2/u[n*i+j];
			if(vx>0) {a=-vx;b[k]=1+vx;c[k]=0.;}
			else {a=0.;b[k]=1-vx;c[k]=vx;}
		}
		else {a=0.;b[k]=1.;c[k]=0.;}
		w=a/b[k-1];
		if((j>0)&&(j<n-1)) {
			vy=-(u[n*i+j+1]-u[n*i+j-1])*fac2/u[n*i+j];
			if(vy>0) {d[k]=x[n*i+j]*(1-vy)+vy*x[n*i+j-1];
					 d2[k]=y[n*i+j]*(1-vy)+vy*y[n*i+j-1];}
			else {d[k]=x[n*i+j]*(1+vy)-vy*x[n*i+j+1];
				 d2[k]=y[n*i+j]*(1+vy)-vy*y[n*i+j+1];}
		}
		else {d[k]=x[n*i+j];d2[k]=y[n*i+j];}
		b[k]-=w*c[k-1];
		d[k]-=w*d[k-1];
		d2[k]-=w*d2[k-1];
	}
	x[mn-1]=d[mn-1]/b[mn-1];
	y[mn-1]=d2[mn-1]/b[mn-1];
	for(k=mn-2;k>=0;k--) {
		i=k%m;j=k/m;
		i2=(k+1)%m;j2=(k+1)/m;
		x[n*i+j]=(d[k]-c[k]*x[n*i2+j2])/b[k];
		y[n*i+j]=(d2[k]-c[k]*y[n*i2+j2])/b[k];
	}
	// Next treat the y derivative implicitly:
	*b=1.;*c=0.;*d=*x;*d2=*y;
	for(k=1;k<mn;k++) {
		i=k/n;j=k%n;
		if((j>0)&&(j<n-1)) {
			vy=-(u[n*i+j+1]-u[n*i+j-1])*fac2/u[n*i+j];
			if(vy>0) {a=-vy;b[k]=1+vy;c[k]=0.;}
			else {a=0.;b[k]=1-vy;c[k]=vy;}
		}
		else {a=0.;b[k]=1.;c[k]=0.;}
		w=a/b[k-1];
		if((i>0)&&(i<m-1)) {
			vx=-(u[n*(i+1)+j]-u[n*(i-1)+j])*fac2/u[n*i+j];
			if(vx>0) {d[k]=x[n*i+j]*(1-vx)+vx*x[n*(i-1)+j];
					 d2[k]=y[n*i+j]*(1-vx)+vx*y[n*(i-1)+j];}
			else {d[k]=x[n*i+j]*(1+vx)-vx*x[n*(i+1)+j];
				 d2[k]=y[n*i+j]*(1+vx)-vx*y[n*(i+1)+j];}
		}
		else {d[k]=x[n*i+j];d2[k]=y[n*i+j];}
		b[k]-=w*c[k-1];
		d[k]-=w*d[k-1];
		d2[k]-=w*d2[k-1];
	}
	x[mn-1]=d[mn-1]/b[mn-1];
	y[mn-1]=d2[mn-1]/b[mn-1];
	for(k=mn-2;k>=0;k--) {
		i=k/n;j=k%n;
		i2=(k+1)/n;j2=(k+1)%n;
		x[n*i+j]=(d[k]-c[k]*x[n*i2+j2])/b[k];
		y[n*i+j]=(d2[k]-c[k]*y[n*i2+j2])/b[k];
	}

	// Finite difference feature map update.
	// First treat x derivative implicitly:
	*b=1+fac;*c=-fac;*d=*u*(1-fac)+fac*u[1];
	for(k=1;k<mn;k++) {
		i=k%m;j=k/m;
		if(i==0) {a=0.;b[k]=1+fac;c[k]=-fac;}
		else if(i==m-1) {a=-fac;b[k]=1+fac;c[k]=0;}
		else {a=-fac;b[k]=1+nu;c[k]=-fac;}
		w=a/b[k-1];
		if(j==0) d[k]=u[n*i+j]*(1-fac)+fac*u[n*i+j+1];
		else if(j==n-1) d[k]=u[n*i+j]*(1-fac)+fac*u[n*i+j-1];
		else d[k]=u[n*i+j]*(1-nu)+fac*(u[n*i+j-1]+u[n*i+j+1]);
		b[k]-=w*c[k-1];
		d[k]-=w*d[k-1];
	}
	u[mn-1]=d[mn-1]/b[mn-1];
	for(k=mn-2;k>=0;k--) {
		i=k%m;j=k/m;
		i2=(k+1)%m;j2=(k+1)/m;
		u[n*i+j]=(d[k]-c[k]*u[n*i2+j2])/b[k];
	}
	// Next treat the y derivative implicitly:
	*b=1+fac;*c=-fac;*d=*u*(1-fac)+fac*u[n];
	for(k=1;k<mn;k++) {
		i=k/n;j=k%n;
		if(j==0) {a=0.;b[k]=1+fac;c[k]=-fac;}
		else if(j==n-1) {a=-fac;b[k]=1+fac;c[k]=0;}
		else {a=-fac;b[k]=1+nu;c[k]=-fac;}
		w=a/b[k-1];
		if(i==0) d[k]=u[n*i+j]*(1-fac)+fac*u[n*(i+1)+j];
		else if(i==m-1) d[k]=u[n*i+j]*(1-fac)+fac*u[n*(i-1)+j];
		else d[k]=u[n*i+j]*(1-nu)+fac*(u[n*(i-1)+j]+u[n*(i+1)+j]);
		b[k]-=w*c[k-1];
		d[k]-=w*d[k-1];
	}
	u[mn-1]=d[mn-1]/b[mn-1];
	for(k=mn-2;k>=0;k--) {
		i=k/n;j=k%n;
		i2=(k+1)/n;j2=(k+1)%n;
		u[n*i+j]=(d[k]-c[k]*u[n*i2+j2])/b[k];
	}
}

int main(int argc,char **argv) {
	if(argc!=3) {
		fputs("Syntax: ./cartogram <map_size> <feature_name>\n",stderr);
		return 1;
	}

	// Allocate space for feature values and label map.
	int N=50,m,n,mn;
	if(strcmp(argv[1],"small")==0) {m=204;n=304;}
	else if(strcmp(argv[1],"med")==0) {m=408;n=608;}
	else if(strcmp(argv[1],"large")==0) {m=816;n=1216;}
	else {
		fprintf(stderr,"Invalid map file '%s'.\n",argv[1]);
        return 1;
	}
	mn=m*n;
	uint8_t *l=new uint8_t[mn];
	double *rh=new double[N];
	char buf[128];

	// Read in labels from binary file.
	sprintf(buf,"maps/%s.bin",argv[1]);
	FILE *fp=open_file(buf,"rb");
	read_file(l,sizeof(uint8_t),mn,fp);
	fclose(fp);

	// Read in feature values from binary file.
	sprintf(buf,"data/%s.bin",argv[2]);
	fp=open_file(buf,"rb");
	read_file(rh,sizeof(double),N,fp);
	fclose(fp);

	// Allocate space for features and tracer positions, and auxiliary memory.
	double *u=new double[7*mn],*x=u+mn,*y=x+mn;

	// Set the features on the map, skipping over background.
	double *up,r,rhavg=0.;
	int ct=0;
	uint8_t *lp;
	for(lp=l,up=u;lp<l+mn;lp++,up++) if(*lp>0) {
		r=rh[*lp-1];
		*up=r;rhavg+=r;ct+=1;
	}
	rhavg/=ct;
	// Compute scaling factor so that average feature value is 1.
	double sc=1./rhavg;

	// Set all background values equal to the average feature value.
	for(lp=l,up=u;lp<l+mn;lp++,up++) {
		if(*lp==0) *up=1.;
		else *up*=sc;
	}

	// Initialize reference map coordinates.
	double h=1.; // grid spacing.
	for(int i=0;i<m;i++) {
		for(int j=0;j<n;j++) {
			x[n*i+j]=h*i;
			y[n*i+j]=h*j;
		}
	}

	// Set timestep size and integration time.
	double dt=0.99*h*h,T=0.01*(m*m+n*n);
	int nsteps=int(ceil(T/dt));
	dt=T/nsteps;

	// Integrate, outputting frames.
	double wt=wtime();
	for(int k=0;k<nsteps;k++) {
		step(m,n,u,x,y,dt,h);
	}

	// Write resulting fields to file.
	sprintf(buf,"out/%s_fields.bin",argv[2]);
	fp=open_file(buf,"wb");
	fwrite(u,sizeof(double),3*mn,fp);
	fclose(fp);

	printf("Integrated to T=%.f using %d steps.\n",nsteps*dt,nsteps);
	printf("Elapsed time: %8.5f min.\n",(wtime()-wt)/60.);

	delete [] u;
	delete [] rh;
	delete [] l;
}
