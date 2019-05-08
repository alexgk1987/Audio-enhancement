#include <cmath>
#include <fftw3.h>

#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <Windows.h>
#include <conio.h>
#include <stdio.h>
#include <fftw3.h>
#include <string.h>
#define M_PIl 3.141592653589793238462643383279502884L
#define REAL 0
#define IMAG 1
#define AMPLV 0
#define ANGLEV 1
#define array_size 16384/2
//#define array_size 1024
#define stride (16384/2)-500 // 703
#define steps (array_size - stride)
#define amplification  (3000.0)/(8.0/80)
#define scale 300000.0

//hp 450 800
using namespace std;

#include <alphabeta.h>
#include <alphabeta.cpp>

typedef struct  WAV_HEADER
{
    char                RIFF[4];        // RIFF Header      Magic header
    unsigned long       ChunkSize;      // RIFF Chunk Size
    char                WAVE[4];        // WAVE Header
    char                fmt[4];         // FMT header
    unsigned long       Subchunk1Size;  // Size of the fmt chunk
    unsigned short      AudioFormat;    // Audio format 1=PCM,6=mulaw,7=alaw, 257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM
    unsigned short      NumOfChan;      // Number of channels 1=Mono 2=Sterio
    unsigned long       SamplesPerSec;  // Sampling Frequency in Hz
    unsigned long       bytesPerSec;    // bytes per second
    unsigned short      blockAlign;     // 2=16-bit mono, 4=16-bit stereo
    unsigned short      bitsPerSample;  // Number of bits per sample
    char                Subchunk2ID[4]; // "data"  string
    unsigned long       Subchunk2Size;  // Sampled data length

} wav_hdr;

struct au_bytes
{
    unsigned char left_l;
    unsigned char left_h;
    unsigned char right_l;
    unsigned char right_h;
} ;

struct au_int
{
    short int left;
    short int right;
};


long double *xk_1=(long double*)  calloc((array_size/2+1),sizeof(long double));
long double *vk_1=(long double*)  calloc((array_size/2+1),sizeof(long double));
long double *xk=(long double*)  calloc((array_size/2+1),sizeof(long double));
long double *vk=(long double*)  calloc((array_size/2+1),sizeof(long double));
long double *arin=(long double*)  calloc((array_size/2+1),sizeof(long double));
long double *arout=(long double*)  calloc((array_size/2+1),sizeof(long double));



long double *xk_1_re=(long double*)  calloc((array_size/2+1),sizeof(long double));
long double *vk_1_re=(long double*)  calloc((array_size/2+1),sizeof(long double));
long double *xk_re=(long double*)  calloc((array_size/2+1),sizeof(long double));
long double *vk_re=(long double*)  calloc((array_size/2+1),sizeof(long double));

long double *xk_1_im=(long double*)  calloc((array_size/2+1),sizeof(long double));
long double *vk_1_im=(long double*)  calloc((array_size/2+1),sizeof(long double));
long double *xk_im=(long double*)  calloc((array_size/2+1),sizeof(long double));
long double *vk_im=(long double*)  calloc((array_size/2+1),sizeof(long double));


long double dt = 0.5, rk=0;


// envelope threshold
#define offset 0.1
// lowpass envelope
//0.05 0.04 0.08 -- 0.05 +1.25
//0.0005
long double  a = 0.1;
long double b=0.8*((2.0-a*a-2.0*sqrt(1.0-a*a) )/(a*a)  )*1.25;
//  f
long double  are = 0.5;//0.1/2.0;
long double bre=2.0*(0.8*((2.0-are*are-2.0*sqrt(1.0-are*are) )/(are*are)  )); //4.0
// lowpass f
long double  aim = 0.025;//0.1/2.0;
long double bim=0.8*((2.0-aim*aim-2.0*sqrt(1.0-aim*aim) )/(aim*aim)  ); //4.0




inline void alphabeta(long double *arin, long double *arout, int infos,long double *a_ar, long double *b_ar )
{
    int i;
    for (i=1; i < infos; i++)
    {

        xk[i] = xk_1[i] + ( vk_1[i] * dt );
        vk[i] = vk_1[i];
        rk = arin[i] - xk[i];
        xk[i] += a_ar[i] * rk;
        vk[i] += ( b_ar[i] * rk ) / dt;

        xk_1[i] = xk[i];
        vk_1[i] = vk[i];

        arout[i]=xk[i];
    }
}


inline void alphabeta_re(long double *arin, long double *arout, int infos)
{
    int i;
    for (i=0; i < infos; i++)
    {

        xk_re[i] = xk_1_re[i] + ( vk_1_re[i] * dt );
        vk_re[i] = vk_1_re[i];
        rk = arin[i] - xk_re[i];
        xk_re[i] += are * rk;
        vk_re[i] += ( bre * rk ) / dt;

        xk_1_re[i] = xk_re[i];
        vk_1_re[i] = vk_re[i];

        arout[i]=xk_re[i];

    }
}


inline void alphabeta_im(long double *arin, long double *arout, int infos)
{
    int i;
    for (i=0; i < infos; i++)
    {

        xk_im[i] = xk_1_im[i] + ( vk_1_im[i] * dt );
        vk_im[i] = vk_1_im[i];

        rk = arin[i] - xk_im[i];


        xk_im[i] += aim * rk;
        vk_im[i] += ( bim * rk ) / dt;

        xk_1_im[i] = xk_im[i];
        vk_1_im[i] = vk_im[i];

        arout[i]=xk_im[i];

    }
}


void Boxcar_setup (int elements, long double *welc)
{
    int count1;
    for (count1=0; count1 < elements; count1++)
    {
        welc[count1]= (long double)1.0;
    }
}


long double tukey_sub (int i,int N, long double alpha)
{
    long double anm12 = (long double)0.5*alpha*((long double)N-1);

    if( (long double)i <= anm12 )
    {
        return (long double)0.5*((long double)1.0+cos(M_PIl*((long double)i/anm12 - (long double)1.0)));
    }
    else if ( i < (N-1)*(1.0-0.5*alpha) )
    {
        return 1;
    }
    else
    {
        return (long double)0.5*((long double)1.0+cos(M_PIl*((long double)i/anm12 - (long double)2.0/alpha + (long double)1.0)));
    }
}



//HANN!
void Hamming_setup (int elements, long double *hamming)
{
    int count1;
    long double hamming_f, high, low;
    low=(long double)elements-(long double)1.0;
    for (count1=0; count1 < elements; count1++)
    {
        high=(long double)2.0*(long double)M_PIl*(long double)count1;
        hamming_f=(long double)0.5*((long double)1.0-(long double)cos(high/low));
        hamming[count1]= hamming_f;
    }
//cout << "\n" << (float)hamming[1] <<"\t" << (float) hamming[elements-2];

}


void Nuttall_setup (int elements, long double *nuttall)
{
    int count1;
    long double a0=0.3635819;
    long double a1=0.4891775;
    long double a2=0.1365995;
    long double a3=0.0106411;

    long double ax, bx, cx;
    for (count1=0; count1 < elements; count1++)
    {
        ax=a1*(long double)cos((long double)2.0*M_PIl*(long double) count1/(long double) elements-(long double)1.0);
        bx=a2*(long double)cos((long double)4.0*M_PIl*(long double) count1/(long double) elements-(long double)1.0);
        cx=a3*(long double)cos((long double)6.0*M_PIl*(long double)count1/(long double)  elements-(long double)1.0);
        nuttall[count1]=a0-ax+bx-cx;
    }
}

void fir_spectre(int ar_size_i, int fil_size_i, long double *in_array, long double *coeff, long double min_val)
{
    long double accumulator=0;
    long double int_array[ar_size_i+200];
    int a,b;
//pad array with f value
    for(a=0; a< fil_size_i/2+1; a++)
        int_array[a] = min_val;
//write other walues
    for(a=fil_size_i/2+1; a< ar_size_i+fil_size_i/2; a++)
        int_array[a] = in_array[a-fil_size_i/2];

    for(a=0; a< fil_size_i/2; a++)
        int_array[ar_size_i+fil_size_i/2+a] = min_val;

    for(a=1; a< ar_size_i; a++)
    {
        accumulator=0.0;
        for(b=0; b < fil_size_i; b++)
            accumulator+=int_array[a+b]*coeff[b];
        in_array[a]=accumulator;
    }

}
inline long double modpi(long double angle )
{
    long double twoPi = 2.0 * 3.141592865358979;
    return angle - twoPi * floor( angle / twoPi );
    //return fmod(angle, twoPi);
}


// 100 400 90 LOWPASS!

void process_vinyl(long int size1, long double *windowf, short int *datain, long double *dataout, long double compress)
{
bool ftime=true;

    int i=0;
    long double *inar=(long double*)  fftwl_alloc_real(array_size);

    long double *bcoeff=(long double*)  calloc(array_size,sizeof(long double));

    fftwl_complex *outar=(fftwl_complex*)  fftwl_alloc_complex(array_size);

    long int position=0;
    long int position2=0;
    long double array_size_d=array_size;
    int cnt1=0;

    long double *ffin=(long double*)   calloc(array_size,sizeof(long double));

    long double *ffout=(long double*)   calloc(array_size,sizeof(long double));

    long double *filtout=(long double*)   calloc(array_size,sizeof(long double));
long double *a_ar=(long double*)   calloc(array_size,sizeof(long double));
long double *b_ar=(long double*)   calloc(array_size,sizeof(long double));

    long double *ampl=(long double*)   calloc(array_size,sizeof(long double));

      long double *array_x=(long double*)   calloc(array_size,sizeof(long double));
        long double *array_y=(long double*)   calloc(array_size,sizeof(long double));

    long double *ampout=(long double*)   calloc(array_size,sizeof(long double));
long double *abin=(long double*)   calloc(array_size,sizeof(long double));
long double *about=(long double*)   calloc(array_size,sizeof(long double));
    long double *intout=(long double*)   calloc(array_size,sizeof(long double));
    fftwl_plan forw = fftwl_plan_dft_r2c_1d(array_size, inar, outar,FFTW_PATIENT);
    fftwl_plan rev = fftwl_plan_dft_c2r_1d(array_size, outar, inar, FFTW_PATIENT);


    int inita, M=30, N=array_size/2+1;

    long double **phases = calloc(N, sizeof (long double*));
    if (phases)
    {
        for (i = 0; i < N; i++)
        {
            phases[i] = calloc(30, sizeof *phases[i]);
        }
    }



    long double **filter_a = calloc(N,sizeof (long double*));
    if (filter_a)
    {
        for (i = 0; i < N; i++)
        {
            filter_a[i] = calloc(M,sizeof *filter_a[i]);
        }
    }

     long double **first_d = calloc(N,sizeof (long double*));
    if (first_d)
    {
        for (i = 0; i < N; i++)
        {
            first_d[i] = calloc(4,sizeof *filter_a[i]);
        }
    }
    long double **second_d = calloc(N,sizeof (long double*));
    if (second_d)
    {
        for (i = 0; i < N; i++)
        {
            second_d[i] = calloc(4,sizeof *filter_a[i]);
        }
    }

        long double **third_d = calloc(N,sizeof (long double*));
    if (third_d)
    {
        for (i = 0; i < N; i++)
        {
            third_d[i] = calloc(4,sizeof *filter_a[i]);
        }
    }

int cnt2;
for(cnt1=0; cnt1<1000; cnt1++) {
for(cnt2=0; cnt2<=23; cnt2++) {
    phases[cnt1][cnt2]=0;
}}

    for(cnt1=0; cnt1 < array_size/2+1; cnt1++ )
        xk_1[cnt1]=vk_1[cnt1]=xk[cnt1]=vk[cnt1]=0;
    for(cnt1=0; cnt1 < array_size/2+1; cnt1++ )
        xk_1_re[cnt1]=vk_1_re[cnt1]=xk_re[cnt1]=vk_re[cnt1]=0;
    for(cnt1=0; cnt1 < array_size/2+1; cnt1++ )
        xk_1_im[cnt1]=vk_1_im[cnt1]=xk_im[cnt1]=vk_im[cnt1]=0;
    long double magnetude[array_size/2+1][2];

    long double temp1=0;
// BUTTERWORTH FILTERING
    /*
    long double butter, bfreq=48000, countd, temp1;
    long double bstop=3500;


    for(cnt1=0; cnt1 < array_size/2+1; cnt1++ )
    {
    countd=(long double) cnt1;
    bfreq=(long double) countd/(((long double)array_size_d/(long double)2.0)+(long double)1.0);
    butter=(long double)1.0/sqrt((long double)1.0+(long double)pow((bfreq/bstop),8.0));
    bcoeff[cnt1]=(long double)1.0-butter;

    }
    */
    int counter=0;
    long double tempc;
    for(cnt1=0 ; cnt1 < array_size/2+1; cnt1++)
    {
    a_ar[cnt1]=tempc=1e-4;//*((long double)cnt1);
    b_ar[cnt1]=0.8*((2.0-tempc*tempc-2.0*sqrt(1.0-tempc*tempc) )/(tempc*tempc)  );
    }
//-------------------------------------------------------------------------------------------------
    while (position <(long)size1/4 )
    {
        counter++;
        for(cnt1=0; cnt1 < array_size; cnt1++)
        {
            temp1=(long double)datain[position]/((long double)scale);
            inar[cnt1]=temp1*windowf[cnt1];
            if(inar[cnt1]==0) inar[cnt1]=1e-20;
            position++;
        }

        position-=stride;



        for(cnt1=0; cnt1 < array_size; cnt1++)
            inar[cnt1]/=array_size_d;

        fftwl_execute(forw);

        for(cnt1=0; cnt1 < array_size/2+1; cnt1++ )
        {
            magnetude[cnt1][AMPLV]=
                (long double)sqrt((long double)outar[cnt1][REAL]* (long double)outar[cnt1][REAL]+
                                  (long double)outar[cnt1][IMAG]*(long double)outar[cnt1][IMAG]);

            magnetude[cnt1][ANGLEV]=(long double)
                                    atan2l((long double)outar[cnt1][IMAG],(long double)outar[cnt1][REAL]);
        }

        int inta;
        long double newphase, coeff1, nrml;
// phases information
long double array__x, array__y, product, dFiK;
        for(cnt1=0 ; cnt1 < 1000; cnt1++)
        {
/*
            phases[cnt1][22]=magnetude[cnt1][ANGLEV]*
            (long double)0.999-phases[cnt1][0]*(long double)0.999;
            phases[cnt1][22]=modpi(phases[cnt1][22]);
            phases[cnt1][21]=phases[cnt1][22]+phases[cnt1][0];
            phases[cnt1][0]=phases[cnt1][21];
*/

            phases[cnt1][23]=magnetude[cnt1][ANGLEV]-phases[cnt1][22];
            if(abs(phases[cnt1][23])> M_PIl)
            {
                if (phases[cnt1][23]> 0)
                    newphase=phases[cnt1][23]-(long double)2.0*M_PIl;
                else
                    newphase=phases[cnt1][23]+(long double)2.0*M_PIl;
            }
            else newphase=phases[cnt1][23];
            phases[cnt1][22]=magnetude[cnt1][ANGLEV];

            if (ftime==false) phases[cnt1][6]+=newphase;
            else phases[cnt1][6]= magnetude[cnt1][ANGLEV];

    ffin[cnt1]=phases[cnt1][6];

    //phases[cnt1][1]=phases[cnt1][0]*0.25+phases[cnt1][1]*(1.0-0.25);
int c,d;
double swap;

//for (c = 1 ; c < 6 - 1; c++)
  //{
    //for (d = 1 ; d < 6 - c - 1; d++)
    //{
      //if (phasesd[cnt1][d] > phasesd[cnt1][d+1]) /* For decreasing order use < */
      //{
       // swap       = phasesd[cnt1][d];
        // phasesd[cnt1][d]   = phasesd[cnt1][d+1];
      //  phasesd[cnt1][d+1] = swap;
    //  }
  //  }
//  }


//            phases[cnt1][2]=newphase;

}

alphabeta_re(ffin, ffout, 1000);

/*
for(cnt1=0 ; cnt1 < array_size/2+1; cnt1++)
{
for(inita=0; inita < 3; inita++) filter_a[cnt1][inita]=filter_a[cnt1][inita+1];
filter_a[cnt1][3]=pow(magnetude[cnt1][AMPLV]+1e-11, 0.5);
}
*/

//for(cnt1=0 ; cnt1 < array_size/2+1; cnt1++) abin[cnt1]=magnetude[cnt1][AMPLV];
//alphabeta(abin, about,array_size/2+1, a_ar, b_ar );



for(cnt1=0 ; cnt1 < array_size/2+1; cnt1++)
{
filter_a[cnt1][3]=filter_a[cnt1][2];
filter_a[cnt1][2]=filter_a[cnt1][1];
filter_a[cnt1][1]=filter_a[cnt1][0];
filter_a[cnt1][0]=magnetude[cnt1][AMPLV];//-about[cnt1];
}

for(cnt1=0 ; cnt1 < array_size/2+1; cnt1++)
{
first_d[cnt1][1]=first_d[cnt1][0];
first_d[cnt1][0]=filter_a[cnt1][0]-filter_a[cnt1][1];
}

for(cnt1=0 ; cnt1 < array_size/2+1; cnt1++)
{
second_d[cnt1][1]=second_d[cnt1][0];
second_d[cnt1][0]=first_d[cnt1][0]-first_d[cnt1][1];

}
for(cnt1=0 ; cnt1 < array_size/2+1; cnt1++)
{
third_d[cnt1][0]=second_d[cnt1][0]-second_d[cnt1][1];
}


long double quickacc, deriv, koef, low,dderiv, masking, maskingLIN;
for(cnt1=1 ; cnt1 < array_size/2; cnt1++)
{
   // quickacc=(abs(third_d[cnt1][0]))/ampout[cnt1];
  deriv= abs(third_d[cnt1][0]);
 dderiv=deriv/(long double) cnt1;
 low=ampout[cnt1]/(dderiv+1e-60);
// masking=dderiv/pow(ampout[cnt1], 1/1.0);
 maskingLIN=dderiv/(ampout[cnt1]+1e-60);
 //if (cnt1==200) cout << "\n" << (float) koef;
if(                        dderiv >= 1e-4   || low < 1e-2 )  { ampout[cnt1]=0; about[cnt1]=0;}
if( low < 5e-3  )             magnetude[cnt1][AMPLV]=0;

}
long double normalize;
for(cnt1=0 ; cnt1 < array_size/2; cnt1++)
{
if (cnt1 >50 )normalize=((long double) cnt1/((long double) array_size/2.0))/2.0;
else normalize=0;
ampout[cnt1]=exp(log(magnetude[cnt1+1][AMPLV]+10e-6)+log(magnetude[cnt1][AMPLV]+10e-6)*(normalize)+normalize*2.0);//+4e-6;
}
//cout << "\n" << (float) abs(third_d[200][0])* (float) (magnetude[200][AMPLV]);

long double stemp, fst, scnd, ka,kb, t1, t2,ts1,ts2;
int av, bv, sum, dif, iter;
//ka=9e-4;
#define thres 5e-2
ftime=false;
if(false)
{
for (cnt1; cnt1 < array_size/2; cnt1++)
{
ampout[cnt1]+=thres;
about[cnt1]+=thres;
}
}
ka=0.25;//1e-2;
//cout << "\nDone " << counter*240;
if(false) {
//sum
for(av=1 ; av < array_size/2; av++)
{
// correction for summary and difference IM components
    for(bv=av+1 ; bv < array_size/2; bv++)
   {
    sum =bv+av;
    dif=bv-av;

    if(sum >  array_size/2   ) goto dif;
    else
    {
        t1=about[av];
        t2=about[bv];

    ampout[sum]-=(t1*t2)*ka;
     if(ampout[sum] < 0  ) ampout[sum]=0;
    }


    dif:
    if( dif <0 /*||dif==av || dif==bv*/ ) goto fin;
else{
        t1=about[av];
        t2=about[bv];
        ampout[dif]-=(t1*t2)*ka;
    if(ampout[dif] < 0) ampout[dif]=0;
    }

   fin:
    sum=sum;
   }
}
}

//for(cnt1=0 ; cnt1 <array_size/2; cnt1++) ampout[cnt1]+=magnetude[cnt1][AMPLV];

for(cnt1=0 ; cnt1 <array_size/2+1; cnt1++) filtout[cnt1]=magnetude[cnt1][ANGLEV];
for(cnt1=0 ; cnt1 <1000; cnt1++) filtout[cnt1]=ffout[cnt1];
//cout << "\n" << (float) filtout[500];

        for(cnt1=0; cnt1 < array_size/2+1; cnt1++ )
        {
            outar[cnt1][REAL]=(long double)ampout[cnt1]*(long double)cos((long double)  filtout[cnt1]);
            outar[cnt1][IMAG]=(long double)ampout[cnt1]*(long double)sin((long double)  filtout[cnt1]);
        }

        for(cnt1=array_size/2+1; cnt1 < array_size; cnt1++ )
        {
            outar[cnt1][REAL]=0.0;
            outar[cnt1][IMAG]=0.0;
        }

        fftwl_execute(rev);

        for(cnt1=0; cnt1 < array_size; cnt1++)
        {
            dataout[position2]+=(long double)inar[cnt1];
            position2++;
        }


        position2-=stride;
    }


    fftwl_free(outar);
    fftwl_free(inar);

    free(filter_a);

    free(ampout);
    free(filtout);
    free(phases);
    free(bcoeff);

    fftwl_destroy_plan(forw);
    fftwl_destroy_plan(rev);





}

//----------------MAIN END-----------------------------------------------------------

void record_to_file(long int size1, long double *out_le, long double *out_ri, FILE *wavFileO)
{
    long int position=0;
    int cnt1;
    au_bytes au_ff[array_size];
    au_int au_array[array_size];
    bool zeroes=true;
    while (position < size1/4 )
    {
        if (position > 44100*5)
            zeroes=false;
        //load block to structure
        for (cnt1=0; cnt1 < array_size; cnt1++)
        {
            //10000
            au_array[cnt1].left= (long double)out_le[position]*(long double)amplification;
            au_array[cnt1].right= (long double)out_ri[position]*(long double)amplification;
            if(zeroes)
            {
                au_array[cnt1].left= 0;
                au_array[cnt1].right= 0;
            }
            position++;
        }
        // write structure to file as a form of bloc
        for (cnt1=0; cnt1 < array_size; cnt1++)
        {
            au_ff[cnt1].left_h= (short int) au_array[cnt1].left >>8;
            au_ff[cnt1].left_l=(short int) au_array[cnt1].left & 255;

            au_ff[cnt1].right_h= (short int) au_array[cnt1].right >> 8;
            au_ff[cnt1].right_l= (short int) au_array[cnt1].right & 255;
        }
        fwrite(&au_ff, 1,  array_size*4, wavFileO);
    }
}
//----------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    char fname1[50];
    long double compress;
    //define structures to hold the data
    au_bytes au_ff[array_size];
    au_int au_array[array_size];

    int cnt1;

    wav_hdr wavHeader;
    FILE * wavFile;
    FILE * wavFileO;
    FILE * LogFile;
    size_t read_bytes;

    int headerSize= sizeof(wav_hdr);
    int filelength=0;
    long int size1, position;
    // array for f window
    long double windowf[array_size];

    if(argc == 1)
    {
        cout << "\nFile name-> ";
        cin >> fname1;
    }
    else
        strcpy(fname1, argv[1]);
    cout << "Step 1\nMake window finction";
Nuttall_setup(array_size, windowf);
//Boxcar_setup(array_size, windowf);
//Hamming_setup(array_size, windowf);

    wavFile=fopen(fname1,"rb");


    if(argc == 1)
    {
        wavFileO=fopen("out.wav","wb");
    }
    else if(argc==3)
        wavFileO=fopen(argv[2],"wb");
    else
        return 0;
    //LogFile=fopen("d:\\spectrum.log", "w");
    // read header
    fread(&wavHeader, headerSize,  1, wavFile);
    //write header to output file
    fwrite(&wavHeader, headerSize,  1, wavFileO);

    //first pass
    read_bytes=array_size*4;
    do
    {
        read_bytes=fread(&au_ff, 1, array_size*4, wavFile);
    }
    while(read_bytes==array_size*4);

    size1=ftell(wavFile);
    size1=size1-headerSize;
    //normalize array

    rewind(wavFile);
    //fclose(LogFile);

    fread(&wavHeader, headerSize,  1, wavFile);
    // allocate arrays

    short int *in_le=(short int*) calloc (size1/3,sizeof(short int));
    short int *in_ri=(short int*) calloc (size1/3,sizeof(short int));

    long double *out_le=(long double*) calloc (size1/2,sizeof(long double ) );
    long double *out_ri=(long double*) calloc (size1/2,sizeof(long double ) );

    if (in_le==NULL || in_ri==NULL || out_le==NULL ||out_ri==NULL)
        return 205;

//calculate(ar_values_i, ar_values_j);
    cout << "\nStep 2";
    long int pointer=0;

    do
    {
        //read bytes
        read_bytes=fread(&au_ff, 1, array_size*4, wavFile);
        //decode to integers
        for (cnt1=0; cnt1 < array_size; cnt1++)
        {
            au_array[cnt1].left  =   ((short int)  au_ff[cnt1].left_h) << 8;
            au_array[cnt1].left = au_array[cnt1].left | (short int)  au_ff[cnt1].left_l ;
            au_array[cnt1].right =   ((short int)  au_ff[cnt1].right_h) << 8;
            au_array[cnt1].right =au_array[cnt1].right |(short int)  au_ff[cnt1].right_l ;
        }
        //write to arrays
        for(cnt1=0; cnt1 < array_size; cnt1++)
        {
            in_le[pointer]=(short int)au_array[cnt1].left;
            in_ri[pointer]=(short int)au_array[cnt1].right;
            pointer++;
        }
    }
    while(read_bytes==array_size*4);


    cout << "\nProcessing left channel";
    process_vinyl(size1, windowf, in_le, out_le, compress);
    free(in_le);
    cout << "\nProcessing right channel";
    process_vinyl(size1, windowf, in_ri, out_ri, compress);

    record_to_file(size1, out_le, out_ri, wavFileO);
    fclose(wavFileO);
    fclose(LogFile);
    free(out_le);
    free(out_ri);

    free(in_ri);
    fftwl_cleanup();
    cout << "\n 16 bit int is ->> bytes" << sizeof(__int16);
    cout << "\n Int is ->> bytes" << sizeof(int);

    return 0;
}
