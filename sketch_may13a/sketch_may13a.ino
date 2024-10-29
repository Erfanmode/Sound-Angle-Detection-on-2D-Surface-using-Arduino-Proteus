#define pin_mic1 A0
#define pin_mic2 A1
#define N1 50
#define peaks_n 8
#define pi 3.14159
int saved_cosine[257];
int theta=0;

bool signal=true;
int V1[N1];
int V2[N1];
int t1,t2,dt;
void setup() {
  Serial.begin(9600);
  analogReference(EXTERNAL);
  pinMode(pin_mic1,INPUT);
  pinMode(pin_mic2,INPUT);
  save_cosine(saved_cosine);
}



void loop() { 
 int read_value=256; 
 while(abs(read_value-256)<130){read_value=analogRead(pin_mic1);  //to detect high amplitude sounds and not the noise
                            }
    
   for (int i=0;i<N1;i++){
          t1=micros();
        
         V1[i]=analogRead(pin_mic1);
         V2[i]=analogRead(pin_mic2);
         int F=500;
        // V1[i]=30*cos(2*3.14*F*(i*0.000224)+0.79)+40*cos(2*3.14*1.5*F*(i*0.000224)+4+0.79)+3; 
         //V2[i]=23*cos(2*3.14*F*(i*0.000224))+43*cos(2*3.14*1.5*F*(i*0.000224)+4)+20;
        //Serial.print(V1[i]); Serial.print(" ");
        //Serial.println(V2[i]);
        
        
       t2=micros();
        dt=t2-t1;
        //Serial.println(dt);
   } 
       int phase1[N1]={0};int phase2[N1]={0};
       int peaks1[3*peaks_n]={0};int peaks2[3*peaks_n]={0};
       int frequency=4460;
          nha(V1,N1,frequency,phase1,peaks1);    
          nha(V2,N1,frequency,phase2,peaks2);
            float diff;
           int n_start=2*N1/40;
           float delay;
 /*         for(int i=n_start;i<N1/5;i++){  diff=phase1[i]-phase2[i]; if(diff<0){diff=diff+360;}
                                  // Serial.print(diff); 
                                  Serial.print(phase1[i]);Serial.print(" ");Serial.print(phase2[i]);
                                   //////float F=frequency*(i+1)/N1;
                                   //////if(diff>0){delay=diff*1.0/360/F;}
                                  Serial.print(" ");Serial.println((long)frequency/10*i/N1 );
                                  }    */
    
        
 /*  for(int j=0;j<peaks_n;j++){
        Serial.print(peaks1[3*j]);Serial.print(" ");
        Serial.print(peaks2[3*j]);Serial.print(" ");
        Serial.print(peaks1[3*j+1]);Serial.print(" ");
        Serial.println(peaks2[3*j+1]);//Serial.print(" ");
       // Serial.print(peaks1[3*j+2]);Serial.print(" ");
       // Serial.println(peaks2[3*j+2]); 
       } */
       
      theta=find_theta(peaks1,peaks2);
      Serial.println(theta);
      
        
}



void nha (int *x,int N,float Frequency,int *out_im,int *peaks)
{
  int out_r[N]={0};
  long th;//To scale theta from 2*pi to 1024
  for(int k=1;k<=N/5;k++){
    for(int j=0;j<N;j++){
      th=(long)j*k*1024/N;
     out_r[k]+=x[j]*cos(th*0.0061359);
     out_im[k]-=x[j]*sin(th*0.0061359);
     //out_r[k]+=fast_cosine(*(x+j) , th);
     // out_im[k]-=fast_sine(*(x+j) , th);
     }
     //float s; s=sqrt(pow(atan(out_im[k]*1.0/out_r[k]),2)+1) * abs(out_r[k]);out_im[k]=s;                   
     /* if(s>=120){
        if(out_r[k]>0)  out_im[k]=atan(out_im[k]*1.0/out_r[k])*57.3;               
        if(out_r[k]<=0) out_im[k]=atan(out_im[k]*1.0/out_r[k])*57.3 + 180;
       }
      else{out_im[k]=0;} */
           
  } 
 
   int temp1;
   float size,temp2;
     for(int i=0;i<N/5;i++){
          bool is_peak=false;
          size=sqrt( pow(atan(out_im[i]*1.0/out_r[i]),2)+1 )* abs(out_r[i]);
          if(size>0){
             int j=0;
             temp1=0;temp2=0;
              while(j<peaks_n){
                        if( size-*(peaks+3*j+2) > temp2+0.0001 ){
                            temp2=size-*(peaks+3*j+2);
                            temp1=j; 
                            is_peak=true;      
                        }
                j++;
              } 
              j=temp1;
                  if(is_peak){ *(peaks+3*j)=(long)Frequency*i/N;  //initial Frequency 
                              if(out_r[i]>0)  *(peaks+3*j+1)=atan(out_im[i]*1.0/out_r[i])*57.3;               
                              if(out_r[i]<=0) *(peaks+3*j+1)=atan(out_im[i]*1.0/out_r[i])*57.3 + 180;
                                                              //initial Phase
                              *(peaks+3*j+2)=size ;  //initial Amplitude with scale  
                                //Serial.print(*(peaks+3*j+2));Serial.print(" ");Serial.println(*(peaks+3*j));
                             }
              
          }
       }
     
  dsm(N,peaks,x,Frequency); 
  sort_f(peaks);// to sort peaks arrays based on their frequency value    
}



void dsm(int N,int *peaks, int *x,int Freq){
  int n1=30;
  float F,dF1,dF2,dF3,m;
 for(int j=0;j<peaks_n;j++){
        float Mu=1;  
    for(int k1=0;k1<n1;k1++){
           F=0;dF1=0;dF3=0;dF2=0;          
           Mu=Mu*0.99;
            for(int i=0;i<N;i++){
                m=*(x+i)-*(peaks+j*3+2)*cos(  2*pi*(*(peaks+3*j))*i /Freq + *(peaks+3*j+1)*2*pi /360 );
                F+=pow(m,2);
               dF1+= m * *(peaks+j*3+2) *sin(  2*pi*(*(peaks+3*j))*i /Freq + *(peaks+3*j+1)*2*pi /360 ) *4*pi*i/Freq;
               dF2+= m * *(peaks+j*3+2) *sin(  2*pi*(*(peaks+3*j))*i /Freq + *(peaks+3*j+1)*2*pi /360 )*2*pi/360;
               dF3+= m * (-2)*cos(  2*pi*(*(peaks+3*j))*i /Freq + *(peaks+3*j+1)*2*pi /360 );
            }
               
          
/*Serial.print(*(peaks+j*3)); Serial.print(" ");  Serial.print(*(peaks+j*3+1));  Serial.print(" ");
             Serial.print(*(peaks+j*3+2));  Serial.print(" ");  Serial.print(F/10000);
             Serial.print(" ");  Serial.println(k1); */
             *(peaks+j*3)-=Mu*dF1 /N  ;  
             *(peaks+j*3+1)-=Mu*dF2 /N  ;
             *(peaks+j*3+2)-=Mu*dF3 /N  ;  
             while(*(peaks+j*3+1)>360){*(peaks+j*3+1)-=360;} while(*(peaks+j*3+1)<=0){*(peaks+j*3+1)+=360;}    
              
    }
    if(*(peaks+j*3+2)<0){*(peaks+j*3+2)= (-*(peaks+j*3+2)); *(peaks+j*3+1)-=180; } 
    *(peaks+3*j)=abs(*(peaks+3*j));     
      
  } 
}


void sort_f(int *p){
  int temp=0;
  for(int i=0;i<peaks_n;i++){
    for(int j=0;j<peaks_n-i;j++){
      if(*(p+3*j+3)<*(p+3*j)){
        temp=*(p+3*j+3);
        *(p+3*j+3)=*(p+3*j);
        *(p+3*j)=temp;

        temp=*(p+3*j+4);
        *(p+3*j+4)=*(p+3*j+1);
        *(p+3*j+1)=temp;

        temp=*(p+3*j+5);
        *(p+3*j+5)=*(p+3*j+2);
        *(p+3*j+2)=temp;
      }
    }
  }
}

float find_theta(int *peaks1,int *peaks2){
  float sum=0,counter=0;
  int diff=0;
  for(int i=0;i<peaks_n;i++){
    for(int j=0;j<peaks_n;j++){
      if((peaks1[3*i]-peaks2[3*j])<50){
        diff=peaks1[3*i+1]-peaks2[3*j+1];
        while(diff<-180){diff+=360;}while(diff>180){diff-=360;}
        sum+= ( ( (float)(diff)/(peaks1[3*i]+peaks2[3*j])*2 / 360 )-0.000224)/0.000382 ;
        counter++;
      }
    }
  }
return((int)(sum/counter * 100));

}



void save_cosine(int *p)
{
 for(int i=0;i<256;i++)
 {
   *(p+i)=1000*cos((float)3.1415/2.0*i/256);
 }
  *(p+256)=0;
}

float fast_cosine(int amp , int th0)
{ // if th is 2*pi the input  for th is 1024
  while (th0>1024){
    th0=th0-1024;
  }
  int th=(int)th0;
  if(th>768){return ((float)amp*saved_cosine[1024-th] / 1000); }
  else if(th>512){return ((float)(-amp)*saved_cosine[th-512] / 1000); }
  else if(th>256){return ((float)(-amp)*saved_cosine[512-th] / 1000);}
  else{return( (float)amp*saved_cosine[th] / 1000 );}

}


float fast_sine(int amp , int th)
{
  return (-fast_cosine(amp,256+th));
}


