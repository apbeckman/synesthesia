

// Shakespeare Quest by eiffie (solving for 2 joints)

#define size RENDERSIZE
#define SHADOWS
//perllin noise from madweedfall
#define firstOctave 3
#define octaves  8
#define persistence 0.6

//Not able to use bit operator like <<, so use alternative noise function from YoYo
//
//https://www.shadertoy.com/view/Mls3RS
//
//And it is a better realization I think
float noise(int x,int y)
{   
    float fx = float(x);
    float fy = float(y);
    
    return 2.0 * fract(sin(dot(vec2(fx, fy) ,vec2(12.9898,78.233))) * 43758.5453) - 1.0;
}

float smoothNoise(int x,int y)
{
    return noise(x,y)/4.0+(noise(x+1,y)+noise(x-1,y)+noise(x,y+1)+noise(x,y-1))/8.0+(noise(x+1,y+1)+noise(x+1,y-1)+noise(x-1,y+1)+noise(x-1,y-1))/16.0;
}

float COSInterpolation(float x,float y,float n)
{
    float r = n*3.1415926;
    float f = (1.0-cos(r))*0.5;
    return x*(1.0-f)+y*f;
    
}

float InterpolationNoise(float x, float y)
{
    int ix = int(x);
    int iy = int(y);
    float fracx = x-float(int(x));
    float fracy = y-float(int(y));
    
    float v1 = smoothNoise(ix,iy);
    float v2 = smoothNoise(ix+1,iy);
    float v3 = smoothNoise(ix,iy+1);
    float v4 = smoothNoise(ix+1,iy+1);
    
   	float i1 = COSInterpolation(v1,v2,fracx);
    float i2 = COSInterpolation(v3,v4,fracx);
    
    return COSInterpolation(i1,i2,fracy);
    
}

float PerlinNoise2D(float x,float y)
{
    float sum = 0.0;
    float frequency =0.0;
    float amplitude = 0.0;
    for(int i=firstOctave;i<octaves + firstOctave;i++)
    {
        frequency = pow(2.0,float(i));
        amplitude = pow(persistence,float(i));
        sum = sum + InterpolationNoise(x*frequency,y*frequency)*amplitude;
    }
    
    return sum;
}
float timeNoise = PerlinNoise2D(smoothTime, smoothTime);
float time;
float time2;
float time3;
/* //"solves" for two joints (there is a better way I'm sure)
void djsolve( vec3 a, vec3 b, vec3 l, vec3 rt, out vec3 j1, out vec3 j2 )//mod from iq's
{//the vec "l" has the segment lengths, rt is axis of bend
 float l2=(l.y+l.z)*sqrt(length(a-b))/sqrt(l.x+l.y+l.z);
 vec3 p=b-a,q=p*(0.5+0.5*(l.x*l.x-l2*l2)/dot(p,p));
 j1=a+q+sqrt(max(0.0,l.x*l.x-dot(q,q)))*normalize(cross(p,rt));
 p=b-j1;q=p*(0.5+0.5*(l.y*l.y-l.z*l.z)/dot(p,p));
 j2=j1+q+sqrt(max(0.0,l.y*l.y-dot(q,q)))*normalize(cross(p,rt));
}
*/

vec4 djsolve( vec2 p )//mod from fizzer's mod of iq's
{//the segment lengths are precalculated as 0.5,0.4,0.3 and starts at point 0,0
 float l2=0.639*pow(dot(p,p),0.25);
 vec2 q=p*(0.5+0.5*(0.25-l2*l2)/dot(p,p));
 vec2 j1=q+sqrt(max(0.0,0.25/dot(q,q)-1.0))*q.yx*vec2(-1.0,1.0);
 p=p-j1;
 q=p*(0.5+0.035/dot(p,p));
 vec2 j2=j1+q+sqrt(max(0.0,0.16/dot(q,q)-1.0))*q.yx*vec2(-1.0,1.0);
 return vec4(j1,j2);
}
mat3 lookat(vec3 fw,vec3 up){
 fw=normalize(fw);vec3 rt=normalize(cross(fw,up));return mat3(rt,cross(rt,fw),fw);
}
float Tube(vec3 pa, vec3 ba, float r){//mod from iq's
 float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
 return length(pa - ba*h)-r+h*0.015;
 //return length(pa-ba*clamp(dot(pa,ba)/dot(ba,ba),0.0,1.0))-r;
}

float finger(vec3 p, float z){
 vec3 j4=vec3(0.7+sin(z+time3+sin(z*1.3))*0.25,abs(sin(z+time3))*0.6,0.0);
 vec4 j=djsolve(j4.xy);
 vec3 j2=vec3(j.xy,0.0),j3=vec3(j.zw,0.0);
 float d=Tube(p,j2,0.145);
 d=min(d,Tube(p-j2,j3-j2,0.13));
 return min(d,Tube(p-j3,j4-j3,0.115));
}
float DE(vec3 p0){
 vec3 p=p0;
 p.x=mod(p.x,2.0)-1.0;
 float d3=length(max(abs(p.xy+vec2(0.0,0.27))-vec2(0.47,0.2),0.0))-0.03;
 float d2=min(length(p.xy-vec2(0.75,0.0)),length(p.xy+vec2(0.75,0.0)));
 float za=0.0;
 if(abs(fract(p.z)-0.5)>0.25){
  p.x=-p.x;
  p.z+=0.5;
  za=5.0;
 }
 p.x+=0.75;
 float z=floor(p.z)+floor(p0.x*0.5)*10.0+za;
 p.z=fract(p.z)-0.5;
 float d=finger(p,z);
 return min(0.16,min(d,min(d2-0.075,d3)));
}

float G(vec2 p){
  return smoothstep(0.0,0.05,min(max(abs(p.x-0.6)-0.15,abs(p.y-0.5)),max(min(p.x-p.y,p.y-0.5),abs(length(p-vec2(0.5))-0.25))));
  }

float A(vec2 p){
  return smoothstep(0.0,0.05,min(max(0.25-p.y,abs(abs(p.x-0.5)+p.y*0.5-0.4)),max(abs(p.x-0.5)-0.1,abs(p.y-0.5))));
  }

float C(vec2 p){
  return smoothstep(0.0,0.05,max(abs(length(p-vec2(0.5))-0.25),p.x-0.6));
  }
  
float T(vec2 p){
  return smoothstep(0.0,0.05,min(max(abs(p.x-0.5),abs(p.y-0.5)-0.25),max(abs(p.x-0.5)-0.25,abs(p.y-0.75))));
}

vec3 getColor(vec3 p0){
 float tim=time-0.2;
 float tim2=time2-0.2;
 float tim3=time3-0.2;
 vec3 col;
 vec3 p=p0;
 p.x=mod(p.x,2.0)-1.0;
 float d3=length(max(abs(p.xy+vec2(0.0,0.27))-vec2(0.47,0.2),0.0))-0.03;
 float d2=min(length(p.xy-vec2(0.75,0.0)),length(p.xy+vec2(0.75,0.0)));
 float ftz=0.1625,za=0.0;
 if(abs(fract(p.z)-0.5)>0.25){
  p.x=-p.x;
  p.z+=0.5;
  ftz=0.125;
  za=5.0;
 }
 p.x+=0.75;
 float z=floor(p.z)+floor(p0.x*0.5)*10.0+za;
 p.z=fract(p.z)-0.5;
 float d=finger(p,z);
 vec4 ft=vec4(0.7+sin(z+tim3+sin(z*1.3))*0.25,abs(sin(z+tim3))*0.6,p0.z+ftz,p.x);
 if(d2<d && d2<d3){col=vec3(0.7,0.65,0.4);}
 else if(d3<d){
  p0.z+=0.125;
  vec2 xz=abs(mod(p0.xz,0.25)-0.125);
  col=(1.0-20.0*dot(xz,xz))*vec3(0.5,0.24,0.3);
  xz=floor(p0.xz*4.0);
  if(ft.y<0.1){//fat fingering :)
   if(mod(floor(ft.w*4.0),4.0)==mod(floor(ft.x*4.0),4.0) && fract(ft.z)<0.25)
    col*=vec3(4.0,3.0,2.0)*(1.0+highhits);
  }
  float i=mod(xz.x-xz.y,4.0);
  xz=p0.xz*4.0-xz;
  if(i<1.0)col*=G(xz);
  else if(i<2.0)col*=A(xz);
  else if(i<3.0)col*=C(xz);
  else col*=T(xz);
 }else col=vec3(0.23,0.123,0.24);
 return col;
}
#ifdef SHADOWS
float linstep(float a, float b, float t){return clamp((t-a)/(b-a),0.,1.);}//from knighty
//random seed and generator
float randSeed,GoldenAngle;
float randStep(){//crappy random number generator
 randSeed=fract(randSeed+GoldenAngle);
 return  (0.8+0.2*randSeed);
}
float FuzzyShadow(vec3 ro, vec3 rd, float coneGrad){
 float t=0.01,d,s=1.0,r;
 ro+=rd*t;
 for(int i=0;i<8;i++){
  r=t*coneGrad;
  d=DE(ro+rd*t)+r*0.5;
  s*=linstep(-r,r,d);
  t+=d*randStep();
 }
 return clamp(s,0.2,1.0);
}
#endif

vec3 Light(vec3 p, vec3 rd, vec3 L){
 vec2 v=vec2(0.001,0.0);
 vec3 mcol=getColor(p);
 vec3 N=normalize(vec3(DE(p+v.xyy)-DE(p-v.xyy),DE(p+v.yxy)-DE(p-v.yxy),DE(p+v.yyx)-DE(p-v.yyx)));
 vec3 col=mix(mcol.yzx,mcol,abs(dot(rd,N)))*max(0.0,dot(N,L));
 col+=vec3(1.0,0.85,0.7)*pow(max(0.0,dot(reflect(rd,N),L)),26.0);
#ifdef SHADOWS
 col*=FuzzyShadow(p,L,0.5);
#endif
 return col;
}
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

 time=smoothTimeC*.75;
 time2=smoothTimeB*1.;
 time3=smoothTimeC*.75;

#ifdef SHADOWS
 GoldenAngle=2.0-0.5*(1.0+sqrt(5.0));
 randSeed=fract(sin(dot(fragCoord.xy,vec2(13.434,77.2378))+time*0.1)*4132.34526);
#endif
 vec3 ro=vec3(2.0*(1.0+sin(smoothTime*0.01)),1.5+direction,smoothTime*0.6+cos(smoothTime*0.025)*0.25);
 mat3 rotCam=lookat(vec3(0.0,0.3-0.5*ro.y,1.0),vec3(0.0,1.0,0.0));
 vec3 rd=rotCam*normalize(vec3(_uvc*PI*FOV+(2.0*fragCoord.xy-size.xy)/size.y,1.5));
 rd.x=-rd.x;
 float t=max(0.0,(0.875-ro.y)/rd.y),d=1.0,dm=d,tm=t;
 for(int i=0;i<24;i++){
  t+=d=DE(ro+rd*t);
  if(d<dm){dm=d;tm=t;}
 }
 vec3 col=vec3(rd.y*rd.y);
 vec3 L=normalize(vec3(0.2,0.6,-0.3));
 bool secondMarch=(dm>0.00012*tm);
 if(secondMarch){//unless we got really close speed up the march
  for(int i=0;i<16;i++){
   t+=d=DE(ro+rd*t)*1.5;//just experimenting
   t+=0.01;
  }
  if(d<0.02*t){
   vec3 p=ro+rd*t;
   vec3 scol=Light(p,rd,L);
   col=mix(scol,col,smoothstep(0.0,0.02*t,d));
   col*=exp(-t*0.1);
  }
 }
 if(dm<0.002*tm){
  vec3 p=ro+rd*tm;
  vec3 scol=Light(p,rd,L);
  col=mix(scol,col,smoothstep(0.0,0.002*tm,dm));
  col*=exp(-tm*0.1);
 }
 
 
 fragColor = vec4(clamp(col*2.0,0.0,1.0),1.0);

	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}