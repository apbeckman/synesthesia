vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


#define time TIME
#define size RENDERSIZE

float pixelSize,focalDistance,aperture,fudgeFactor=0.6;//,shadowCone=0.5;
mat2 rmat(float a){float sa=sin(a),ca=cos(a);return mat2(ca,sa,-sa,ca);}
vec3 mcol=vec3(0);
const float mr=0.16, mxr=1.0; 
const vec4 scale=vec4(-2.0,-2.0,-2.0,2.0); 
vec4 p0=vec4(3.,0.76,1.12,0.2);//0.32,.76
float lightPos; 
vec2 DE(in vec3 z0){//amazing surface by kali/tglad with mods
 p0.x=cos(((smoothTimeC*0.1)+z0.y)*0.05)*3.5;
 z0.xz=z0.xz*rmat(z0.y*0.07);
 z0.y=abs(mod(z0.y,4.0)-2.0);
 vec4 z = vec4(z0,1.0); float dL=100.;
 for (int n = 0; n < 4; n++) { 
  if(z.x<z.z)z.xz=z.zx; 
  z.xy=clamp(z.xy, -1.0, 1.0) *2.0-z.xy; 
  z*=scale/clamp(max(dot(z.xy,z.xy),dot(z.xz,z.xz)),mr,mxr); 
  z+=p0; 
  if(n==1)dL=length(z.xyz+vec3(0.5,lightPos,0.5))/z.w;
 } 
 if(mcol.x>0.)mcol+=vec3(0.6)+sin(z.xyz*0.1)*0.4; 
 z.xyz=abs(z.xyz)-vec3(1.4,32.8,0.7); 
 return vec2(max(z.x,max(z.y,z.z))/z.w,dL); 
} 

float CircleOfConfusion(float t){//calculates the radius of the circle of confusion at length t
 return max(abs(focalDistance-t)*aperture,pixelSize*(1.0+t));
}
mat3 lookat(vec3 fw,vec3 up){
 fw=normalize(fw);vec3 rt=normalize(cross(fw,normalize(up)));return mat3(rt,cross(rt,fw),fw);
}
float linstep(float a, float b, float t){return clamp((t-a)/(b-a),0.,1.);}// i got this from knighty and/or darkbeam
//random seed and generator
vec2 randv2;
float rand2(){// implementation derived from one found at: lumina.sourceforge.net/Tutorials/Noise.html
 randv2+=vec2(1.0,1.0);
 return fract(sin(dot(randv2 ,vec2(12.9898,78.233))) * 43758.5453);
}
vec3 bg(vec3 rd){
  float d=max(0.,rd.x+rd.y+rd.z);
  return vec3(d*d*.25)+rd*.05;
}
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

 randv2=fract(cos((fragCoord.xy+fragCoord.yx*vec2(100.0,100.0))+vec2(smoothTimeC*0.1)*10.0)*1000.0);
 pixelSize=1.0/size.y;
 float tim=smoothTime*0.1;//camera, lighting and object setup
 lightPos=sin(smoothTimeC*0.5)*5.5; 
 vec3 ro=vec3(cos(tim*0.5),tim*2.0,sin(tim*0.5))*5.0; 
 vec3 rd=lookat(vec3(-ro.x,5.0,-ro.z),vec3(0.0,1.0,1.0))*normalize(vec3((2.0*gl_FragCoord.xy-size.xy)/size.y,2.0)); 
 focalDistance=min(length(ro)+0.001,1.0);
 aperture=0.007*focalDistance;
 vec3 rt=normalize(cross(vec3(0,1,0),rd)),up=cross(rd,rt);//just need to be perpendicular
 vec3 lightColor=vec3(1.0,0.5,0.25)*2.0;
 vec4 col=vec4(0.0);vec3 blm=vec3(0);//color accumulator, .w=alpha, bloom accum
 vec2 D;//for surface and light dist
 float t=0.0,mld=100.0,od,d=1.,old,ld=100.,dt=0.,ot;//distance traveled, minimum light distance
 for(int i=1;i<70;i++){//march loop
  if(col.w>0.9 || t>15.0)break;//bail if we hit a surface or go out of bounds
  float rCoC=CircleOfConfusion(t);//calc the radius of CoC
  od=D.x;old=D.y,dt=t-ot;ot=t;//save old distances for normal, light direction calc
  D=DE(ro+rd*t);
  d=D.x+0.33*rCoC;
  ld=D.y;//the distance estimate to light
  mld=min(mld,ld);//the minimum light distance along the march
  if(d<rCoC){//if we are inside the sphere of confusion add its contribution
   vec3 p=ro+rd*(t-dt);//back up to previos checkpoint
   mcol=vec3(0.01);//collect color samples with normal deltas
   vec2 Drt=DE(p+rt*dt),Dup=DE(p+up*dt);
   vec3 N=normalize(rd*(D.x-od)+rt*(Drt.x-od)+up*(Dup.x-od));
   if(N!=N)N=-rd;//if no gradient assume facing us
   vec3 L=-normalize(rd*(D.y-old)+rt*(Drt.y-old)+up*(Dup.y-old));
   if(L!=L)L=up;
   float lightStrength=1.0/(1.0+ld*ld*20.0);
   vec3 scol=mcol*(0.4*(1.0+dot(N,L)+.2))*lightStrength;//average material color * diffuse lighting * attenuation
   scol+=pow(max(0.0,dot(reflect(rd,N),L)),8.0)*lightColor;//specular lighting
   mcol=vec3(0);//clear the color accumulator before shadows
   //scol*=FuzzyShadow(p,L,ld,shadowCone,rCoC);//now stop the shadow march at light distance
   blm+=lightColor*exp(-mld*t*10.)*(1.0-col.w);//add a bloom around the light
   mld=100.0;//clear the minimum light distance for the march
   float alpha=fudgeFactor*(1.0-col.w)*linstep(-rCoC,rCoC,-d);//calculate the mix like cloud density
   col=vec4(col.rgb+scol*alpha,clamp(col.w+alpha,0.0,1.0));//blend in the new color 
  }//move the minimum of the object and light distance
  d=abs(fudgeFactor*min(d,ld+0.33*rCoC)*(0.8+0.2*rand2()));//add in noise to reduce banding and create fuzz
  t+=d;
 }//mix in background color and remaining bloom
 t=min(15.,t);
 blm+=lightColor*exp(-mld*t*10.)*(1.0-col.w);///(1.0+mld*mld*3000.0
 col.rgb=mix(col.rgb,bg(rd),t/15.);
 fragColor = vec4(clamp(col.rgb+blm,0.0,1.0),1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}