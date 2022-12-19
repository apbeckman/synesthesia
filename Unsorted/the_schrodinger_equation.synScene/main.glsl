

			//******** Common Code Begins ********

#define size RENDERSIZE.xy
#define SAMPLR(a, pos) texture((a), (pos)/size)
#define SAMPLRs(a, pos,sz) texture((a), (pos)/sz)
#define A BuffA
#define B BuffB
#define C BuffC
#define D BuffD


//scales
#define dt sqrt(2.)/7.8
#define dx 1.
#define border 0.01
#define PI 3.14159265


vec4 Laplacian(sampler2D F, vec2 pos, vec2 sz)
{
    vec4 a = SAMPLRs(F, pos,sz);
    vec3 k = vec3(-1./12., 4./3., -5./2.);
	vec4 x_l1 = SAMPLRs(F, pos + vec2(-1,0),sz);
    vec4 x_l2 = SAMPLRs(F, pos + vec2(-2,0),sz);
    vec4 x_r1 = SAMPLRs(F, pos + vec2(1,0),sz);
    vec4 x_r2 = SAMPLRs(F, pos + vec2(2,0),sz);
    vec4 y_l1 = SAMPLRs(F, pos + vec2(0,-1),sz);
    vec4 y_l2 = SAMPLRs(F, pos + vec2(0,-2),sz);
    vec4 y_r1 = SAMPLRs(F, pos + vec2(0,1),sz);
    vec4 y_r2 = SAMPLRs(F, pos + vec2(0,2),sz);
 
	return ((x_l2+x_r2+y_l2+y_r2)*k.s + (x_l1+x_r1+y_l1+y_r1)*k.y + a*2.*k.z)/(dx*dx);
}


vec4 adam(sampler2D F, sampler2D Fp, vec2 pos, vec2 sz)
{
    return Laplacian(F, pos, sz);
}

float gaussian(vec2 pos, float r)
{
    pos /= r;
    return exp(-dot(pos,pos));
}

			//******** BuffA Code Begins ********

vec4 renderPassA() {
	vec4 field = vec4(0.0);
	vec2 pos = _xy;
   
    field = SAMPLR(D, pos);
    
    //mouse interaction
    if(_mouse.z>0.)
        field.x += gaussian(pos-_mouse.xy, 10.);
    
    field.x += dt*adam(D, B, pos,size).y;
    
    if(FRAMECOUNT <=1)
    {
        field = 5.*vec4(gaussian(pos-size*0.4, 10.)+gaussian(pos-size*0.6, 10.));
    }
	return field; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB() {
	vec4 field = vec4(0.0);
	vec2 pos = _xy;
   
    field = SAMPLR(A, pos);
    
    field.y += -dt*adam(A, C, pos,size).x;
      if(FRAMECOUNT <=1)
    {
        field = 5.*vec4(gaussian(pos-size*0.4, 10.)+gaussian(pos-size*0.6, 10.));
    }
	return field; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC() {
	vec4 field = vec4(0.0);
	vec2 pos = _xy;
   
    field = SAMPLR(B, pos);
    
    //mouse interaction
    if(_mouse.z>0.)
        field.x += gaussian(pos-_mouse.xy, 10.);
    
    field.x += dt*adam(B, D, pos,size).y;
      if(FRAMECOUNT <= 1)
    {
        field = 5.*vec4(gaussian(pos-size*0.4, 10.)+gaussian(pos-size*0.6, 10.));
    }
	return field; 
 } 


			//******** BuffD Code Begins ********

vec4 renderPassD() {
	vec4 field = vec4(0.0);
	vec2 pos = _xy;
   
    field = SAMPLR(C, pos);
    
    field.y += -dt*adam(C, A, pos,size).x;
    
      if(FRAMECOUNT <= 1)
    {
        field = 5.*vec4(gaussian(pos-size*0.4, 10.)+gaussian(pos-size*0.6, 10.));
    }
	return field; 
 } 



vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 pos = _xy;

    vec4  field = SAMPLR(C, pos);
  
    vec2 psi = field.xy; //probablility amplitude
  
     vec2 red = vec2(1,0), green = vec2(-0.5,sqrt(3.)*0.5), blue = vec2(-0.5,-sqrt(3.)*0.5);
   
    vec3 RGB =  vec3(dot(psi, red),dot(psi, green),dot(psi, blue));
    fragColor.xyz = RGB*RGB*2.;
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return renderPassC();
	}
	if(PASSINDEX == 3){
		return renderPassD();
	}
	if(PASSINDEX == 4){
		return renderMainImage();
	}
}