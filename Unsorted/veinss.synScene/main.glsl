

			//******** BuffA Code Begins ********

#define A(u) texture(BuffA,(u)/RENDERSIZE.xy)
#define B(u) texture(BuffB,(u)/RENDERSIZE.xy)
#define C(u) texture(BuffC,(u)/RENDERSIZE.xy)
#define D(u) texture(BuffD,(u)/RENDERSIZE.xy)
#define D(u) texture(BuffD,(u)/RENDERSIZE.xy)
#define T(u) texture(image46,(u)/RENDERSIZE.xy)
float growthFactor = normalize(pow((syn_BassLevel*0.5)+(syn_MidLevel*0.35)+(syn_Level*0.15), 2.0));

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;

    vec4  r = vec4(0);
    vec4  a = D(u);
    float z = 6.*K;
    float t = 0.;
    for(float i=-z; i<=z; ++i){
    for(float j=-z; j<=z; ++j){
        vec2  ij= vec2(i,j);
        vec2  c = ij*(3./z);
        float l = length(ij);
        ij /= l; if(l==0.){ij= vec2(0);}
        float e = exp(-dot(c,c)); t+=e;
        vec4 b = D(u+vec2(i,j))-a;
        r.x += e*dot(b.xy,ij);
        r.y += abs(e*dot(b.xy,ij.yx*vec2(-1,1)));
    }}
    fragColor = r/t;
	return fragColor; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;
    u += PI*MoveXY;

    vec4  r = D(u);
    vec4  a = A(u);
    a.xy += (_uvc*PI*Flow);

    float z = 6.;
    for(float i=-z; i<=z; ++i){
    for(float j=-z; j<=z; ++j){
        vec2  ij= vec2(i,j);
        vec2  c = ij*(3./z);
        float l = length(ij);
        ij /= l; if(l==0.){ij= vec2(0);}
        float e = exp(-dot(c,c));
        vec4 b = A(u+vec2(i,j));
        r.xy += +(b.x)    *ij*e*.1
                +(b.y-a.y)*ij*e*.5;
    }}

    fragColor = r;
	return fragColor; 
 } 


			//******** BuffC Code Begins ********
vec4 renderPassC() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;
    u -= PI*_uvc*Zoom;

    vec4 t = B(u);
    vec2 m = +t.xy
             +B(u).xy*(t.z-.5)*0.
             +t.z*vec2(0,.0)
             -T(u).x*t.xy*.0;
    float s = 0.;
    float z    = 6.;//kernel convolution size
    for(float i=-z; i<=z; ++i){
    for(float j=-z; j<=z; ++j){
      vec2 c = (m+vec2(i,j))*1.;
      s += exp(-dot(c,c));
    }}
    if(s==0.){s = 1.;}
    s = 1./s;
    
    fragColor = vec4(m,s,0);
	return fragColor; 
 } 


			//******** BuffD Code Begins ********

vec4 renderPassD() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;

    vec4 a = vec4(0);
    float z    = 6.;//kernel convolution size
    for(float i=-z; i<=z; ++i){
    for(float j=-z; j<=z; ++j){
      vec4 t = C(u+vec2(i,j)); t.z = 1.;
      vec4 m = C(u+vec2(i,j));
      vec2 c = (m.xy-vec2(i,j))*1.;
      float z = exp(-dot(c,c))*m.z*t.z;
      a.xy += z*m.xy;
      a.z  += z;
    }}
    float tz = 1./a.z; if(a.z==0.){tz = 0.;}
    a.xy *= tz;
    //a = A(u);
    if(_mouse.z>0.)
    {
        vec2 m = 42.*(u-_mouse.xy)/RENDERSIZE.y;
        a += vec4(m,0,0)*exp(-dot(m,m))*-.2;
    }
    if(FRAMECOUNT<=1 || Reset > 0.)
    {
        vec2 m = (42.*Size)*(u-RENDERSIZE.xy*.5)/RENDERSIZE.y;
        a = vec4(m,1,1)*exp(-dot(m,m));
    }
    fragColor = a;
	return fragColor; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 u = fragCoord/RENDERSIZE.xy;
    vec4 a = texture(BuffD,u);
    fragColor =+sin(a.x*4.+vec4(1,3,5,4))*.25
               +sin(a.y*4.+vec4(1,3,2,4))*.25
               +.5;
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