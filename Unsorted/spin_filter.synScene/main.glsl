

			//******** BuffA Code Begins ********

#define D(u) texture(BuffD,(u)/RENDERSIZE.xy)
#define A(u) texture(BuffA,(u)/RENDERSIZE.xy)
#define B(u) texture(BuffB,(u)/RENDERSIZE.xy)
#define C(u) texture(BuffC,(u)/RENDERSIZE.xy)
#define T(u) texture(image46,(u)/RENDERSIZE.xy)

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;

    vec4  r = vec4(0);
    vec4  a = D(u);
    float z    = 8.;//kernel convolution size
    float blur = 1./z;
    for(float i=-z; i<=z; ++i){
    for(float j=-z; j<=z; ++j){
        if((i==0.&&j==0.)||length(vec2(i,j))>z+.1){continue;}
        vec2  c = vec2(i,j)*blur;
        vec2  n = normalize(c);
        vec4  a2= D(u+vec2(i,j));
        vec4  b = a2-a;
        vec2  b2= b.yx*vec2(-1,1);
        r.xy += length(b.xy)*sin(dot(c,b2)*6.);
        r.z  += dot(n,b2);
    }}
    fragColor = r;
	return fragColor; 
 } 


			//******** BuffB Code Begins ********

#define D(u) texture(BuffD,(u)/RENDERSIZE.xy)
#define A(u) texture(BuffA,(u)/RENDERSIZE.xy)
vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;

    vec4  r = D(u);
    float z    = 8.;//kernel convolution size
    float blur = 1./z;
    for(float i=-z; i<=z; ++i){
    for(float j=-z; j<=z; ++j){
      if((i==0.&&j==0.)||length(vec2(i,j))>z+.1){continue;}
      vec4 b = A(u+vec2(i,j));
      vec2 c = (-vec2(i,j))*blur;
           c = normalize(c)*exp(-dot(c,c));//*r.z;
      r.xy += c*b.xy*b.z*.000001;
    }}
    fragColor = r;
	return fragColor; 
 } 


			//******** BuffC Code Begins ********

#define B(u) texture(BuffB,(u)/RENDERSIZE.xy)
#define A(u) texture(BuffA,(u)/RENDERSIZE.xy)
//#define C(u) texture(image46,(u)/RENDERSIZE.xy)
vec4 renderPassC() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;

    vec4 t = B(u);
    vec2 m = +t.xy
             +A(u).xy*(t.z-.5)*0.
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

#define B(u) texture(BuffB,(u)/RENDERSIZE.xy)
#define C(u) texture(BuffC,(u)/RENDERSIZE.xy)
vec4 renderPassD() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;

    float tz = 0.;
    vec4 a = vec4(0);
    float z    = 6.;//kernel convolution size
    for(float i=-z; i<=z; ++i){
    for(float j=-z; j<=z; ++j){
      vec4 t = A(u+vec2(i,j)); //t.z = 1.;
      vec4 m = C(u+vec2(i,j));
      vec2 c = (m.xy-vec2(i,j))*1.;
      float z = exp(-dot(c,c));
      a.xy += z*m.xy;
      a.z  += z*m.z*t.z;
      tz   += z;
    }}
    if(tz==0.){tz = 1.;}
    a.xy /= tz;
    //a = A(u);
    if(_mouse.z>0.)
    {
        vec2 m = 16.*(u-_mouse.xy)/RENDERSIZE.y;
        a += vec4(-1,0,0,0)*exp(-dot(m,m))*.1;
    }
    if(FRAMECOUNT==0)
    {
        vec2 m;
        m = 6.*(u-RENDERSIZE.xy*.5)/RENDERSIZE.y-vec2(2,1);
        a+= vec4(+2,0,0,0)*exp(-dot(m,m));
        m = 6.*(u-RENDERSIZE.xy*.5)/RENDERSIZE.y+vec2(2,1);
        a+= vec4(-2,0,0,0)*exp(-dot(m,m));
    }
    fragColor = a;
	return fragColor; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 u = fragCoord/RENDERSIZE.xy;
    vec4 a = texture(BuffD,u);
    fragColor = .5+vec4(0,a.xy,0)*.7;
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