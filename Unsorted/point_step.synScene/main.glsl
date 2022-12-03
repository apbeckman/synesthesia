

			//******** BuffA Code Begins ********
#define A(u) texture(BuffA,(u)/RENDERSIZE.xy)
#define B(u) texture(BuffB,(u)/RENDERSIZE.xy)
#define C(u) texture(BuffC,(u)/RENDERSIZE.xy)

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;

    vec4  r = vec4(0);
    vec4  a = C(u);
    float z = 6.;//kernel convolution size
    for(float i=-z; i<=z; ++i){
    for(float j=-z; j<=z; ++j){
        float l = length(vec2(i,j));
        if((i==0.&&j==0.)|| l>z-.1){continue;}
        vec2  c = normalize(vec2(i,j));
        vec4  a2= C(u+vec2(i,j));
        vec4  b = a2-a;
        r.xy += c   *(sin(dot(b.xy,vec2(i,j))))*exp(-l*l*.2)*-.06;
        r.xy += a.xy*(sin(dot(b.xy,vec2(i,j))))*exp(-l*l*.3)*.45;
        r.xy += a.xy*(cos(dot(b.xy,vec2(i,j))))*exp(-l*l*.3)*-.1;
        
    }}
    fragColor = r;
	return fragColor; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;

    vec4 t = C(u);
    vec2 m = +t.xy +A(u).xy;
    float s = 0.;
    float z = 6.;//kernel convolution size
    for(float i=-z; i<=z; ++i){
    for(float j=-z; j<=z; ++j){
        if(length(vec2(i,j))>z-.1){continue;}
        vec2 c = (m+vec2(i,j));
        s += exp(-dot(c,c));
    }}
    if(s==0.){s = 1.;}
    s = 1./s;
    
    fragColor = vec4(m,s,0);
	return fragColor; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;

    float tz = 0.;
    vec4  a  = vec4(0);
    float z  = 6.;//kernel convolution size
    for(float i=-z; i<=z; ++i){
    for(float j=-z; j<=z; ++j){
        if(length(vec2(i,j))>z-.1){continue;}
        vec4 t = C(u+vec2(i,j));
        vec4 m = B(u+vec2(i,j));
        vec2 c = (m.xy-vec2(i,j));
        float z = t.z*exp(-dot(c,c));
        a.xy += z*m.xy;
        a.z  += z*m.z;
        tz   += z;
    }}
    if(tz==0.){tz = 1.;}
    a.xy /= tz;
    //a = A(u);
    if(_mouse.z>0.)
    {
        vec2 m = 22.*(u-_mouse.xy)/RENDERSIZE.y;
        a += vec4(-normalize(m),0,0)*exp(-dot(m,m));
    }
    //float keyW = texture( iChannel2, vec2(87.5/256.,.25) ).x;
    if(int(FRAMECOUNT)==1/*||Reset == 1.*/)
    {
        vec2 m = 8.*(u-RENDERSIZE.xy*.5)/450.;
        a = vec4(0,0,.5,0)-vec4(8,8,0,0)*exp(-dot(m,m));
    }
    fragColor = a;
	return fragColor; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 u = fragCoord/RENDERSIZE.xy;
    vec4 a = texture(BuffC,u);
    fragColor = a.zzzz;
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
		return renderMainImage();
	}
}