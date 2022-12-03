

			//******** BuffA Code Begins ********
#define A(u) texture(BuffA,(u)/RENDERSIZE.xy)
#define B(u) texture(BuffB,(u)/RENDERSIZE.xy)
#define C(u) texture(BuffC,(u)/RENDERSIZE.xy)

vec2 hash( vec2 p ) // replace this by something better
{
	p = vec2( dot(p,vec2(127.1,311.7)), dot(p,vec2(269.5,183.3)) );
	return -1.0 + 2.0*fract(sin(p)*43758.5453123);
}
// Compact, self-contained version of IQ's 3D value noise function. I have a transparent noise
// example that explains it, if you require it.
float n3D(vec3 p){
    
	const vec3 s = vec3(7, 157, 113);
	vec3 ip = floor(p); p -= ip; 
    vec4 h = vec4(0., s.yz, s.y + s.z) + dot(ip, s);
    p = p*p*(3. - 2.*p); //p *= p*p*(p*(p * 6. - 15.) + 10.);
    h = mix(fract(sin(h)*43758.5453), fract(sin(h + s.x)*43758.5453), p.x);
    h.xy = mix(h.xz, h.yw, p.y);
    return mix(h.x, h.y, p.z); // Range: [0, 1].
}


float n2D( in vec2 p )
{
    const float K1 = 0.366025404; // (sqrt(3)-1)/2;
    const float K2 = 0.211324865; // (3-sqrt(3))/6;

	vec2  i = floor( p + (p.x+p.y)*K1 );
    vec2  a = p - i + (i.x+i.y)*K2;
    float m = step(a.y,a.x); 
    vec2  o = vec2(m,1.0-m);
    vec2  b = a - o + K2;
	vec2  c = a - 1.0 + 2.0*K2;
    vec3  h = max( 0.5-vec3(dot(a,a), dot(b,b), dot(c,c) ), 0.0 );
	vec3  n = h*h*h*h*vec3( dot(a,hash(i+0.0)), dot(b,hash(i+o)), dot(c,hash(i+1.0)));
    return dot( n, vec3(70.0) );
}


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
        r.xy += c*b.z*cos(dot(b.yx*vec2(-1,1),vec2(i,j))*.3)*exp(-l*l*.1)*.5;
        //r.xy += c*min(b.z,0.)*cos(dot(b.yx*vec2(-1,1),vec2(i,j))*.2)*exp(-l*l*.1)*.5;
    }}
    fragColor = r;
	return fragColor; 
 } 

float nT = n3D(vec3(cos(smoothTimeC*0.5), sin(smoothTimeC*0.25), cos(smoothTimeC*0.01)));
			//******** BuffB Code Begins ********

//#define C(u) texture(BuffC,(u)/RENDERSIZE.xy)
//#define B(u) texture(BuffA,(u)/RENDERSIZE.xy)
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

//#define C(u) texture(BuffC,(u)/RENDERSIZE.xy)
//#define B(u) texture(BuffB,(u)/RENDERSIZE.xy)
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
    //float keyW = texture( _mouse.zzzz, vec2(87.5/256.,.25) ).x;
    if(int(FRAMECOUNT)==1||Reset == 1.)
    {
        vec2 m = 22.*(u-RENDERSIZE.xy*.5)/450.;
        a = vec4(0,0,.5,0)-vec4(8,8,0,0)*exp(-dot(m,m));
    }
    fragColor = a;
	return fragColor; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 u = fragCoord/RENDERSIZE.xy;
    u += _uvc*PI;
    vec4 a = texture(BuffC,u);
    fragColor = (+sin(a.x+vec4(0,1,2,4)+5.)*.5
                 +sin(a.y+vec4(5,3,1,4)+5.)*.5+1.)*a.z;
    //fragColor = a.zzzz;
    //fragColor = a;
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