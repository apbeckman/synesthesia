



			//******** BuffA Code Begins ********

float a = 1.;

float o = a*.7071;

vec2 xy_noise = vec2(_noise(vec2(sin(TIME*0.2), cos(TIME*0.2))));

const int samples = 45,
          LOD = 4,         // gaussian done on MIPmap at scale LOD
          sLOD = 2 << LOD; // tile size = 2^LOD
const float sigma = float(samples) * .25;

float gaussian(vec2 i) {
    return exp( -.5* dot(i/=sigma,i) ) / ( 6.28 * sigma*sigma );
}

vec4 blur(sampler2D sp, vec2 U, vec2 scale) {
    vec4 O = vec4(0);  
    int s = samples/sLOD;
    
    for ( int i = 0; i < s*s; i++ ) {
        vec2 d = vec2(i%s, i/s)*float(sLOD) - float(samples)/2.;
        O += gaussian(d) * textureLod( sp, U + scale * d , float(LOD) );
    }
    
    return O / O.a;
}

#define R RENDERSIZE.xy

vec4 D (vec2 U) {return texture(BuffD,U/R);}

vec4 renderPassA() {

	vec4 Q = vec4(0.0);

	vec2 U = _xy;
    
    // U.x*=(1.0-U.y*U.y)*U.x*0.2;
    // U.y=(1.0-U.x*U.x)*U.y*0.2;

    U += vec2(xy_noise)*0.15;

    U += vec2(dot(cos(_uvc.x),sin(_uvc.y)))*Pull*PI;

    U -= _uvc*Zoom*(1.0+2.0*low+syn_Intensity*1.5);

    U += _uvc*Stretch*(1.0+low);

    U += Drift*(1.0+low);

    Q = D(U);

    //float a = TIME+3.*length(U), c = cos(a), s = sin(a);

    float g = smoothTime+3.*length(U), c = cos(g), s = sin(g);



    mat2 m = mat2(c,-s,s,c);

    vec4 n = D(U+vec2(0,a)*m)+D(U+vec2(a,0)*m)+D(U+vec2(0,-a)*m)+D(U+vec2(-a,0)*m)+D(U+vec2(-o,o)*m)+D(U+vec2(o,-o)*m)+D(U+vec2(-o,-o)*m)+D(U+vec2(o,o)*m);

    n *= .125;

    vec4 dx = n-Q;

    Q += dx*vec4(1.,.3,1.,1);

    

    float x = (.3+complexity)*Q.x*Q.y*(1.-Q.y);

    

    Q.y = Q.y+x-0.025-.1*Q.z;

    Q.x = Q.x+0.11*Q.x*(1.-Q.x)-x;

    Q.z = Q.z*0.98+ .6*dx.y+.001*Q.y+Split*0.005;

    

    Q = clamp(Q,0.,1.);

    if (_mouse.z> 0. && length(U-_mouse.xy) < 40.) Q.y=1.;

    // if(_exists(syn_UserImage)) Q.y += _luminance(_edgeDetectSobel(syn_UserImage))*syn_Level*0.2;
    if(_exists(syn_UserImage)) {
        Q.y += mix(0., (_luminance(mix(_edgeDetectSobel(syn_UserImage), _loadMedia(), 0.35))), syn_Hits*0.025*syn_Intensity +syn_Presence*syn_Level*0.05);
    }

    if (FRAMECOUNT <= 1 || Reset != 0.) {

        

        Q = vec4(1,0,0,0);

        if (length(U-0.5*R)<20.) Q.y += 1.;

    }

	return Q; 

 } 





			//******** BuffB Code Begins ********



#define R RENDERSIZE.xy

vec4 A (vec2 U) {return texture(BuffA,U/R);}

vec4 renderPassB() {

	vec4 Q = vec4(0.0);

	vec2 U = _xy;

    U += vec2(dot(cos(_uvc.x),sin(_uvc.y)))*Pull*PI;


    Q = A(U);

    float g = smoothTime+3.*length(U), c = cos(g), s = sin(g);

    mat2 m = mat2(c,-s,s,c);

    o = a*.7071;

    vec4 n = A(U+vec2(0,a)*m)+A(U+vec2(a,0)*m)+A(U+vec2(0,-a)*m)+A(U+vec2(-a,0)*m)+A(U+vec2(-o,o)*m)+A(U+vec2(o,-o)*m)+A(U+vec2(-o,-o)*m)+A(U+vec2(o,o)*m);

    n *= .125;

    vec4 dx = n-Q;

    Q += dx*vec4(1.,.3,1.,1);

    

    float x = .53*Q.x*Q.y*(1.-Q.y);

    

    Q.y = Q.y+x-(0.025-Thiccness)-.1*Q.z;

    Q.x = Q.x+0.1*Q.x*(1.-Q.x)-x;

    Q.z = Q.z*0.98+ .6*dx.y+.001*Q.y;

    

    Q = clamp(Q,0.,1.);

    

    if (_mouse.z> 0. && length(U-_mouse.xy) < 40.) Q.y=1.;

    if(_exists(syn_UserImage)) {
        Q.y += mix(0., (_luminance(mix(_edgeDetectSobel(syn_UserImage), _loadMedia(), 0.5))), syn_Hits*0.05*syn_Intensity +syn_Presence*syn_Level*0.095);
    }


    if (FRAMECOUNT <= 1) {

        

        Q = vec4(1,0,0,0);

        if (length(U-0.5*R)<10.) Q.y += 1.;

    }

	return Q; 

 } 





			//******** BuffC Code Begins ********



#define R RENDERSIZE.xy

vec4 B (vec2 U) {return texture(BuffB,U/R);}

vec4 renderPassC() {

	vec4 Q = vec4(0.0);

	vec2 U = _xy;



    Q = B(U);

    float g = smoothTime+3.*length(U), c = cos(g), s = sin(g);

    mat2 m = mat2(c,-s,s,c);

    a += 1.0;

    o = a*.7071;

    vec4 n = B(U+vec2(0,a)*m)+B(U+vec2(a,0)*m)+B(U+vec2(0,-a)*m)+B(U+vec2(-a,0)*m)+B(U+vec2(-o,o)*m)+B(U+vec2(o,-o)*m)+B(U+vec2(-o,-o)*m)+B(U+vec2(o,o)*m);

    n *= .125;

    vec4 dx = n-Q;

    Q += dx*vec4(1.,.3,1.,1);

    

    float x = (.3+Test)*Q.x*Q.y*(1.-Q.y);

    

    Q.y = Q.y+x-(0.025)-.1*Q.z;

    Q.x = Q.x+0.1*Q.x*(1.-Q.x)-x;

    Q.z = Q.z*0.98+ .6*dx.y+.001*Q.y;

    

    Q = clamp(Q,0.,2.);

    

    if (_mouse.z> 0. && length(U-_mouse.xy) > 20. && length(U-_mouse.xy) < 40.) Q.y=1.;
    if(_exists(syn_UserImage)) {
        Q.y += mix(0., sin(_luminance(mix(_edgeDetectSobel(syn_UserImage), _loadMedia(), 0.5))), syn_Hits*0.05*syn_Intensity +syn_Presence*syn_Level*0.095);
    }

    

    if (FRAMECOUNT <= 1) {

        

        Q = vec4(1,0,0,0);

        if (length(U-0.5*R)<10.) Q.y += 1.;

    }

	return Q; 

 } 





			//******** BuffD Code Begins ********



vec4 mediaPass() {
    vec4 media = vec4(0.);
    vec2 U = _xy;
    media = _edgeDetectSobel(syn_UserImage);
    return media;
}
vec4 C (vec2 U) {return texture(BuffC,U/R);}

vec4 renderPassD() {

	vec4 Q = vec4(0.0);

	vec2 U = _xy;


    Q = C(U);

    float g = smoothTime+3.*length(U), c = cos(g), s = sin(g);

    mat2 m = mat2(c,-s,s,c);

    a += 1.0;

    o = a*.7071;

    vec4 n = C(U+vec2(0,a)*m)+C(U+vec2(a,0)*m)+C(U+vec2(0,-a)*m)+C(U+vec2(-a,0)*m)+C(U+vec2(-o,o)*m)+C(U+vec2(o,-o)*m)+C(U+vec2(-o,-o)*m)+C(U+vec2(o,o)*m);

    n *= .125;

    vec4 dx = n-Q;

    Q += dx*vec4(1.,.3,1.,1);

    

    float x = .3*Q.x*Q.y*(1.-Q.y);

    

    //Q.y = Q.y+x-0.025-.1*Q.z;

    Q.y = Q.y+x-(0.025-Thiccness)-.1*Q.z;

    Q.x = Q.x+0.1*Q.x*(1.-Q.x)-x;

    Q.z = Q.z*0.98+ .6*dx.y+.001*Q.y;

    

    Q = clamp(Q,0.,1.);

    

    if (_mouse.z> 0. && length(U-_mouse.xy) > 20. && length(U-_mouse.xy) < 40.) Q.y=1.;
    if(_exists(syn_UserImage)) {
        Q.y += mix(0., (_luminance(mix(_edgeDetectSobel(syn_UserImage), _loadMedia(), 0.5))), syn_Hits*0.05*syn_Intensity +syn_Presence*syn_Level*0.095);
    }

    

    if (FRAMECOUNT <= 1) {

        

        Q = vec4(1,0,0,0);

        if (length(U-0.55*R)<10.) Q.y += 1.;


    }

	return Q; 

 } 




float ln (vec3 p, vec3 q, vec3 b) {return length(p-q-(b-q)*dot(p-q,b-q)/dot(b-q,b-q));}

vec4 renderMainImage() {

	vec4 Q = vec4(0.0);

	vec2 U = _xy;


    vec4

        n = A(U+vec2(0,1)),

        e = A(U+vec2(1,0)),

        s = A(U-vec2(0,1)),

        w = A(U-vec2(1,0));

    vec2 h = vec2(e.w-w.w,n.w-s.w);

    vec2 g = vec2(e.y-w.y,n.y-s.y);

    vec4 b = A(U);

   	vec4 dx = 0.25*(n+e+s+w) - b;

	Q = (b.w+b.y)*normalize(abs(sin(6.+.75*(b.z-b.x)*vec4(3,2,1,4))));

    

    vec3 no = normalize(vec3(g-h,-1));

    vec3 r = reflect(vec3(0,0,1),no);

    vec3 l = vec3(2.*R.x,2.*R.y,10.*R.y);

    vec3 u = vec3(U,0);

    vec3 lu = l-u;

    float _o = length(r-lu*dot(r,lu)/dot(lu,lu));

    Q = (Q)*(exp(-o)+5.*exp(-17.5*o)+1e2*exp(-30.*o))*(1.0+syn_HighLevel*0.5);

    //Q *= 1.+texture(iChannel1,r);

    Q *= 1.+texture(prism,r.xz);

    no = normalize(vec3(g-h,1));

    r = refract(vec3(0,0,1),no,.2);

    l = vec3(3.*R.x,3.*R.y,10.*R.y);

    u = vec3(U,0);

    lu = l-u;

    o = length(r-lu*dot(r,lu)/dot(lu,lu));

    Q = 2.*dx+(Q*0.9+0.1)*(0.4+exp(-o)+10.*exp(-9.5*o)*(1.0+syn_HighLevel*0.5));

   // Q *= 1.+0.2*sin(Q+3.*texture(iChannel1,r).x*vec4(1,2,3,4));

    //Q *= 1.+0.3*sin(Q+3.*texture(prism,r.yz).x*vec4(1,2,3,4));
    Q *= 1.+0.3*sin(Q+3.*texture(prism,r.yz).x*vec4(1,2,3,4));

   //Q *= 1.+0.2*sin(Q+3.*vec4(1,2,3,4));

	return Q; 

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

		return mediaPass();

	}
	if(PASSINDEX == 5){

		return renderMainImage();

	}

}